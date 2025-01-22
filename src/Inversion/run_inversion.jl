#run inversion 

#function run_inversion!(simulation::Simulation, surf_speeds::Array{Float64,1},accumulation_rate::Array{Float64,2},dhdt::Array{Float64,2})

function run_inversion!(simulation,inversion)
   
    @unpack model,timestepping_params, output_params, clock = simulation
   
    # get_neumann_velocs
    update_velocities!(model)
    update_velocities_on_h_grid!(model)

    #get_dirichlet_velocs
    solve_dirichelt_velocities!(model,inversion)

    #Will then need to update lots of things here - see matlab
    # get_heating_rates
    # update_beta_inverse
    # update_preBfactor_inverse
    # update viscosity and F0/1/2 here
    # get_JKV

    #   update_surface_elevation!(model)
    #   update_geometry_on_uv_grids!(model)
    #   update_height_above_floatation!(model)
    #   update_grounded_fraction_on_huv_grids!(model)
    #   update_accumulation_rate!(model)
    #   update_basal_melt!(model, clock)
    #   update_weertman_c!(model)
    #   update_dsdh!(model)
    #   update_model_velocities!(model)
    #   update_velocities_on_h_grid!(model)
    #   update_dhdt!(model)
    #   update_model_wavelets!(model)
        
    return model
end

#function solve_dirichelt_velocities!(model::AbstractModel{T,N}, surf_speeds::Array{Float64,1}, accumulation_rate::Array{Float64,2},dhdt::Array{Float64,2}) where {T,N}
function solve_dirichelt_velocities!(model, inversion) 

    #    function solve_dirichelt_velocities!(model::Model, dhdt::Array{Float64,2}) 
    @unpack params,solver_params=model
    @unpack gu,gv,wu,wv = model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 

    x0=get_start_guess_dirichlet(model,inversion)

    #get the ops needed:
    Op_combo=get_op_combo(model,inversion)
  #  println("Size of Op_combo is," ,size(Op_combo),length(Op_combo))

    #get rhs f1, and get rhsf2:
    #using f1 same as in forward problem: this is b1 in Arthern at al 2015 (A3) 
    #NEED TO CHECK SCHUR HERE!!!
    f1=get_rhs_dirichlet_inversion(model)
  #  println("Size of f1 (rhs) is," ,size(f1))

    #get f2
#    f2=get_rhs_dirichlet_inverse_data(model,surf_speeds,accumulation_rate,dhdt)
    f2=get_rhs_dirichlet_inverse_data(model,inversion)
   # println("Size of f2 (rhs) is," ,size(f2))

    f=[f1;f2]

    println("Size of Op_combo is " ,size(Op_combo))
    println("size of x0 is " ,size(x0))
    println( "and size of f is ",size(f))

    resid=get_resid(x0,Op_combo,f)
  
    println("Initial residual norm: ", norm(resid))
    relative_residual = norm(resid) / norm(f)
    println("Initial relative Residual: ", relative_residual)

     ni=gu.ni+gv.ni+gudata.ni + gvdata.ni + ghdata.n
     xm=Vector{Float64}(undef,ni);
    read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_x0guess.bin",xm)
    xm.=ntoh.(xm)
    nirhsd1=gu.ni+gv.ni
    rhsd1=Vector{Float64}(undef,nirhsd1);
    read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_RHSD1.bin",rhsd1)
    rhsd1.=ntoh.(rhsd1)
    nirhsd2=gudata.ni + gvdata.ni + ghdata.n
    rhsd2=Vector{Float64}(undef,nirhsd2);
    read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_RHSD2.bin",rhsd2)
    rhsd2.=ntoh.(rhsd2)
    residm=get_resid(xm,Op_combo,f)
    println("Matlab residual norm: ", norm(residm))
    relative_residualm = norm(residm) / norm(f)
    println("Matlab relative Residual: ", relative_residualm)

  #  println("f  element: ",f2[7000:7005])
  #  println("f  element: ",f2[19000:19005])

    fm=[rhsd1;rhsd2]
    f_diff=(f1-rhsd1)

    index = argmax(f_diff) 

    println("max f diff is" ,f_diff[index])
    println("mean f diff is" ,sum(f_diff)/length(f_diff))
    println("rhsd1 at max is " ,rhsd1[index])
    println("f1 at max is " ,f1[index])

    ff=gg
    
    reltol = 1e-6  # Tolerance for convergence
    max_iter = 60000  # Maximum number of iterations

    x = gmres!(x0,Op_combo, f, abstol=reltol, maxiter=max_iter, restart=400, log=true,verbose=true )
    x_output=x[1]
    println("Check convergence:"  ,x[2])
    
   # x_output = cg!(x0, Op_combo, f,maxiter=max_iter)
    
    resid=get_resid(x_output,Op_combo,f)
    println("Residual norm: ", norm(resid))

    relative_residual = norm(resid) / norm(f)
    println("Relative Residual: ", relative_residual)

    ff=gg

   # u_d_p_d = cg!(x0, Op_combo, f)
   # u_d = u_d_p_d[1:n]
   # p_d = u_d_p_d[n+1:end] 

    #    cg!(z2, op, z1; PL=p_dirichlet)

    #    update_preconditioner!(model)
    
    println("Solved for Dirichlet velocites")
    
    return model
end
   
"""
    get_start_guess_dirichlet(model::AbstractModel)

 Return starting guess used to begin iterative solution of Dirichlet problem.

"""
function get_start_guess_dirichlet(model,inversion)
    #still needs editing!!
    @unpack gu,gv,gh=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
    @assert eltype(gu.u)==eltype(gv.v)
   #    x=[gu.samp_inner*gu.u[:];gv.samp_inner*gv.v[:]]
    x=[gu.samp_inner*gu.u[:];gv.samp_inner*gv.v[:]; ones(gudata.ni);ones(gvdata.ni); zeros(ghdata.n)]
#  ni=gu.ni+gv.ni+gudata.ni + gvdata.ni + ghdata.n
#  x=Vector{Float64}(undef,ni);
#read!("/data/hpcdata/users/chll1/WAVI_Initial_Data_github/WAVI-WAIS-setups/inversion_data/bedmachinev3/full_stripe_fix_8km/Inverse_x0guess.bin",x)
#x.=ntoh.(x)

    return x
end

function get_op_combo(model::AbstractModel{T,N}, inversion) where {T,N}
    @unpack gu,gv,gh=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
    ni = gu.ni + gv.ni+gudata.ni + gvdata.ni+ghdata.n
    op_fun_combo! = get_op_fun_combo(model,inversion)
    op_combo=LinearMap{T}(op_fun_combo!,ni;issymmetric=true,ismutating=true,ishermitian=true,isposdef=true)
end

function get_op_fun_combo(model::AbstractModel{T,N}, inversion) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 

        
    #Preallocate intermediate variables used by op_fun_combo
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    uvsampi :: Vector{T} = zeros(gu.ni+gv.ni);                               @assert length(uvsampi) == gu.ni+gv.ni
    pisampi :: Vector{T} = zeros(gudata.ni+gvdata.ni+ghdata.n);              @assert length(pisampi) == gudata.ni+gvdata.ni+ghdata.n
    f1 :: Vector{T} = zeros(gu.ni+gv.ni);                                    @assert length(f1) == gu.ni+gv.ni
    f2 :: Vector{T} = zeros(gudata.ni+gvdata.ni+ghdata.n);                   @assert length(f2) == gudata.ni+gvdata.ni+ghdata.n
    opvecprod :: Vector{T} = zeros(gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n);         @assert length(opvecprod) == gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n

    function op_fun_combo!(opvecprod::AbstractArray,inputVector::AbstractVector;vecSampled::Bool=true)
        @assert length(inputVector)==(gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n)

         #Split vector into u- and v-components for tau and h- for sigma
         uvsampi .= @view inputVector[1:(gu.ni+gv.ni)]
         pisampi .= @view inputVector[(gu.ni+gv.ni+1):(gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n)]

            #Contruct A:
        op_A=get_op(model)
     #   println("Size of op (A) is," ,size(op_A))
     
        #Constuct B:
        op_B=get_op_B(model,inversion)
     #   println("Size of op_B (B) is," ,size(op_B))
      #  println("Type of op_B is " ,typeof(op_B))

        #Constuct B-transpose:
        op_BT=get_op_BT(model,inversion)
     #   println("Size of op_B (BT) is," ,size(op_BT))
     #   println("Type of op_BT is " ,typeof(op_BT))

        # Construct C:
        op_C=get_op_C(model,inversion)
      #  println("Size of op_C is," ,size(op_C))
      #  println("Type of op_C is " ,typeof(op_C))

     #   println("Size of uvsampi is," ,size(uvsampi))
     #   println("Size of pisampi is," ,size(pisampi))
 
        f1 = op_A*uvsampi + op_BT*pisampi 
      #  println("Size of f1 is," ,size(f1))
     #   println("Type of f1 is " ,typeof(f1))
       # println(typeof(f1))
        f2 = op_B*uvsampi + op_C*pisampi
      #  println("Size of f2 is," ,size(f2))
      #  println("Type of f2 is " ,typeof(f2))

        opvecprod[1:gu.ni+gv.ni] .= f1
        opvecprod[(gu.ni+gv.ni+1):(gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n)] .= f2

        return opvecprod
    end
    #Return op_fun_combo as a closure
    return op_fun_combo!
end


function get_op_A(model::AbstractModel{T,N}) where {T,N}
    @unpack gu,gv,gh=model.fields
    n_input = gu.ni + gv.ni
    n_output=gu.ni + gv.ni + gh.n
    op_fun_A! = get_op_fun_A(model)
    op_A=LinearMap{T}(
        op_fun_A!,           # The mutating function for A
        n_output,                  # Output size
        n_input,                 # Input size
        issymmetric=false,    # B^T is not symmetric in general
        ismutating=true,      # Function modifies input/output
        ishermitian=false,    # B^T is not Hermitian
        isposdef=false        # B^T is not positive definite
    )
end


function get_op_fun_A(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    
    #Preallocate intermediate variables used by op_fun
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    usampi :: Vector{T} = zeros(gu.ni);                             @assert length(usampi) == gu.ni
    vsampi :: Vector{T} = zeros(gv.ni);                             @assert length(vsampi) == gv.ni
    uspread :: Vector{T} = zeros(nxnyu);                            @assert length(uspread) == nxnyu
    vspread :: Vector{T} = zeros(nxnyv);                            @assert length(vspread) == nxnyv
    dudx :: Vector{T} = zeros(nxnyh);                               @assert length(dudx) == nxnyh
    dvdy :: Vector{T} = zeros(nxnyh);                               @assert length(dvdy) == nxnyh
    r_xx_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xx_strain_rate_sum) == nxnyh
    r_yy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_yy_strain_rate_sum) == nxnyh
    r_xx :: Vector{T} = zeros(nxnyh);                               @assert length(r_xx) == nxnyh
    r_yy :: Vector{T} = zeros(nxnyh);                               @assert length(r_yy) == nxnyh
    dudy_c :: Vector{T} = zeros(nxnyc);                             @assert length(dudy_c) == nxnyc
    dvdx_c :: Vector{T} = zeros(nxnyc);                             @assert length(dvdx_c) == nxnyc
    r_xy_strain_rate_sum_c :: Vector{T} = zeros(nxnyc);             @assert length(r_xy_strain_rate_sum_c) == nxnyc
    r_xy_strain_rate_sum_crop_c :: Vector{T} = zeros(nxnyc);        @assert length(r_xy_strain_rate_sum_crop_c) == nxnyc
    r_xy_c :: Vector{T} = zeros(nxnyc);                             @assert length(r_xy_c) == nxnyc
    r_xy_crop_c :: Vector{T} = zeros(nxnyc);                        @assert length(r_xy_crop_c) == nxnyc
    d_rxx_dx :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxx_dx) == nxnyu
    d_rxy_dy :: Vector{T} = zeros(nxnyu);                           @assert length(d_rxy_dy) == nxnyu
    d_ryy_dy :: Vector{T} = zeros(nxnyv);                           @assert length(d_ryy_dy) == nxnyv
    d_rxy_dx :: Vector{T} = zeros(nxnyv);                           @assert length(d_rxy_dx) == nxnyv
    taubx :: Vector{T} = zeros(nxnyu);                              @assert length(taubx) == nxnyu
    tauby :: Vector{T} = zeros(nxnyv);                              @assert length(tauby) == nxnyv
    qx :: Vector{T} = zeros(nxnyu);                                 @assert length(qx) == nxnyu
    qx_crop :: Vector{T} = zeros(nxnyu);                            @assert length(qx_crop) == nxnyu
    dqxdx :: Vector{T} = zeros(nxnyh);                              @assert length(dqxdx) == nxnyh
    qy :: Vector{T} = zeros(nxnyv);                                 @assert length(qy) == nxnyv
    qy_crop :: Vector{T} = zeros(nxnyv);                            @assert length(qy_crop) == nxnyv
    dqydy :: Vector{T} = zeros(nxnyh);                              @assert length(dqydy) == nxnyh
    divq :: Vector{T} = zeros(nxnyh);                               @assert length(divq) == nxnyh
    extra :: Vector{T} = zeros(nxnyh);                              @assert length(extra) == nxnyh
    d_extra_dx :: Vector{T} = zeros(nxnyu);                         @assert length(d_extra_dx) == nxnyu
    d_extra_dy :: Vector{T} = zeros(nxnyv);                         @assert length(d_extra_dy) == nxnyv
    h_d_extra_dx :: Vector{T} = zeros(nxnyu);                       @assert length(h_d_extra_dx) == nxnyu
    h_d_extra_dy :: Vector{T} = zeros(nxnyv);                       @assert length(h_d_extra_dy) == nxnyv
    fx :: Vector{T} = zeros(nxnyu);                                 @assert length(fx) == nxnyu
    fy :: Vector{T} = zeros(nxnyv);                                 @assert length(fy) == nxnyv
    fx_sampi :: Vector{T} = zeros(gu.ni);                           @assert length(fx_sampi) == gu.ni
    fy_sampi :: Vector{T} = zeros(gv.ni);                           @assert length(fy_sampi) == gv.ni
    opvecprod :: Vector{T} = zeros(gu.ni+gv.ni);                    @assert length(opvecprod) == gu.ni + gv.ni

    function op_fun_A!(opvecprod::AbstractVector,inputVector::AbstractVector;vecSampled::Bool=true)
        if vecSampled
            @assert length(inputVector)==(gu.ni+gv.ni)

            #Split vector into u- and v- components
            usampi .= @view inputVector[1:gu.ni]
            vsampi .= @view inputVector[(gu.ni+1):(gu.ni+gv.ni)]

            #Spread to vectors that include all grid points within rectangular domain.
            @!  uspread = gu.spread_inner*usampi
            @!  vspread = gv.spread_inner*vsampi

        else
            #Vector already includes all grid points within rectangular domain.
            @assert length(inputVector)==(gu.nxu*gu.nyu+gv.nxv*gv.nyv)
            
            uspread .= @view inputVector[1:gu.nxu*gu.nyu]
            vspread .= @view inputVector[(gu.nxu*gu.nyu+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv)]

        end
            #Extensional resistive stresses
        @!  dudx = gu.∂x*uspread
        @!  dvdy = gv.∂y*vspread
        @.  r_xx_strain_rate_sum = 2dudx + dvdy
        @.  r_yy_strain_rate_sum = 2dvdy + dudx
        @!  r_xx = gh.dneghηav[]*r_xx_strain_rate_sum
        @.  r_xx = -2r_xx
        @!  r_yy = gh.dneghηav[]*r_yy_strain_rate_sum
        @.  r_yy = -2r_yy

            #Shearing resistive stresses
        @!  dudy_c = gu.∂y*uspread
        @!  dvdx_c = gv.∂x*vspread
        @.  r_xy_strain_rate_sum_c = dudy_c + dvdx_c
        @!  r_xy_strain_rate_sum_crop_c = gc.crop*r_xy_strain_rate_sum_c
        @!  r_xy_c = gc.dneghηav[]*r_xy_strain_rate_sum_crop_c
        @.  r_xy_c = -r_xy_c
        @!  r_xy_crop_c = gc.crop*r_xy_c

            #Gradients of resisitve stresses
        @!  d_rxx_dx = gu.∂xᵀ*r_xx
        @.  d_rxx_dx = - d_rxx_dx 
        @!  d_rxy_dy = gu.∂yᵀ*r_xy_crop_c
        @.  d_rxy_dy = - d_rxy_dy 
        @!  d_ryy_dy = gv.∂yᵀ*r_yy
        @.  d_ryy_dy = -d_ryy_dy
        @!  d_rxy_dx = gv.∂xᵀ*r_xy_crop_c
        @.  d_rxy_dx = -d_rxy_dx

            #Basal drag
        @!  taubx = gu.dnegβeff[]*uspread
        @.  taubx = -taubx
        @!  tauby = gv.dnegβeff[]*vspread
        @.  tauby = -tauby
            
            #Extra terms arising from Schur complement of semi implicit system (Arthern et al. 2015).
            qx .= vec(gu.h).*uspread
        @!  qx_crop = gu.crop*qx
        @!  dqxdx = gu.∂x*qx_crop
            qy .=  vec(gv.h).* vspread
        @!  qy_crop = gv.crop*qy
        @!  dqydy = gv.∂y*qy_crop
            divq .= dqxdx .+ dqydy
        @!  extra = gh.dimplicit[]*divq
        @!  d_extra_dx = gu.∂xᵀ*extra
        @.  d_extra_dx = -d_extra_dx
        @!  d_extra_dy = gv.∂yᵀ*extra
        @.  d_extra_dy = -d_extra_dy
            h_d_extra_dx .= vec(gu.h).*d_extra_dx
            h_d_extra_dy .= vec(gv.h).*d_extra_dy

            #Resistive forces resolved in x anf y directions
            fx .= d_rxx_dx .+ d_rxy_dy .- taubx
            # .- h_d_extra_dx
            fy .= d_ryy_dy .+ d_rxy_dx .- tauby 
            #.- h_d_extra_dy

            #Resistive forces sampled at valid grid points
        @!  fx_sampi = gu.samp_inner*fx
        @!  fy_sampi = gv.samp_inner*fy

            opvecprod[1:gu.ni] .= fx_sampi
            opvecprod[(gu.ni+1):(gu.ni+gv.ni)] .= fy_sampi

            return opvecprod
    end

    #Return op_fun as a closure
    return op_fun_A!
end



function get_op_B(model::AbstractModel{T,N},inversion) where {T,N}
        @unpack gu,gv,gh=model.fields
        @unpack gudata, gvdata, ghdata=inversion.data_fields 
        n_input = gu.ni + gv.ni
        n_output=gudata.ni + gvdata.ni+ghdata.n
        op_fun_B! = get_op_fun_B(model,inversion)
        op_B=LinearMap{T}(
            op_fun_B!,           # The mutating function for B^T
            n_output,                  # Output size
            n_input,                 # Input size
            issymmetric=false,    # B^T is not symmetric in general
            ismutating=true,      # Function modifies input/output
            ishermitian=false,    # B^T is not Hermitian
            isposdef=false        # B^T is not positive definite
        )
end



function get_op_fun_B(model::AbstractModel{T,N},inversion) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
        
    #Preallocate intermediate variables used by op_fun
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    usampi :: Vector{T} = zeros(gu.ni);                             @assert length(usampi) == gu.ni
    vsampi :: Vector{T} = zeros(gv.ni);                             @assert length(vsampi) == gv.ni
    uspread :: Vector{T} = zeros(nxnyu);                            @assert length(uspread) == nxnyu
    vspread :: Vector{T} = zeros(nxnyv);                            @assert length(vspread) == nxnyv
    R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
    R_crop :: Vector{T} = zeros(nxnyh);                             @assert length(R_crop) == nxnyh
    u_on_h :: Vector{T} = zeros(nxnyh);                             @assert length(u_on_h) == nxnyh
    Ru_on_h :: Vector{T} = zeros(nxnyh);                            @assert length(Ru_on_h) == nxnyh
    Ru_on_u :: Vector{T} = zeros(nxnyu);                            @assert length(Ru_on_u) == nxnyu
 #   Ru_on_u_i:: Vector{T} = zeros(gu.ni);                           @assert length(Ru_on_u_i) == gu.ni
    v_on_h :: Vector{T} = zeros(nxnyh);                             @assert length(v_on_h) == nxnyh
    Rv_on_h :: Vector{T} = zeros(nxnyh);                            @assert length(Rv_on_h) == nxnyh
    Rv_on_v :: Vector{T} = zeros(nxnyv);                            @assert length(Rv_on_v) == nxnyv
 #   Rv_on_v_i:: Vector{T} = zeros(gv.ni);                           @assert length(Rv_on_v_i) == gv.ni
    fx1 :: Vector{T} = zeros(nxnyu);                                @assert length(fx1) == nxnyu
    fy1 :: Vector{T} = zeros(nxnyv);                                @assert length(fy1) == nxnyv
    fx1_sampi :: Vector{T} = zeros(gudata.ni);                      @assert length(fx1_sampi) == gudata.ni
    fy1_sampi :: Vector{T} = zeros(gvdata.ni);                      @assert length(fy1_sampi) == gvdata.ni
    hu :: Vector{T} = zeros(nxnyu);                                 @assert length(hu) == nxnyu
    hu_crop :: Vector{T} = zeros(nxnyu);                            @assert length(hu_crop) == nxnyu
    divhu :: Vector{T} = zeros(nxnyh);                              @assert length(divhu) == nxnyh
    divhu_crop :: Vector{T} = zeros(nxnyh);                         @assert length(divhu_crop) == nxnyh
    hv :: Vector{T} = zeros(nxnyv);                                 @assert length(hv) == nxnyv
    hv_crop :: Vector{T} = zeros(nxnyv);                            @assert length(hv_crop) == nxnyv
    divhv :: Vector{T} = zeros(nxnyh);                              @assert length(divhv) == nxnyh
    divhv_crop :: Vector{T} = zeros(nxnyh);                         @assert length(divhv_crop) == nxnyh
  #  div_huv :: Vector{T} = zeros(nxnyh);                            @assert length(div_huv) == nxnyh 
  #  div_huv_samp :: Vector{T} = zeros(gh.n);                        @assert length(div_huv_samp) == gh.n 
    fx :: Vector{T} = zeros(nxnyh);                                 @assert length(fx) == nxnyh
    fx_sampi :: Vector{T} = zeros(ghdata.n);                            @assert length(fx_sampi) == ghdata.n
    opvecprod :: Vector{T} = zeros(gudata.ni+gvdata.ni+ghdata.n);               @assert length(opvecprod) == gudata.ni + gvdata.ni + ghdata.n
#    opvecprod :: Vector{T} = zeros(udi + vdi + hdi);                @assert length(opvecprod) == udi + vdi + hdi

    function op_fun_B!(opvecprod::AbstractVector,inputVector::AbstractVector;vecSampled::Bool=true)
       #vecSampled=false
       if vecSampled
        @assert length(inputVector)==(gu.ni+gv.ni)

        #Split vector into u- and v- components
        usampi .= @view inputVector[1:gu.ni]
        vsampi .= @view inputVector[(gu.ni+1):(gu.ni+gv.ni)]

        #Spread to vectors that include all grid points within rectangular domain.
        @!  uspread = gu.spread_inner*usampi
        @!  vspread = gv.spread_inner*vsampi

        else
        #Vector already includes all grid points within rectangular domain.
        @assert length(inputVector)==(gu.nxu*gu.nyu+gv.nxv*gv.nyv)
        uspread .= @view inputVector[1:gu.nxu*gu.nyu]
        vspread .= @view inputVector[(gu.nxu*gu.nyu+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv)]
        end
    

   #     R=S_u* S_u^T (c_h^u)^T R (c_h^u)
        
    #    move u onto h grid using gu.cent*uspread
    #   multiply R by this R*gu.cent*uspread
    #   then move back onto u grid gu.cent^T R*gu.cent*uspread
    #   then select only points in u mask. 
    #   Then select only points with u* data there using new mask.

        #Compute R from Arthern et al 2015 and build the R term from equation (28):
        @. R = (1.0 + gh.β[:]*gh.quad_f1[:])/ (1.0 + gh.β[:]*gh.quad_f2[:])
        @! R_crop=gh.crop*R

        @! u_on_h=gu.cent*(uspread)
        @. Ru_on_h=R_crop*u_on_h
        @! Ru_on_u=gu.centᵀ*Ru_on_h
      #  @! Ru_on_u_i=gu.samp_inner*Ru_on_u
        @! fx1=gu.crop*Ru_on_u

        @! v_on_h=gv.cent*(vspread)
        @. Rv_on_h=R_crop*v_on_h
        @! Rv_on_v=gv.centᵀ*Rv_on_h
    #    @! Rv_on_v=gv.samp_inner*Rv_on_v
        @! fy1=gv.crop*Rv_on_v
       
        #sample these to only include data masked points:
        @! fx1_sampi = gudata.samp_inner*fx1
        @! fy1_sampi = gvdata.samp_inner*fy1



        #Then look at the div term in the B matrix - this will multiply the velocties:
        #-div(hu)

        @.  hu=gu.h[:]*uspread
        @!  hu_crop=gu.crop*hu
        @!  divhu=gu.∂x*hu_crop
        @!  divhu_crop=gh.crop*divhu

        @.  hv=gv.h[:]*vspread
        @!  hv_crop=gv.crop*hv
        @!  divhv=gv.∂y*hv_crop
        @!  divhv_crop=gh.crop*divhv

        @! fx = divhu_crop .+ divhv_crop
        #sample these to only include data masked points:
        @! fx_sampi = ghdata.samp*fx

        opvecprod[1:gudata.ni] .= fx1_sampi
        opvecprod[(gudata.ni+1):(gudata.ni+gvdata.ni)] .= fy1_sampi
        opvecprod[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)] .= fx_sampi

        return opvecprod
    end

    #Return op_fun as a closure
    return op_fun_B!
end


function get_op_BT(model::AbstractModel{T,N},inversion) where {T,N}
        @unpack gu,gv,gh=model.fields
        @unpack gudata, gvdata, ghdata=inversion.data_fields 
     #   n_input = udi + vdi + hdi
        n_input = gudata.ni + gvdata.ni+ghdata.n
        n_output=gu.ni + gv.ni
        op_fun_BT! = get_op_fun_BT(model,inversion)
        op_BT = LinearMap{T}(
            op_fun_BT!,           # The mutating function for B^T
            n_output,                  # Output size
            n_input;                 # Input size
            issymmetric=false,    # B^T is not symmetric in general
            ismutating=true,      # Function modifies input/output
            ishermitian=false,    # B^T is not Hermitian
            isposdef=false        # B^T is not positive definite
        )
#        op_BT = LinearMap{T}((input, output) -> op_fun_BT!(output, input), n_u, n_pi;
 #                        issymmetric=false, ismutating=true, ishermitian=false, isposdef=false)
        #    op_BT=LinearMap{T}(op_fun_BT!,ni;issymmetric=true,ismutating=true,ishermitian=true,isposdef=true)
end


function get_op_fun_BT(model::AbstractModel{T,N}, inversion) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
        
    #Preallocate intermediate variables used by op_fun
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    tausampui :: Vector{T} = zeros(gudata.ni);                          @assert length(tausampui) == gudata.ni
    tausampvi :: Vector{T} = zeros(gvdata.ni);                          @assert length(tausampvi) == gvdata.ni
    sigmassampi :: Vector{T} = zeros(ghdata.n);                        @assert length(sigmassampi) == ghdata.n
    tausspreadu :: Vector{T} = zeros(nxnyu);                          @assert length(tausspreadu) == nxnyu
    tausspreadv :: Vector{T} = zeros(nxnyv);                          @assert length(tausspreadv) == nxnyv
    sigmasspread :: Vector{T} = zeros(nxnyh);                         @assert length(sigmasspread) == nxnyh
    R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
  #  RT :: Array{T,2} = zeros(gh.nyh,gh.nxh);                        @assert length(RT) == gh.nxh*gh.nyh
    tauspreadu_on_h :: Vector{T} = zeros(nxnyh);                      @assert length(tauspreadu_on_h) == nxnyh
    tauspreadu_c :: Vector{T} = zeros(nxnyu);                         @assert length(tauspreadu_c) == nxnyu
    tauspreadu_c_on_h :: Vector{T} = zeros(nxnyh);                      @assert length(tauspreadu_c_on_h) == nxnyh
    Rtauu_on_h :: Vector{T} = zeros(nxnyh);                           @assert length(Rtauu_on_h) == nxnyh
    Rtauu_on_u :: Vector{T} = zeros(nxnyu);                           @assert length(Rtauu_on_u) == nxnyu
   # RTtauu_on_u_crop :: Vector{T} = zeros(nxnyu);                   @assert length(RTtauu_on_u_crop) == nxnyu
    tauspreadv_on_h :: Vector{T} = zeros(nxnyh);                        @assert length(tauspreadv_on_h) == nxnyh
    tauspreadv_c :: Vector{T} = zeros(nxnyv);                         @assert length(tauspreadv_c) == nxnyv
    tauspreadv_c_on_h :: Vector{T} = zeros(nxnyh);                      @assert length(tauspreadv_c_on_h) == nxnyh
    Rtauv_on_h :: Vector{T} = zeros(nxnyh);                           @assert length(Rtauv_on_h) == nxnyh
    Rtauv_on_v :: Vector{T} = zeros(nxnyv);                           @assert length(Rtauv_on_v) == nxnyv 
   # RTtauv_on_v_crop :: Vector{T} = zeros(nxnyv);                   @assert length(RTtauv_on_v_crop) == nxnyv
    #
    sigmasspread_crop :: Vector{T} = zeros(nxnyh);                    @assert length(sigmasspread_crop) == nxnyh 
    dsigmadx :: Vector{T} = zeros(nxnyu);                           @assert length(dsigmadx) == nxnyu
    dsigmady :: Vector{T} = zeros(nxnyv);                           @assert length(dsigmady) == nxnyv
    dsigmadx_crop :: Vector{T} = zeros(nxnyu);                      @assert length(dsigmadx_crop) == nxnyu
    dsigmady_crop :: Vector{T} = zeros(nxnyv);                      @assert length(dsigmady_crop) == nxnyv
    h_dsigmadx :: Vector{T} = zeros(nxnyu);                         @assert length(h_dsigmadx) == nxnyu
    h_dsigmady :: Vector{T} = zeros(nxnyv);                         @assert length(h_dsigmady) == nxnyv
    h_dsigmadx_crop :: Vector{T} = zeros(nxnyu);                    @assert length(h_dsigmadx_crop) == nxnyu
    h_dsigmady_crop :: Vector{T} = zeros(nxnyv);                    @assert length(h_dsigmady_crop) == nxnyv
    fx1 :: Vector{T} = zeros(nxnyu);                                 @assert length(fx1) == nxnyu
    fy1 :: Vector{T} = zeros(nxnyv);                                 @assert length(fy1) == nxnyv
    fx2 :: Vector{T} = zeros(nxnyu);                                 @assert length(fx2) == nxnyu
    fy2 :: Vector{T} = zeros(nxnyv);                                 @assert length(fy2) == nxnyv
    fx1_sampi :: Vector{T} = zeros(gu.ni);                           @assert length(fx1_sampi) == gu.ni
    fy1_sampi :: Vector{T} = zeros(gv.ni);                           @assert length(fy1_sampi) == gv.ni
    fx2_sampi :: Vector{T} = zeros(gu.ni);                           @assert length(fx2_sampi) == gu.ni
    fy2_sampi :: Vector{T} = zeros(gv.ni);                           @assert length(fy2_sampi) == gv.ni
    fx_sampi :: Vector{T} = zeros(gu.ni);                           @assert length(fx_sampi) == gu.ni
    fy_sampi :: Vector{T} = zeros(gv.ni);                           @assert length(fy_sampi) == gv.ni
    opvecprod :: Vector{T} = zeros(gu.ni+gv.ni);                    @assert length(opvecprod) == gu.ni + gv.ni

    function op_fun_BT!(opvecprod::AbstractVector,inputVector::AbstractVector;vecSampled::Bool=true)
       #vecSampled=false
        if vecSampled

            @assert length(inputVector)==(gudata.ni+gvdata.ni+ghdata.n)

            #Split vector into u- and v-components for tau and h- for sigma
            tausampui .= @view inputVector[1:gudata.ni]
            tausampvi .= @view inputVector[(gudata.ni+1):(gudata.ni+gvdata.ni)]
            sigmassampi .= @view inputVector[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)]

            #Spread to vectors that include all grid points within rectangular domain.
            @!  tausspreadu = gudata.spread_inner*tausampui
            @!  tausspreadv = gvdata.spread_inner*tausampvi
            @!  sigmasspread =  ghdata.spread*sigmassampi

            else
            #Vector already includes all grid points within rectangular domain.
            @assert length(inputVector)==(gu.nxu*gu.nyu+gv.nxv*gv.nyv+gh.nxh*gh.nyh)

            tausspreadu .= @view inputVector[1:gu.nxu*gu.nyu]
            tausspreadv .= @view inputVector[(gu.nxu*gu.nyu+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv)]
            sigmasspread .= @view inputVector[(gu.nxu*gu.nyu+gv.nxv*gv.nyv+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv+gh.nxh*gh.nyh)]
    
        end
    
        #Compute R^T from Arthern et al 2015 and build the R^T term from equation (28):
        @. R = (1.0 + gh.β[:]*gh.quad_f1[:])/ (1.0 + gh.β[:]*gh.quad_f2[:])
       # R_mat=(1.0 .+ gh.β.*gh.quad_f1)./ (1.0 .+ gh.β.*gh.quad_f2)
     #   RT = transpose(R_mat)

##########SHOULD THESE CROPS BE DATA.CROPS????
        @!  tauspreadu_c=gu.crop*(tausspreadu)
        @!  tauspreadu_on_h=gu.cent*(tauspreadu_c)
        @!  tauspreadu_c_on_h=gh.crop*(tauspreadu_on_h)
        @.  Rtauu_on_h=R*tauspreadu_c_on_h
        @!  Rtauu_on_u=gu.centᵀ*Rtauu_on_h

        @!  tauspreadv_c=gv.crop*(tausspreadv)
        @!  tauspreadv_on_h=gv.cent*(tauspreadv_c)
        @!  tauspreadv_c_on_h=gh.crop*(tauspreadv_on_h)
        @.  Rtauv_on_h=R*tauspreadv_c_on_h
        @!  Rtauv_on_v=gv.centᵀ*Rtauv_on_h
  
        fx1 = Rtauu_on_u
        fy1 = Rtauv_on_v
  
        #sample these to only include masked points:
        @!  fx1_sampi = gu.samp_inner*fx1
        @!  fy1_sampi = gv.samp_inner*fy1
        
#=       #  tauspreadu_on_h=gu.cent*(gu.crop*(tausspreadu))
        RTtauu_on_h=RT[:].*tauspreadu_on_h
        RTtauu_on_u=gu.centᵀ*RTtauu_on_h
    #    tauspreadv_on_h=gv.cent*(gv.crop*(tausspreadv))
        RTtauv_on_h=RT[:].*tauspreadv_on_h
        RTtauv_on_v=gv.centᵀ*RTtauv_on_h
         fx1_sampi_n = gu.samp_inner*RTtauu_on_u
         fy1_sampi_n = gv.samp_inner*RTtauv_on_v
         diff_if=fx1_sampi-fx1_sampi_n
         println("diff on L606 is " ,count(!iszero ,diff_if)) =#


        #Then calculate grad term in the B matrix: hgrad(sigma_s)

        @!  sigmasspread_crop=gh.crop*sigmasspread
        @!  dsigmadx=gu.∂xᵀ*sigmasspread_crop
        @!  dsigmadx_crop=gu.crop*dsigmadx
        @!  h_dsigmadx=gu.h[:].*dsigmadx_crop
        @!  h_dsigmadx_crop=gu.crop*h_dsigmadx
  
        @!  dsigmady=gv.∂yᵀ*sigmasspread_crop
        @!  dsigmady_crop=gv.crop*dsigmady
        @!  h_dsigmady=gv.h[:].*dsigmady_crop
        @!  h_dsigmady_crop=gv.crop*h_dsigmady
        
        fx2 =  h_dsigmadx_crop
        fy2 =  h_dsigmady_crop

        #sample these to only include masked points:
        @!  fx2_sampi = gu.samp_inner*fx2
        @!  fy2_sampi = gv.samp_inner*fy2

        #Putting together both parts of BT operator using addition:
        fx_sampi.= fx1_sampi.+ fx2_sampi
        fy_sampi.= fy1_sampi.+ fy2_sampi

        opvecprod[1:gu.ni] .= fx_sampi
        opvecprod[(gu.ni+1):(gu.ni+gv.ni)] .= fy_sampi

        return opvecprod
    end

    #Return op_fun as a closure
    return op_fun_BT!
end

function get_op_C(model::AbstractModel{T,N},inversion) where {T,N}
    @unpack gu,gv,gh=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
    ni = gudata.ni + gvdata.ni + ghdata.n
    op_fun_C! = get_op_fun_C(model,inversion)
    op_C=LinearMap{T}(op_fun_C!,ni;issymmetric=true,ismutating=true,ishermitian=true,isposdef=true)
end

function get_op_fun_C(model::AbstractModel{T,N},inversion) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
        
    #Preallocate intermediate variables used by op_fun_C
    nxnyh :: N = gh.nxh*gh.nyh
    nxnyu :: N = gu.nxu*gu.nyu
    nxnyv :: N = gv.nxv*gv.nyv
    nxnyc :: N = gc.nxc*gc.nyc
    tausampui :: Vector{T} = zeros(gudata.ni);                          @assert length(tausampui) == gudata.ni
    tausampvi :: Vector{T} = zeros(gvdata.ni);                          @assert length(tausampvi) == gvdata.ni
    sigmassampi :: Vector{T} = zeros(ghdata.n);                        @assert length(sigmassampi) == ghdata.n
    tausspreadu :: Vector{T} = zeros(nxnyu);                          @assert length(tausspreadu) == nxnyu
    tausspreadv :: Vector{T} = zeros(nxnyv);                          @assert length(tausspreadv) == nxnyv
    tausspreadu_crop :: Vector{T} = zeros(nxnyu);                     @assert length(tausspreadu_crop) == nxnyu
    tausspreadv_crop :: Vector{T} = zeros(nxnyv);                     @assert length(tausspreadv_crop) == nxnyv
    sigmasspread :: Vector{T} = zeros(nxnyh);                         @assert length(sigmasspread) == nxnyh
    sigmasspread_crop :: Vector{T} = zeros(nxnyh);                    @assert length(sigmasspread_crop) == nxnyh 
    R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
    P :: Vector{T} = zeros(nxnyh);                                  @assert length(P) == nxnyh
    tauspreadu_on_h :: Vector{T} = zeros(nxnyh);                      @assert length(tauspreadu_on_h) == nxnyh
    Ptauu_on_h :: Vector{T} = zeros(nxnyh);                           @assert length(Ptauu_on_h) == nxnyh
    Ptauu_on_u :: Vector{T} = zeros(nxnyu);                           @assert length(Ptauu_on_u) == nxnyu
    Ptauu_on_h_crop :: Vector{T} = zeros(nxnyh);                      @assert length(Ptauu_on_h_crop) == nxnyh
    Ptauu_on_u_crop :: Vector{T} = zeros(nxnyu);                      @assert length(Ptauu_on_u_crop) == nxnyu
    tauspreadv_on_h :: Vector{T} = zeros(nxnyh);                      @assert length(tauspreadv_on_h) == nxnyh
    Ptauv_on_h :: Vector{T} = zeros(nxnyh);                           @assert length(Ptauv_on_h) == nxnyh
    Ptauv_on_v :: Vector{T} = zeros(nxnyv);                           @assert length(Ptauv_on_v) == nxnyv 
    Ptauv_on_h_crop :: Vector{T} = zeros(nxnyh);                      @assert length(Ptauv_on_h_crop) == nxnyh
    Ptauv_on_v_crop :: Vector{T} = zeros(nxnyv);                      @assert length(Ptauv_on_v_crop) == nxnyv 
   # Psigma :: Vector{T} = zeros(nxnyh);                             @assert length(Psigma) == nxnyh
    fx1 :: Vector{T} = zeros(nxnyu);                                 @assert length(fx1) == nxnyu
    fy1 :: Vector{T} = zeros(nxnyv);                                 @assert length(fy1) == nxnyv
    fh :: Vector{T} = zeros(nxnyh);                                  @assert length(fh)== nxnyh
    fx1_sampi :: Vector{T} = zeros(gudata.ni);                       @assert length(fx1_sampi) == gudata.ni
    fy1_sampi :: Vector{T} = zeros(gvdata.ni);                       @assert length(fy1_sampi) == gvdata.ni
    fh_sampi :: Vector{T} = zeros(ghdata.n);                         @assert length(fh_sampi) == ghdata.n
    opvecprod :: Vector{T} = zeros(gudata.ni+gvdata.ni+ghdata.n);    @assert length(opvecprod) == gudata.ni + gvdata.ni + ghdata.n
   
    function op_fun_C!(opvecprod::AbstractArray,inputVector::AbstractVector;vecSampled::Bool=true)
     
        if vecSampled

            @assert length(inputVector)==(gudata.ni+gvdata.ni+ghdata.n)

            #Split vector into u- and v-components for tau and h- for sigma
            tausampui .= @view inputVector[1:gudata.ni]
            tausampvi .= @view inputVector[(gudata.ni+1):(gudata.ni+gvdata.ni)]
            sigmassampi .= @view inputVector[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)]

            #Spread to vectors that include all grid points within rectangular domain.
            @!  tausspreadu = gudata.spread_inner*tausampui
            @!  tausspreadv = gvdata.spread_inner*tausampvi
            @!  sigmasspread =  ghdata.spread*sigmassampi

            else
            #Vector already includes all grid points within rectangular domain.
            @assert length(inputVector)==(gu.nxu*gu.nyu+gv.nxv*gv.nyv+gh.nxh*gh.nyh)

            tausspreadu .= @view inputVector[1:gu.nxu*gu.nyu]
            tausspreadv .= @view inputVector[(gu.nxu*gu.nyu+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv)]
            sigmasspread .= @view inputVector[(gu.nxu*gu.nyu+gv.nxv*gv.nyv+1):(gu.nxu*gu.nyu+gv.nxv*gv.nyv+gh.nxh*gh.nyh)]
    
        end

            #Compute R from Arthern et al 2015 and build the R term from equation (28):
        @.  R = (1.0 + gh.β[:]*gh.quad_f1[:])/ (1.0 + gh.β[:]*gh.quad_f2[:])
        
        #THIS STILL NEEDS CHANGING TO F0!!!!!
       # P=gh.quad_f1[:].-(R.+1).*gh.quad_f1[:].+R.*gh.quad_f2[:]
        @.  P=gh.quad_f0[:]-(R+1)*gh.quad_f1[:]+R*gh.quad_f2[:]
        

        @!  tausspreadu_crop=gu.crop*tausspreadu
        @!  tauspreadu_on_h=gu.cent*tausspreadu_crop
        @.  Ptauu_on_h=P[:]*tauspreadu_on_h
        @!  Ptauu_on_h_crop=gh.crop*Ptauu_on_h
        @!  Ptauu_on_u=gu.centᵀ*Ptauu_on_h_crop
        @!  Ptauu_on_u_crop=gu.crop*Ptauu_on_u

        @!  tausspreadv_crop=gv.crop*tausspreadv
        @!  tauspreadv_on_h=gv.cent*tausspreadv_crop
        @.  Ptauv_on_h=P[:]*tauspreadv_on_h
        @!  Ptauv_on_h_crop=gh.crop*Ptauv_on_h
        @!  Ptauv_on_v=gv.centᵀ*Ptauv_on_h_crop
        @!  Ptauv_on_v_crop=gv.crop*Ptauv_on_v
        #still need to add sampler for selecting only points where there is data

        fx1 = Ptauu_on_u_crop
        fy1 = Ptauv_on_v_crop

        #sample these to only include masked points:
        @!  fx1_sampi = gudata.samp_inner*fx1
        @!  fy1_sampi = gvdata.samp_inner*fy1

        fh_sampi=zeros(ghdata.n)

        opvecprod[1:gudata.ni] .= fx1_sampi
        opvecprod[(gudata.ni+1):(gudata.ni+gvdata.ni)] .= fy1_sampi
        opvecprod[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)] .= fh_sampi
        return opvecprod
    end
    #Return op_fun_C as a closure
    return op_fun_C!
end



function get_rhs_dirichlet_inverse_data(model, inversion)
    @unpack gh,gu,gv,gc=model.fields
    @unpack params, solver_params = model
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
       
     #Preallocate intermediate variables used 
     nxnyh = gh.nxh*gh.nyh
     nxnyu = gu.nxu*gu.nyu
     nxnyv = gv.nxv*gv.nyv
     
     us_data= zeros(nxnyv);                               
     vs_data= zeros(nxnyv);                               
     us_data_sampi  = zeros(gudata.ni);                      
     vs_data_sampi  = zeros(gvdata.ni);                    
     dhdt_data = zeros(nxnyh);                              
     accumulation_data = zeros(nxnyh);                     
     dhdtacc_data= zeros(nxnyh);                         
     dhdtacc_data_sampi= zeros(ghdata.n);                 
     f1 = zeros(gudata.ni+gvdata.ni+ghdata.n);            
     rhs_dirichlet = zeros(gudata.ni+gvdata.ni+ghdata.n);
    
     us_data = gu.centᵀ*(gu.cent*(gudata.speed_u[:]))
     vs_data =  gv.centᵀ*(gv.cent*(gvdata.speed_v[:]))

     us_data_sampi = gudata.samp_inner*us_data
     vs_data_sampi =gvdata.samp_inner*vs_data
          
     f1[1:gudata.ni] .= us_data_sampi
     f1[(gudata.ni+1):(gudata.ni+gvdata.ni)] .= vs_data_sampi
    
     dhdt_data = ghdata.dhdt
     accumulations_data = ghdata.accumulation_rate
      dhdtacc_data=dhdt_data[:] + accumulations_data[:]

     dhdtacc_data_sampi=ghdata.samp*dhdtacc_data
    
    f1[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)] .=  dhdtacc_data_sampi
        
    rhs_dirichlet .= f1 
    
    return rhs_dirichlet
end
    

"""
    get_rhs_dirichlet_inversion(model::AbstractModel)

 Return right hand side vector of momentum equations.

"""
function get_rhs_dirichlet_inversion(model::AbstractModel{T,N}) where {T,N}
    @unpack gh,gu,gv,gc=model.fields
    @unpack params, solver_params = model
    onesvec=ones(T,gh.nxh*gh.nyh)
    surf_crop= gh.crop*(gh.s[:])
    
    rhs = zeros(T,gu.ni+gv.ni)
    f1 = zeros(T,gu.ni+gv.ni)
    f2 = zeros(T,gu.ni+gv.ni)
    f3 = zeros(T,gu.ni+gv.ni)
    tmph = zeros(T,gh.nxh*gh.nyh)
    tmpu = zeros(T,gu.nxu*gu.nyu)
    tmpui = zeros(T,gu.ni)
    tmpv = zeros(T,gv.nxv*gv.nyv)
    tmpvi = zeros(T,gv.ni)
    sui = zeros(T,gu.ni)
    hui = zeros(T,gu.ni)
    dui = zeros(T,gu.ni)
    svi = zeros(T,gv.ni)
    hvi = zeros(T,gv.ni)
    dvi = zeros(T,gv.ni)


@!  tmpu = gu.∂xᵀ*surf_crop
@.  tmpu = -tmpu
@!  tmpui = gu.samp_inner*tmpu
@.  tmpui = (params.density_ice*params.g*gu.h[gu.mask_inner]).* tmpui

@!  tmpv = gv.∂yᵀ*surf_crop
@.  tmpv = -tmpv
@!  tmpvi = gv.samp_inner*tmpv
@.  tmpvi = (params.density_ice*params.g*gv.h[gv.mask_inner]).*tmpvi
              
    f1[1:gu.ni] .= tmpui
    f1[(gu.ni+1):(gu.ni+gv.ni)] .= tmpvi

    sui .= gu.s[gu.mask_inner]
    hui .= gu.h[gu.mask_inner]
    dui .= icedraft.(sui,hui,params.sea_level_wrt_geoid)
@!  tmph = gh.crop*onesvec
@!  tmpu = gu.∂xᵀ*tmph
@.  tmpu = -tmpu
@!  tmpui = gu.samp_inner*tmpu
@.  tmpui = tmpui * params.g*(0.5*params.density_ice*hui^2
                            - 0.5*params.density_ocean*dui^2
                            - params.density_ice*hui*sui)

    svi .= gv.s[gv.mask_inner]
    hvi .= gv.h[gv.mask_inner]
    dvi .= icedraft.(svi,hvi,params.sea_level_wrt_geoid)
@!  tmph = gh.crop*onesvec
@!  tmpv = gv.∂yᵀ*tmph
@.  tmpv = -tmpv
@!  tmpvi = gv.samp_inner*tmpv
@.  tmpvi = tmpvi * params.g*(0.5*params.density_ice*hvi^2
                            - 0.5*params.density_ocean*dvi^2
                            - params.density_ice*hvi*svi)

    f2[1:gu.ni] .= tmpui
    f2[(gu.ni+1):(gu.ni+gv.ni)] .= tmpvi

 #   f2=[
 #       (0.5*params.density_ice*params.g*gu.h[gu.mask_inner].^2
 #       .- 0.5*params.density_ocean*params.g*(icedraft.(gu.s[gu.mask_inner],gu.h[gu.mask_inner],params.sea_level_wrt_geoid)).^2
 #       .- params.density_ice*params.g*gu.h[gu.mask_inner].*gu.s[gu.mask_inner]).*gu.samp_inner*(-gu.∂xᵀ*(gh.crop*onesvec))
 #       ;
 #       (0.5*params.density_ice*params.g*gv.h[gv.mask_inner].^2
 #       .- 0.5*params.density_ocean*params.g*(icedraft.(gv.s[gv.mask_inner],gv.h[gv.mask_inner],params.sea_level_wrt_geoid)).^2
 #       .- params.density_ice*params.g*gv.h[gv.mask_inner].*gv.s[gv.mask_inner]).*gv.samp_inner*(-gv.∂yᵀ*(gh.crop*onesvec))
 #       ]

        get_rhs_dirichlet!(f3,model)

        rhs .= f1 .+ f2 .+ f3

    return rhs
end
