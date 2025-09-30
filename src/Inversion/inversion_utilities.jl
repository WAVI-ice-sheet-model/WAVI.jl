#run inversion 

"""
     update_betas_dirichlet!(model,inversion)

update beta and beta_eff for the dirichlet solution by setting as equal to the Neumann (model) fields
"""

function update_betas_dirichlet!(model,inversion)
    inversion.fields.gh.β[inversion.fields.gh.mask] .= model.fields.gh.β[model.fields.gh.mask]
    inversion.fields.gh.βeff[inversion.fields.gh.mask]  =model.fields.gh.βeff[model.fields.gh.mask]
    inversion.fields.gu.βeff[inversion.fields.gu.mask] =model.fields.gu.βeff[model.fields.gu.mask]
    inversion.fields.gv.βeff[inversion.fields.gv.mask] =model.fields.gv.βeff[model.fields.gv.mask]
 return inversion
end

"""
    start_guess_β_inversion!(model::AbstractModel)

set the drag coefficient at the start of the inversion.
"""
function start_guess_β_inversion!(model::AbstractModel, inversion)
    @unpack gh=model.fields
    @unpack inversion_params=inversion
    aground = (gh.haf .>= 0)
    gh.β .= inversion_params.βfloating_start*ones(gh.nxh,gh.nyh)
    gh.β[aground].=inversion_params.βgrounded_start
    return model
end

"""
    start_guess_η_inversion!(model::AbstractModel)

set  η at the start of the inversion.
"""
function start_guess_η_inversion!(model::AbstractModel, inversion)
    @unpack gh,g3d=model.fields
    @unpack inversion_params=inversion
    for k=1:g3d.nσs
        for j=1:g3d.nys
            for i=1:g3d.nxs
                if gh.mask[i,j]
                        g3d.η[i,j,k] = inversion_params.ηstart_guess
                end
            end
        end
    end
    return model
end


function solve_dirichelt_neumann_velocities!(model, inversion,clock) 

    @unpack params,solver_params=model
    @unpack gu,gv,wu,wv,gh = model.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
    @unpack inversion_params=inversion

        #get rhs f1, and get rhsf2:
        #using f1 same as in forward problem: this is b1 in Arthern at al 2015 (A3) 
        f1=get_rhs_dirichlet_inversion(model)

        #Contruct A:
        op_A=get_op_A(model)
        #Constuct B:
        op_B=get_op_B(model,inversion)
        #Constuct B-transpose:
        op_BT=get_op_BT(model,inversion)
        op_C=get_op_C(model,inversion)

        #using inbuilt function op_diag:
        A_diag_vals=get_op_diag(model,op_A)
        M1 = Diagonal(A_diag_vals)
        M2=[]

        #Neumann solve:
        uN=get_start_guess(model)
        
        uN, ch1 = cg!(uN, op_A, f1,  abstol=inversion_params.abstol, maxiter=inversion_params.maxiter, log=true, verbose=false, Pl=M1)
        if !@isdefined residuN
            residuN = similar(f1)  # Allocate once
        end
        get_resid!(residuN, uN, op_A, f1)  # In-place operation
        total_iters = ch1.iters
        true_rel_resid_norm = norm(residuN)/norm(f1)
        println("uN solved after $total_iters iterations, with True Relative Residual Norm = $true_rel_resid_norm")
        #set Neumann velocites:
        set_velocities!(model,uN)

         #get f2
        f2=get_rhs_dirichlet_inverse_data(model,inversion)

        # Define MSchur as a linear operator
        t_start = time()  # Get start time
       
        mi=gudata.ni+gvdata.ni+ghdata.n

        println("mi is " ,mi)
        println("gudata.ni is" ,gudata.ni)
        println("gvdata.ni is" ,gvdata.ni)
        println("ghdata.ni is" ,gudata.n)
        println("size of f2 is" ,size(f2))

        # Define the LinearMap with the modified schur_apply
        MSchur_op = LinearMap(x -> schur_apply(x, op_B, op_BT, A_diag_vals, op_C), mi, mi, ismutating = false)
        #  MSchur_op = LinearMap(x -> schur_apply(x, op_B, op_BT, A_diag_vals, op_C), mi, mi, ismutating = false)
        MSchur_new=get_op_diag_inversion(inversion, model, MSchur_op)
        t_end = time()  # Get end time
     #   println("MSchur_new linearmap computation took: ", t_end - t_start, " seconds")

        xD_guess=get_start_guess_dirichlet(inversion)

        #Dirichlet solve
        #Define right hand side for Schur (pressure) solve
        u0=xD_guess[1:gu.ni+gv.ni]
      #  u0 = gmres!(u0,op_A, f1,  abstol=inversion_params.abstol, maxiter=inversion_params.maxiter, restart=inversion_params.gmres_restart, log=false,verbose=false,Pl=M1)
        u0, chu0 = cg!(u0,op_A, f1,  abstol=inversion_params.abstol, maxiter=inversion_params.maxiter, log=true,verbose=false,Pl=M1)
        if !@isdefined residu0
            residu0 = similar(f1)  # Allocate once
        end
        get_resid!(residu0,u0,op_A,f1)
        relative_residual = norm(residu0) / norm(f1)
        total_iters_u0 = chu0.iters
        true_rel_resid_norm_u0 = norm(residuN)/norm(f1)
        println("u0 solved after $total_iters_u0 iterations, with True Relative Residual Norm = $true_rel_resid_norm_u0")

        # Pass tolerance for inner iterative solver
        inner_tol = inversion_params.inner_tol
        inner_maxiters=inversion_params.inner_maxiters

     #   bSchur=-op_B*u0+f2;
        bSchur = similar(f2,size(op_B, 1))  # Preallocate bSchur with the same size as f2
        mul!(bSchur, op_B, vec(u0))  # Compute bSchur = op_B * u0
        bSchur .*= -1           # In-place negation (bSchur = -op_B * u0)
        bSchur .+= f2           # In-place addition (bSchur = -op_B * u0 + f2)
        schur_output = similar(bSchur)  
        u = similar(bSchur,size(op_BT, 1))  # Buffer for PCG solve
        b = similar(bSchur,size(op_BT, 1))  # Buffer for right-hand side
        residBu = similar(bSchur,size(op_B, 1))  # Buffer for intermediate residuals
 
      SchurOp = LinearMap(x -> (schurVec!(schur_output, x, op_A, op_B, op_BT, op_C, M1, M2, inversion_params, u, b, residBu, inversion_params.inner_tol, inversion_params.inner_maxiters)), size(op_C, 1), size(op_C, 2))

    
       ni=gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n
        pD = similar(xD_guess, ni - (gu.ni + gv.ni))  # Preallocate correct size
       # pD=xD_guess[gu.ni+gv.ni+1:ni]
        fill!(pD,0)

        ch_pD = nothing
        t_start=time()  
        num_restarts = inversion_params.maxiter ÷ inversion_params.gmres_restart
        Pl = Diagonal(MSchur_new)  # Store once outside loop

        if inversion_params.gmres
                  # this one is to try and control stopping criteria:
            for i in 1:num_restarts
                #  mem_usage = @allocated pD, ch = gmres!(pD,SchurOp, bSchur,  abstol=inversion_params.abstol, maxiter=inversion_params.gmres_restart, restart=inversion_params.gmres_restart, Pl=Pl,verbose=false,log=true)
                pD, ch_pD =  gmres!(pD,SchurOp, bSchur, abstol=inversion_params.abstol, maxiter=inversion_params.gmres_restart, restart=inversion_params.gmres_restart, Pl=Pl,verbose=false,log=true)
                # println("Memory allocated: ", mem_usage, " bytes") 
                 total_iters_pD = ch_pD.iters
                # Compute true residual norm
                if !@isdefined residpD
                residpD = similar(bSchur)  # Allocate once
                end 
                get_resid!(residpD, pD, SchurOp, bSchur)  # In-place operation
                true_rel_resid_norm_pD = norm(residpD)/norm(bSchur)
                #compute total number of iterations:
                iter_total_pD=total_iters_pD+(i-1)*inversion_params.gmres_restart
                println("Loop $i : True Relative Residual Norm = $true_rel_resid_norm_pD")
                if  true_rel_resid_norm  < inversion_params.reltol
                println("Stopping early: Relative Residual below tolerance.")
                println("pD solved after $iter_total_pD iterations, with True Relative Residual Norm = $true_rel_resid_norm_pD")
                    break
                end
                
            end
            elseif inversion_params.cg              
                pD, ch_pD =  cg!(pD,SchurOp, bSchur,  abstol=inversion_params.abstol, maxiter=inversion_params.maxiter, Pl=Pl,verbose=false,log=true)
                total_iters = ch_pD.iters
                # Compute true residual norm
                if !@isdefined residpD
                residpD = similar(bSchur)  # Allocate once
                end 
                get_resid!(residpD, pD, SchurOp, bSchur)  # In-place operation
                true_rel_resid_norm_pD = norm(residpD)/norm(bSchur)
                total_iters_pD = ch_pD.iters
                println("pD solved after $total_iters_pD iterations, with True Relative Residual Norm = $true_rel_resid_norm_pD")
            else
                error("Error: iterative solver for pD isn't specified. Exiting..")
        end
        t_end = time()  # Get end time
        println("pD computation took: ", t_end - t_start, " seconds")

        #Define right hand side for velocity solve
#        b=f1-op_BT*pD;
        brhs = similar(f1, size(op_BT, 1))  # Ensure b has the same size as f1
        mul!(brhs, op_BT, vec(pD))  # Computes b = op_BT * pD in-place
        brhs .= f1 .- b # b=f1-b          Element-wise subtraction

        #Solve for Dirichlet velocity
        uD=xD_guess[1:gu.ni+gv.ni]
#        uD, ch_uD= cg!(uD,op_A,brhs,reltol=inversion_params.abstol, maxiter=inversion_params.maxiter,Pl=M1,log=true)
        uD, ch_uD= cg!(uD,op_A,brhs,abstol=inversion_params.abstol, maxiter=inversion_params.maxiter,Pl=M1,log=true)
        if !@isdefined residuD
            residuD = similar(brhs)  # Allocate once
        end
        residuD=get_resid!(residuD,uD,op_A,brhs)
        true_rel_resid_norm_uD = norm(residuD) / norm(b)
        total_iters_uD = ch_uD.iters
        println("uD solved after $total_iters_uD iterations, with True Relative Residual Norm = $true_rel_resid_norm_uD")

        x_output=[uD;pD]

        #set Dirichelt velocites:
        set_velocities!(inversion,x_output)

        #set pressues:
        set_inversion_pressure!(model,inversion,x_output)
end

"""
    get_op_diag(wavi::AbstractModel,op::LinearMap)

 Get diagonal of operator for use in preconditioner.

"""
function get_op_diag_inversion(inversion::AbstractModel,model::AbstractModel,op::LinearMap)
    @unpack gudata,gvdata, ghdata=inversion.data_fields
    @unpack params,solver_params=model
    mi,ni = size(op)
    
    @assert mi == ni == gudata.ni + gvdata.ni + ghdata.n
    op_diag=zeros(eltype(op),ni)
    ope_tmp=zeros(eltype(op),ni)
    sm=solver_params.stencil_margin

    sweep=[[1+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gudata.nxu, j=1:gudata.nyu][gudata.mask_inner];
           [1+sm^2+mod((i-1),sm)+sm*mod((j-1),sm) for i=1:gvdata.nxv, j=1:gvdata.nyv][gvdata.mask_inner];
           [1+2*sm^2+mod(i-1,sm)+sm*mod(j-1,sm) for i=1:ghdata.nxh, j=1:ghdata.nyh][ghdata.mask]]  # New `ghdata` part
          # sweep = vcat(sweep, zeros(Int, ghdata.n))  # Extend with zeros to match `gh` size
    e=zeros(Bool,ni)
    for i = unique(sweep)
        e .= sweep .== i
        mul!(ope_tmp,op,Float64.(e))
        op_diag[e] .= ope_tmp[e]
    end
    return op_diag
end


 # Define the schur_apply! function to accept the preallocated arrays
function schur_apply(x, op_B, op_BT, A_diag_vals, op_C)

    # Preallocate the arrays based on the size of `x` (input vector)
    temp = similar(x,size(op_BT, 1))         # Preallocate temp, the same size as x
    result = similar(x, size(op_B, 1))       # Preallocate result, the same size as x
    A_inv = similar(A_diag_vals)  # Preallocate A_inv, same size as A_diag_vals
    temp_scaled = similar(temp)  # Preallocate temp_scaled, the same size as x

    # Compute B^T * x
    mul!(temp, op_BT, vec(x))
    
    # Inverse of diagonal (A_diag_vals) and apply element-wise scaling
    A_inv .= 1.0 ./ A_diag_vals  # Inverse of diagonal, A_inv is updated in place
    
    # Compute the final result: -B * (A_inv * (B^T * x))
    mul!(temp_scaled, Diagonal(A_inv), temp)  # In-place scaling

    # Compute -B * (A_inv * (B^T * x))
    mul!(result, op_B, temp_scaled, -1.0, 0.0)  # Scale and overwrite

    # Add diagonal contribution of op_C directly to result
    mul!(result, op_C, vec(x), 1.0, 1.0)  # Efficient accumulation

    return result
end 
 

    function schurVec!(uc, x, A, B, BT, C, M1, M2, inversion_params, u, b, residBu, inner_tol, inner_maxiters)
    # Solve A * u = B' * x 
    mem_BT_x = @allocated mul!(b, BT, x)
  #  println("Memory allocated inside BT * x mul: ", mem_BT_x, " bytes")

    fill!(u,0)
    # Solve using CG with preconditioners M1 and M2 (if provided)
    cg!(u,A, b, reltol=inner_tol, maxiter=inner_maxiters,Pl=M1)

    # Compute final Schur complement operation
    #  return -B * u + C * x
        mul!(residBu, B, u)   # resid = B * u
        residBu .*= -1        # resid = -B * u
    #    uc=similar(x,  size(C, 1))
         mul!(uc, C, x)       # uC = C * x
        uc .+= residBu         # u = -B * u + C * x
    
end


"""
    get_start_guess_dirichlet(model::AbstractModel)

 Return starting guess used to begin iterative solution of Dirichlet problem.

"""
function get_start_guess_dirichlet(inversion)
    #still needs editing!!
    @unpack gu,gv,gh=inversion.fields
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
    @assert eltype(gu.u)==eltype(gv.v)
    x=[gu.samp_inner*gu.u[:];gv.samp_inner*gv.v[:]; gudata.samp_inner*gu.τsurf[:]; gvdata.samp_inner*gv.τsurf[:]; ghdata.samp*gh.σzzsurf[:]]
    return x
end


function get_op_A(model::AbstractModel{T,N}) where {T,N}
    @unpack gu,gv,gh=model.fields
    n_input = gu.ni + gv.ni
    n_output=gu.ni + gv.ni 
    op_fun_A! = get_op_fun_A(model)
    op_A=LinearMap{T}(
        op_fun_A!,           # The mutating function for A
        n_output,                  # Output size
        n_input,                 # Input size
        issymmetric=false,    # A is not symmetric in general
        ismutating=true,      # Function modifies input/output
        ishermitian=false,    # A is not Hermitian
        isposdef=false        # A is not positive definite
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
    ###added to match matlab...
    r_xy_strain_rate_sum :: Vector{T} = zeros(nxnyh);               @assert length(r_xy_strain_rate_sum) == nxnyh
    r_xy :: Vector{T} = zeros(nxnyh);                               @assert length(r_xy) == nxnyh
    ######
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
        ##EDIT TO MATCH WITH MATLAB!!!!
     #   @!  r_xy_strain_rate_sum = gc.cent*r_xy_strain_rate_sum_crop_c
     #   @!  r_xy = gh.dneghηav[]*r_xy_strain_rate_sum
     #   @.  r_xy = -r_xy
     #   @!  r_xy_c = gc.centᵀ*r_xy
     #CORRECT code for julia:
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

            #Resistive forces resolved in x and y directions
            fx .= d_rxx_dx .+ d_rxy_dy .- taubx 
            #.- h_d_extra_dx
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
#    R :: Array{T,2} = zeros(gh.nxh,gh.nyh);                               @assert length(R) == nxnyh
    beta_f1 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f1) == nxnyh
    beta_f2 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f2) == nxnyh
    R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
    R_crop :: Vector{T} = zeros(nxnyh);                             @assert length(R_crop) == nxnyh
    uspread_cent :: Vector{T} = zeros(nxnyh);                        @assert length(uspread_cent) == nxnyh
    u_on_h :: Vector{T} = zeros(nxnyh);                             @assert length(u_on_h) == nxnyh
    Ru_on_h :: Vector{T} = zeros(nxnyh);                            @assert length(Ru_on_h) == nxnyh
    Ru_on_u :: Vector{T} = zeros(nxnyu);                            @assert length(Ru_on_u) == nxnyu
 #   Ru_on_u_i:: Vector{T} = zeros(gu.ni);                           @assert length(Ru_on_u_i) == gu.ni
    vspread_cent :: Vector{T} = zeros(nxnyh);                       @assert length(vspread_cent) == nxnyh
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
    fx_sampi :: Vector{T} = zeros(ghdata.n);                        @assert length(fx_sampi) == ghdata.n
    opvecprod :: Vector{T} = zeros(gudata.ni+gvdata.ni+ghdata.n);   @assert length(opvecprod) == gudata.ni + gvdata.ni + ghdata.n
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
    
        #Compute R from Arthern et al 2015 and build the R term from equation (28):
     #   R[gh.mask] .= (1.0 .+ gh.β[gh.mask].*gh.quad_f1[gh.mask])./ (1.0 .+ gh.β[gh.mask].*gh.quad_f2[gh.mask])
      #  @! R_crop=gh.crop*vec(R)
         #R .= (1.0 .+ vec(gh.β).*vec(gh.quad_f1))./ (1.0 .+ vec(gh.β).*vec(gh.quad_f2))

         beta_f1 .= vec(gh.β).*vec(gh.quad_f1)
         beta_f1 .+=1
         beta_f2 .= vec(gh.β).*vec(gh.quad_f2)
         beta_f2 .+=1
        @. R = beta_f1 / beta_f2
        @! R_crop=gh.crop*R

        @! uspread_cent = gu.cent*uspread
        @! u_on_h=gh.crop*uspread_cent
        @. Ru_on_h=R_crop*u_on_h
        @! Ru_on_u=gu.centᵀ*Ru_on_h
      #  @! Ru_on_u_i=gu.samp_inner*Ru_on_u
        @! fx1=gu.crop*Ru_on_u
       #  fx1=Ru_on_u

       @! vspread_cent = gv.cent*vspread
        @! v_on_h=gh.crop*vspread_cent
        @. Rv_on_h=R_crop*v_on_h
        @! Rv_on_v=gv.centᵀ*Rv_on_h
    #    @! Rv_on_v=gv.samp_inner*Rv_on_v
        @! fy1=gv.crop*Rv_on_v
    #      fy1=Rv_on_v
       
        #sample these to only include data masked points:
        @! fx1_sampi = gudata.samp_inner*fx1
        @! fy1_sampi = gvdata.samp_inner*fy1

        #Then look at the div term in the B matrix - this will multiply the velocties:
        #-div(hu)

            hu.=vec(gu.h).*uspread
        @!  hu_crop=gu.crop*hu
        @!  divhu=gu.∂x*hu_crop
        @!  divhu_crop=gh.crop*divhu

            hv.=vec(gv.h).*vspread
        @!  hv_crop=gv.crop*hv
        @!  divhv=gv.∂y*hv_crop
        @!  divhv_crop=gh.crop*divhv

        @. fx = divhu_crop + divhv_crop
        @. fx = -fx 
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
    beta_f1 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f1) == nxnyh
    beta_f2 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f2) == nxnyh
  #  R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
    R_crop :: Vector{T} = zeros(nxnyh);                             @assert length(R_crop) == nxnyh
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
    vec_gvh :: Vector{T} = zeros(nxnyv);                      @assert length(vec_gvh) == nxnyv
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
       # @. R = (1.0 + gh.β[:]*gh.quad_f1[:])/ (1.0 + gh.β[:]*gh.quad_f2[:])

        beta_f1 .= vec(gh.β).*vec(gh.quad_f1)
        beta_f1 .+=1
        beta_f2 .= vec(gh.β).*vec(gh.quad_f2)
        beta_f2 .+=1
       @. R = beta_f1 / beta_f2
       @! R_crop=gh.crop*R

        @!  tauspreadu_c=gu.crop*(tausspreadu)
        @!  tauspreadu_on_h=gu.cent*(tauspreadu_c)
        @!  tauspreadu_c_on_h=gh.crop*(tauspreadu_on_h)
        @.  Rtauu_on_h=R_crop*tauspreadu_c_on_h
        @!  Rtauu_on_u=gu.centᵀ*Rtauu_on_h

        @!  tauspreadv_c=gv.crop*(tausspreadv)
        @!  tauspreadv_on_h=gv.cent*(tauspreadv_c)
        @!  tauspreadv_c_on_h=gh.crop*(tauspreadv_on_h)
        @.  Rtauv_on_h=R_crop*tauspreadv_c_on_h
        @!  Rtauv_on_v=gv.centᵀ*Rtauv_on_h

        @!  fx1_sampi = gu.samp_inner*Rtauu_on_u
        @!  fy1_sampi = gv.samp_inner*Rtauv_on_v

        #Then calculate grad term in the B matrix: hgrad(sigma_s)
        @!  sigmasspread_crop=gh.crop*sigmasspread
        @!  dsigmadx=gu.∂xᵀ*sigmasspread_crop
   #    @.  dsigmadx = -dsigmadx
        @!  dsigmadx_crop=gu.crop*dsigmadx
            h_dsigmadx.=vec(gu.h).*dsigmadx_crop
        @!  h_dsigmadx_crop=gu.crop*h_dsigmadx
        @.   h_dsigmadx_crop =-h_dsigmadx_crop

        @!  dsigmady=gv.∂yᵀ*sigmasspread_crop
    #    @.  dsigmady=-dsigmady
        @!  dsigmady_crop=gv.crop*dsigmady
        @views vec_gvh .= gv.h[:]  # no new memory allocation
        @.  h_dsigmady.=vec_gvh*dsigmady_crop
        @!  h_dsigmady_crop=gv.crop*h_dsigmady
        @.  h_dsigmady_crop = -h_dsigmady_crop
    
        @!  fx2_sampi = gu.samp_inner*h_dsigmadx_crop
        @!  fy2_sampi = gv.samp_inner*h_dsigmady_crop

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
    tauspreadu_on_h_crop :: Vector{T} = zeros(nxnyh);                     @assert length(tauspreadu_on_h_crop) == nxnyh
    tauspreadv_on_h_crop :: Vector{T} = zeros(nxnyh);                     @assert length(tauspreadv_on_h_crop) == nxnyh
    sigmasspread :: Vector{T} = zeros(nxnyh);                         @assert length(sigmasspread) == nxnyh
    sigmasspread_crop :: Vector{T} = zeros(nxnyh);                    @assert length(sigmasspread_crop) == nxnyh 
    R :: Vector{T} = zeros(nxnyh);                                  @assert length(R) == nxnyh
    R_plus_1 :: Vector{T} = zeros(nxnyh);                                  @assert length(R_plus_1) == nxnyh
    beta_f1 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f1) == nxnyh
    beta_f2 :: Vector{T} = zeros(nxnyh);                            @assert length(beta_f2) == nxnyh
    term_1 :: Vector{T} = zeros(nxnyh);                            @assert length(term_1) == nxnyh
    term_2 :: Vector{T} = zeros(nxnyh);                            @assert length(term_2) == nxnyh
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
     #   @.  R = (1.0 + gh.β[:]*gh.quad_f1[:])/ (1.0 + gh.β[:]*gh.quad_f2[:])
        
        beta_f1 .= vec(gh.β).*vec(gh.quad_f1)
        beta_f1 .+=1
        beta_f2 .= vec(gh.β).*vec(gh.quad_f2)
        beta_f2 .+=1
       @. R = beta_f1 / beta_f2

         R_plus_1 .= R .+ 1  # Add 1 to each element of R

         term_1 .= R_plus_1 .* vec(gh.quad_f1)  # (R + 1) * gh.quad_f1
         term_2 .= R .* vec(gh.quad_f2)  # R * gh.quad_f2
          P .= vec(gh.quad_f0) .- term_1 .+ term_2  # Compute final P
            
        @!  tausspreadu_crop=gu.crop*tausspreadu
        @!  tauspreadu_on_h=gu.cent*tausspreadu_crop
        @!  tauspreadu_on_h_crop=gh.crop*tauspreadu_on_h
        @.  Ptauu_on_h=P*tauspreadu_on_h_crop
        @!  Ptauu_on_h_crop=gh.crop*Ptauu_on_h
        @!  Ptauu_on_u=gu.centᵀ*Ptauu_on_h_crop
        @!  Ptauu_on_u_crop=gu.crop*Ptauu_on_u

        @!  tausspreadv_crop=gv.crop*tausspreadv
        @!  tauspreadv_on_h=gv.cent*tausspreadv_crop
        @!  tauspreadv_on_h_crop=gh.crop*tauspreadv_on_h
        @.  Ptauv_on_h=P*tauspreadv_on_h_crop
        @!  Ptauv_on_h_crop=gh.crop*Ptauv_on_h
        @!  Ptauv_on_v=gv.centᵀ*Ptauv_on_h_crop
        @!  Ptauv_on_v_crop=gv.crop*Ptauv_on_v

        fx1 = Ptauu_on_u_crop
        fy1 = Ptauv_on_v_crop

        #sample these to only include masked points:
        @!  fx1_sampi = gudata.samp_inner*fx1
        @!  fy1_sampi = gvdata.samp_inner*fy1

        fill!(fh_sampi, 0)

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
     
     us_data= zeros(gu.nxu,gu.nyu);                          
     vs_data= zeros(gv.nxv,gv.nyv);                               
     us_data_sampi  = zeros(gudata.ni);                      
     vs_data_sampi  = zeros(gvdata.ni);     

    
     dhdt_data = zeros(gh.nxh,gh.nyh);                              
     accumulation_data = zeros(gh.nxh,gh.nyh);     
     dhdtacc_data= zeros(nxnyh);                         
     dhdtacc_data_sampi= zeros(ghdata.n);                 
     f1 = zeros(gudata.ni+gvdata.ni+ghdata.n);            
     rhs_dirichlet = zeros(gudata.ni+gvdata.ni+ghdata.n);
     
     us_data = gudata.speed_u
     vs_data = gvdata.speed_v  

     us_data_vec=us_data[gudata.mask]
     vs_data_vec=vs_data[gvdata.mask]

#     us_data_spread=gudata.spread*us_data_vec
#    us_data_crop = gu.crop * us_data_spread[:]  
#    us_data_cent = gu.cent * us_data_crop 
#    us_data_cent_samp = gh.samp  * us_data_cent
#    us_data_cent_samp_spread= gh.spread*us_data_cent_samp
    #
 #   us_data_centT = gu.centᵀ*us_data_cent_samp_spread
 #   us_data_centT_crop = gu.crop*us_data_centT

  #  vs_data_spread=gvdata.spread*vs_data_vec
  #  vs_data_crop = gv.crop * vs_data_spread[:]    
  #  vs_data_cent = gv.cent * vs_data_crop 
  #  vs_data_cent_samp = gh.samp  * vs_data_cent 
  #  vs_data_cent_samp_spread= gh.spread*vs_data_cent_samp
   # vs_data_centT = gv.centᵀ*vs_data_cent_samp_spread
  #  vs_data_centT_crop = gv.crop*vs_data_centT

   
    # us_data_sampi = gudata.samp_inner*us_data_centT_crop
    # vs_data_sampi =gvdata.samp_inner*vs_data_centT_crop

    us_data_sampi = us_data_vec
     vs_data_sampi = vs_data_vec

          f1[1:gudata.ni] .= us_data_sampi
          f1[(gudata.ni+1):(gudata.ni+gvdata.ni)] .= vs_data_sampi



        
     dhdt_data = ghdata.dhdt
     accumulation_data = ghdata.accumulation_rate
      dhdtacc_data=dhdt_data[:] .- accumulation_data[:]
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
    surf_crop= gh.crop*(gh.s[:])

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

        get_rhs_dirichlet!(f3,model)

        rhs .= f1 .+ f2 .+ f3

    return rhs
end

function get_rhs_neumann_data_check(model, inversion)
    @unpack gh,gu,gv,gc=model.fields
    @unpack params, solver_params = model
    @unpack gudata, gvdata, ghdata=inversion.data_fields 
       
     #Preallocate intermediate variables used 
     nxnyh = gh.nxh*gh.nyh
     nxnyu = gu.nxu*gu.nyu
     nxnyv = gv.nxv*gv.nyv
     
     us_data= zeros(gu.nxu,gu.nyu);                          
     vs_data= zeros(gv.nxv,gv.nyv);                               
     us_data_sampi  = zeros(gudata.ni);                      
     vs_data_sampi  = zeros(gvdata.ni);                    
     dhdt_data = zeros(gh.nxh,gh.nyh);                              
     accumulation_data = zeros(gh.nxh,gh.nyh);     

     dhdtacc_data= zeros(nxnyh);                         
     dhdtacc_data_sampi= zeros(ghdata.n);                 
     f1 = zeros(gudata.ni+gvdata.ni+ghdata.n);            
     rhs_dirichlet = zeros(gudata.ni+gvdata.ni+ghdata.n);
    
     
     us_data = gudata.speed_u
     vs_data = gvdata.speed_v  

     us_data_vec=us_data[gudata.mask]
     vs_data_vec=vs_data[gvdata.mask]

     us_data_spread=gudata.spread*us_data_vec
    us_data_crop = gu.crop * us_data_spread[:]  
    us_data_cent = gu.cent * us_data_crop 
    us_data_cent_samp = gh.samp  * us_data_cent
    us_data_cent_samp_spread= gh.spread*us_data_cent_samp
    #
    us_data_centT = gu.centᵀ*us_data_cent_samp_spread
    us_data_centT_crop = gu.crop*us_data_centT

    vs_data_spread=gvdata.spread*vs_data_vec
    vs_data_crop = gv.crop * vs_data_spread[:]    
    vs_data_cent = gv.cent * vs_data_crop 
    vs_data_cent_samp = gh.samp  * vs_data_cent 
    vs_data_cent_samp_spread= gh.spread*vs_data_cent_samp
    vs_data_centT = gv.centᵀ*vs_data_cent_samp_spread
    vs_data_centT_crop = gv.crop*vs_data_centT


     us_data_sampi = gudata.samp_inner*us_data[:]
     vs_data_sampi =gvdata.samp_inner*vs_data[:]

   
     f1[1:gudata.ni] .= us_data_sampi
     f1[(gudata.ni+1):(gudata.ni+gvdata.ni)] .= vs_data_sampi
    
     dhdt_data = ghdata.dhdt
     accumulation_data = ghdata.accumulation_rate
      dhdtacc_data=dhdt_data[:] .- accumulation_data[:]
     dhdtacc_data_sampi=ghdata.samp*dhdtacc_data
    
    f1[(gudata.ni+gvdata.ni+1):(gudata.ni+gvdata.ni+ghdata.n)] .=  dhdtacc_data_sampi
        
    rhs_dirichlet .= f1 

    return rhs_dirichlet
end
    
"""
    set_velocities!(model::AbstractModel,x)

Set velocities to particular values. Input vector x represents stacked u and v components at valid grid points.
"""
function set_inversion_pressure!(model::AbstractModel,inversion::AbstractModel,x)
    @unpack gh,gu,gv,gc=model.fields
    @views inversion.fields.gu.τsurf[inversion.data_fields.gudata.mask] .= x[(gu.ni+gv.ni+1):(gu.ni+gv.ni+inversion.data_fields.gudata.ni)]
    @views inversion.fields.gv.τsurf[inversion.data_fields.gvdata.mask] .= x[(gu.ni+gv.ni+inversion.data_fields.gudata.ni+1):(gu.ni+gv.ni+inversion.data_fields.gudata.ni+inversion.data_fields.gvdata.ni)]
    @views inversion.fields.gh.σzzsurf[inversion.data_fields.ghdata.mask] .= x[(gu.ni+gv.ni+inversion.data_fields.gudata.ni+inversion.data_fields.gvdata.ni+1):(gu.ni+gv.ni+inversion.data_fields.gudata.ni+inversion.data_fields.gvdata.ni+inversion.data_fields.ghdata.n)]
    return inversion
end

"""
    update_velocities_on_h_grid_dirichlet!(model::AbstractModel)

Update the velocities (depth averaged, surface and bed) on the h grid dirichelt 
"""
function update_velocities_on_h_grid_dirichlet!(model)
    @unpack gh,gu,gv = model.fields
    #depth averaged velocities
    gh.u[:] .= gu.cent*gu.u[:] #(gu.u[1:end-1,:] + gu.u[2:end,:])./2
    gh.v[:] .= gv.cent*gv.v[:] #(gv.v[:,1:end-1] + gv.v[:, 2:end])./2

    #bed velocities
#    gh.ub[gh.mask] .=(gh.u[gh.mask].-(gh.quad_f1[gh.mask] .- gh.quad_f2[gh.mask]).*gh.τx_surf[gh.mask])./(1.0 .+ gh.quad_f2[gh.mask] .* gh.β[gh.mask])
    gh.ub .=(gh.u.-(gh.quad_f1 .- gh.quad_f2).*gh.τx_surf)./(1.0 .+ gh.quad_f2 .* gh.β)
    gh.vb .= (gh.v.-(gh.quad_f1 .- gh.quad_f2).*gh.τy_surf)./(1.0 .+ gh.quad_f2 .* gh.β);
    gh.bed_speed  .= sqrt.(gh.ub.^2 .+gh.vb.^2);

    #surface velocities
    gh.us .=  gh.ub.*(1.0 .+ gh.quad_f1 .* gh.β) .+ (gh.quad_f0 .- gh.quad_f1).*gh.τx_surf
    gh.vs .=  gh.vb.*(1.0 .+ gh.quad_f1 .* gh.β) .+ (gh.quad_f0 .- gh.quad_f1).*gh.τy_surf

    return model
end

"""
    update_surf_stress_dirichelt!(model::AbstractModel)

Find the surface stress dirichelt on the h-grid.
"""
function update_surf_stress_dirichelt!(model::AbstractModel)
    @unpack gv,gu,gh=model.fields
    gh.τx_surf[gh.mask] .= gh.samp*(gu.cent*gu.τsurf[:])
    gh.τy_surf[gh.mask] .= gh.samp*(gv.cent*gv.τsurf[:])
    gh.τsurf[gh.mask]  .= sqrt.(gh.τx_surf[gh.mask].^2 .+gh.τx_surf[gh.mask].^2)
    return model
end

"""
    update_surf_speed!(model::AbstractModel)

Find the sliding speed on the h-grid using the speed components.
"""
function update_surf_speed!(model::AbstractModel)
    @unpack gh=model.fields
    gh.surf_speed  .= sqrt.(gh.us.^2 .+gh.vs.^2);
    return model
end

"""
    update_basal_drag_components!(model::AbstractModel)

Find the shear stress at the bed.
"""
function update_basal_drag_components!(model::AbstractModel)
    @unpack gh,gu,gv=model.fields
    gh.τx_bed .= gh.β .* gh.ub
    gh.τy_bed .= gh.β .* gh.vb
    gh.τbed .= gh.β .* gh.bed_speed
    return model
end

"""
    update_rheological_operators_inverse!(model::AbstractModel)

Precompute various diagonal matrices used in defining the momentum operator.
"""
function update_rheological_operators_inverse!(model::AbstractModel,inversion)
    @unpack gh,gu,gv,gc = inversion.fields
    @unpack params, solver_params = model
    gh.dneghηav[] .= gh.crop*Diagonal(-gh.h[:].*gh.ηav[:])*gh.crop
    gc.dneghηav[] .= gc.crop*Diagonal(-gh.cent_xy*(gh.h[:].*gh.ηav[:]))*gc.crop
    gu.dnegβeff[] .= gu.crop*Diagonal(-gu.βeff[:])*gu.crop
    gv.dnegβeff[] .= gv.crop*Diagonal(-gv.βeff[:])*gv.crop
    gh.dimplicit[] .= gh.crop*Diagonal(-params.density_ice * params.g * solver_params.super_implicitness .* params.dt * gh.dsdh[:])*gh.crop
    return model
end

"""
    update_shelf_heating!(model::AbstractModel)

Find the shelf_heating on the h-grid.
"""
function update_shelf_heating!(model)
    @unpack gv,gu,gh=model.fields
    gh.shelf_heating[gh.mask] .= 4. *gh.h[gh.mask].*(gh.ηav[gh.mask].*(gh.shelf_strain_rate[gh.mask].^2));
    return model
end


"""
    update_shelf_heating_dirichlet!(model::AbstractModel)

Find the shelf_heating on the h-grid for dirichelt.
"""
function update_shelf_heating_dirichlet!(model,inversion)
    @unpack gv,gu,gh=model.fields
    inversion.fields.gh.shelf_heating[gh.mask] .= 4. *gh.h[gh.mask].*(gh.ηav[gh.mask].*(inversion.fields.gh.shelf_strain_rate[gh.mask].^2));
    return inversion
end

"""
    update_vert_shear_heating!(model::AbstractModel)

Find the vert_shear_heating on the h-grid.
"""
function update_vert_shear_heating!(model::AbstractModel)
    @unpack gv,gu,gh=model.fields
    gh.vert_shear_heating[gh.mask] .= ((gh.τbed[gh.mask]).^2).*gh.quad_f2[gh.mask]
    return model
end

"""
    update_vert_shear_heating_dirichlet!(model::AbstractModel)

Find the vert_shear_heating on the h-grid dirichlet
"""
function update_vert_shear_heating_dirichlet!(model,inversion)
    @unpack gv,gu,gh=model.fields
    τxy_surf_f0=(inversion.fields.gh.τx_surf[gh.mask].^2 .+ inversion.fields.gh.τy_surf[gh.mask].^2).*gh.quad_f0[gh.mask]
    τx_surfbed_f1=2. *inversion.fields.gh.τx_surf[gh.mask].*(inversion.fields.gh.τx_bed[gh.mask].-inversion.fields.gh.τx_surf[gh.mask]).*gh.quad_f1[gh.mask]
    τy_surfbed_f1=2. *inversion.fields.gh.τy_surf[gh.mask].*(inversion.fields.gh.τy_bed[gh.mask].-inversion.fields.gh.τy_surf[gh.mask]).*gh.quad_f1[gh.mask]
    τxy_surfbed_f2=((inversion.fields.gh.τx_bed[gh.mask] .- inversion.fields.gh.τx_surf[gh.mask]).^2 .+ (inversion.fields.gh.τy_bed[gh.mask] .- inversion.fields.gh.τy_surf[gh.mask]).^2).*gh.quad_f2[gh.mask]
    inversion.fields.gh.vert_shear_heating[gh.mask] .= τxy_surf_f0 .+ τx_surfbed_f1 .+ τy_surfbed_f1 .+ τxy_surfbed_f2
    return inversion
end

"""
    update_vert_shear_heating!(model::AbstractModel)

Find the vert_shear_heating on the h-grid.
"""
function update_drag_heating!(model::AbstractModel)
    @unpack gv,gu,gh=model.fields
    gh.drag_heating[gh.mask] .= (gh.τx_bed[gh.mask].*gh.ub[gh.mask]) .+ (gh.τy_bed[gh.mask].*gh.vb[gh.mask]);
    return model
end


"""
    update_β_inversion!(model::AbstractModel)

Update the drag coefficient at the bed in the inversion
"""
function update_β_inversion!(model::AbstractModel,inversion)
    @unpack gh=model.fields
    @unpack inversion_params=inversion
    gh.β[gh.mask] .=  gh.β[gh.mask].*(gh.drag_heating[gh.mask]./inversion.fields.gh.drag_heating[gh.mask]).^inversion_params.βpower
    return model
end

"""
    update_preBfactor_inversion!(model::AbstractModel)

Update the preBfactor in the inversion
"""
function update_preBfactor_inversion!(model::AbstractModel,inversion)
    @unpack gh=model.fields
    @unpack inversion_params=inversion
    BPower=inversion_params.Bpower_shelf.*ones(gh.nxh,gh.nyh)
    #### WHAT ABOUT PARTIALLY FLOATING?!
    aground = (gh.haf .>= 0)
    BPower[aground] .= inversion_params.Bpower_grounded

    gh.preBfactor[gh.mask] .= gh.preBfactor[gh.mask].*((gh.shelf_heating[gh.mask] .+ gh.vert_shear_heating[gh.mask])./(inversion.fields.gh.shelf_heating[gh.mask] .+ inversion.fields.gh.vert_shear_heating[gh.mask])).^BPower[gh.mask];
    return model
end


function update_preBfactor_3d!(model::AbstractModel)
    @unpack gh, g3d = model.fields
    
    # Iterate through all grid points in the horizontal directions (i,j)
    for j = 1:g3d.nys
        for i = 1:g3d.nxs
            if gh.mask[i,j]
                # Calculate the total weight for the vertical layers (integral over sigma)
               # total_weight = sum(g3d.quadrature_weights)  # This is the sum of the quadrature weights
                # Loop through the vertical layers (k)
                for k = 1:g3d.nσs
                    # Adjust preBfactor for each (i, j, k) so that the weighted sum gives the desired depth-averaged value
                 #   g3d.preBfactor[i,j,k] = gh.preBfactor[i,j] * g3d.quadrature_weights[k] / total_weight
                    g3d.preBfactor[i,j,k] = gh.preBfactor[i,j]/(g3d.nσs*g3d.quadrature_weights[k])  
                end
            end
        end
    end
    
    preBfactor_av = zeros(gh.nxh,gh.nyh)
    for k=1:g3d.nσs
       for j = 1:g3d.nys
          for i = 1:g3d.nxs
            if gh.mask[i,j]
                preBfactor_av[i,j] += g3d.quadrature_weights[k] * g3d.preBfactor[i,j,k]
            end
          end
       end
    end

    return model
end



"""
    update_damage!(model::AbstractModel)

update to damage in the inversion on the 3d grid at all sigma levels.
"""
function update_damage!(model::AbstractModel)
    @unpack gh,g3d=model.fields
  #  @unpack params,solver_params=model 
      for k=1:g3d.nσs
        for j=1:g3d.nys
            for i=1:g3d.nxs
                if gh.mask[i,j]
              #  g3d.Φ[i,j,k] = 1.0 .- g3d.preBfactor[i,j,k]
                g3d.Φ[i,j,k] = 1.0 .- gh.preBfactor[i,j]
                end
            end
        end
    end

    return model
end


"""
    update_glen_b!(model::AbstractModel)

Update  glen_b
"""
function update_glen_b!(model::AbstractModel)
    @unpack gh,g3d=model.fields
    @unpack params=model
 #   @unpack gh_inversion=inversion.fields
 #   @unpack inversion_params=inversion
    #####MASKS OR NOT HERE??
    for i = 1:g3d.nxs
        for j = 1:g3d.nys
            for k = 1:g3d.nσs
                if gh.mask[i,j]
                g3d.glen_b[i,j,k] = glen_b.(g3d.θ[i,j,k],g3d.Φ[i,j,k],params.glen_a_ref[i,j], params.glen_n, params.glen_a_activation_energy, params.glen_temperature_ref, params.gas_const)
                end
            end
        end
    end
    
    return model
end

"""
    update_viscosities_quadratures_dirichlet(model::AbstractModel,inversion)

Update the viscosities and quadrates for the dirichlet solution by setting as equal to the Neumann (model) fields
"""

function update_viscosities_quadratures_dirichlet(model::AbstractModel,inversion)
    inversion.fields.gh.ηav[inversion.fields.gh.mask] =model.fields. gh.ηav[model.fields.gh.mask]
    inversion.fields.gh.quad_f0[inversion.fields.gh.mask] =model.fields.gh.quad_f0[model.fields.gh.mask]
    inversion.fields.gh.quad_f1[inversion.fields.gh.mask] =model.fields.gh.quad_f1[model.fields.gh.mask]
    inversion.fields.gh.quad_f2[inversion.fields.gh.mask] =model.fields.gh.quad_f2[model.fields.gh.mask]
 for k=1:model.fields.g3d.nσs
    inversion.fields.g3d.η[inversion.fields.gh.mask,k] =model.fields.g3d.η[model.fields.gh.mask,k]
 end
    return inversion
end


"""
    update_JKV!(model::AbstractModel)

Find JKV and store
"""
function update_JKV!(model::AbstractModel,inversion,clock)
    @unpack gh=model.fields
    @unpack grid=model
    if inversion.data_fields.ghdata.n==0
        σzzsurf_comp=0
     #   println("Not using dhdt so setting  σzzsurf_comp=0 ")
    else
    σzzsurf_comp=inversion.fields.gh.σzzsurf[gh.mask]'*(inversion.data_fields.ghdata.dhdt[gh.mask].-gh.dhdt[gh.mask]).*(grid.dx*grid.dy) 
    end
    τx_surf_comp=inversion.fields.gh.τx_surf[gh.mask]'*(inversion.fields.gh.us[gh.mask].-gh.us[gh.mask]).*(grid.dx*grid.dy)
    τy_surf_comp=inversion.fields.gh.τy_surf[gh.mask]'*(inversion.fields.gh.vs[gh.mask].-gh.vs[gh.mask]).*(grid.dx*grid.dy)

    JKV=σzzsurf_comp .+ τx_surf_comp .+ τy_surf_comp
  #  inversion.inversion_output.JKV[clock.n_iter+1]=JKV
    push!(inversion.inversion_output.JKV, JKV)
    println("   JKV is " ,inversion.inversion_output.JKV)
    return model
end

"""
    update_JRMS!(model::AbstractModel)

Find JRMS and store
"""
function update_JRMS!(model::AbstractModel,inversion,clock)
    @unpack gu,gv,gh=model.fields
    @unpack gudata,gvdata,ghdata=inversion.data_fields
    @unpack grid=model

    surf_speed_data_on_h=zeros(gh.nxh,gh.nyh)

    us_data = gudata.speed_u[gudata.mask]
    vs_data = gvdata.speed_v[gvdata.mask] 

    us_data_spread =  gudata.crop*(gudata.spread*us_data)
    vs_data_spread =  gvdata.crop*(gvdata.spread*vs_data)

    us_data_cent =gh.crop*(gu.cent*(us_data_spread))
    vs_data_cent =gh.crop*(gv.cent*(vs_data_spread))

    surf_speed_data_on_h[:]=sqrt.(us_data_cent[:].^2 .+ vs_data_cent[:].^2)

    #save into ghdata structure:
    ghdata.us[:] .= us_data_cent
    ghdata.vs[:] .= vs_data_cent
    ghdata.surf_speed .= surf_speed_data_on_h

    # Assume u_mask and v_mask are boolean arrays (true = valid, false = invalid)
    surf_speed_data_on_h_mask = falses(size(gh.mask))  # Initialize h_mask with all false

    # Check validity of the surrounding u and v masks and update h_mask accordingly
    surf_speed_data_on_h_mask[:, :] .= 
    (gudata.mask[2:end, :] .& gudata.mask[1:end-1,:]) .&  # Check valid u-mask for left and right neighbors
    (gvdata.mask[:, 2:end] .& gvdata.mask[:, 1:end-1])       # Check valid v-mask for top and bottom neighbors

    data_and_model_mask = surf_speed_data_on_h_mask .& gh.mask
    data_and_model_mask = convert(Array{Bool,2},data_and_model_mask)

    #save mask on h for which surf data is valid in the ghdata structure:
    ghdata.surf_speed_mask .= surf_speed_data_on_h_mask

    JRMS=sqrt(mean((surf_speed_data_on_h[data_and_model_mask] .- gh.surf_speed[data_and_model_mask]).^2));
   # inversion.inversion_output.JRMS[clock.n_iter+1]=JRMS
    push!(inversion.inversion_output.JRMS, JRMS)
    println("   JRMS is " ,inversion.inversion_output.JRMS)
    return model
end

function update_surface_velocities_on_uv_grid!(model)
    @unpack gh,gu,gv = model.fields
    #surface  velocities
    gu.us[:].=gu.crop*(gu.centᵀ*gh.crop*(gh.us[:]))
    gv.vs[:].=gv.crop*(gv.centᵀ*gh.crop*(gh.vs[:]))

    return model
end

"""
    set_residual_inversion!(model::AbstractModel,residual)

Set residuals to particular values.
"""
function set_residual_inversion!(model::AbstractModel,inversion,residual)
    @unpack gu,gv=model.fields
    @unpack gudata,gvdata,ghdata=inversion.data_fields
#### WILL ERROR if dhdt not provided!!!!! 

    @views inversion.fields.gu.residual[gu.mask_inner] .= residual[1:gu.ni]
    @views inversion.fields.gv.residual[gv.mask_inner] .= residual[(gu.ni+1):(gu.ni+gv.ni)]
    @views gudata.residual[gudata.mask_inner] .= residual[(gu.ni+gv.ni+1):(gu.ni+gv.ni+gudata.ni)]
    @views gvdata.residual[gvdata.mask_inner] .= residual[(gu.ni+gv.ni+gudata.ni+1):(gu.ni+gv.ni+gudata.ni+gvdata.ni)]
    @views ghdata.residual[ghdata.mask] .= residual[(gu.ni+gv.ni+gudata.ni+gvdata.ni+1):(gu.ni+gv.ni+gudata.ni+gvdata.ni+ghdata.n)]
    return model, inversion
end
