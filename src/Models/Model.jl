struct Model{T <: Real, N <: Integer,A, G, M <:AbstractMeltRate, PS <: AbstractParallelSpec, SL <:AbstractSlidingLaw, BH <:AbstractBasalHydrology, TD <: AbstractThermoDynamics} <: AbstractModel{T,N,M,PS,SL,BH,TD}
    grid::Grid{T,N}
    params::Params{T,A,G}
    solver_params::SolverParams{T,N}
    initial_conditions::InitialConditions{T}
    fields::Fields{T,N}
    shelf_melt_rate::M
    parallel_spec::PS
    sliding_law::SL
    basal_hydrology::BH
    thermo_dynamics::TD
end

"""
    Model(;
        grid = nothing, 
        bed_elevation = nothing,
        params = Params(),
        solver_params = SolverParams(),
        initial_conditions = InitialConditions(),
        shelf_melt_rate = UniformMeltRate(),
        parallel_spec = BasicParallelSpec(),
        sliding_law = SlidingLaw(),
        basal_hydrology = BasalHydrology(),
        thermo_dynamics = ThermoDynamics())

Construct a WAVI.jl model object.

Keyword arguments
=================
- `grid`: (required) an instance of a `Grid` object, which defines the computational grid
- `bed_elevation`: (required) an array of size `grid.nx` x `grid.ny` which defines the bed elevation
- `params`: a `Params` object that defines physical parameters 
- `solver_params`: a `SolverParams` object that defines parameters relating to the numerical scheme
- `initial_conditions`: an `InitialConditions` object that (optionally) defines the initial ice thickness, temperature, viscosity, damage, 
                        basal water thickness, effective pressure, grounded basal melt rate, and depth-averaged temperature
- `shelf_melt_rate`: a shelf melt rate model, responsible for controlling/setting the basal melt rate under ice shelves
- `parallel_spec`: specification of parallel computation method
- `sliding_law`: a sliding law model, responsible for controlling/setting the basal friction.
- `basal_hydrology`: a basal hydrology model, responsible for calculating the basal water thickness and effective pressure.
- `thermo_dynamics`: a thermodynamics model, responsible for calculating the ice temperature and grounded basal melt rate.

"""
function Model(;
    grid = nothing, 
    bed_elevation = nothing,
    params = Params(),
    solver_params = SolverParams(),
    initial_conditions = InitialConditions(),
    shelf_melt_rate = UniformMeltRate(),
    parallel_spec = BasicParallelSpec(),
    sliding_law = WeertmanSlidingLaw(),
    basal_hydrology = NoHydrology(),
    thermo_dynamics = NoThermoDynamics())

    #check that a grid and bed has been inputted
    ~(grid === nothing) || throw(ArgumentError("You must specify an input grid"))
    ~(bed_elevation === nothing) || throw(ArgumentError("You must input a bed elevation"))
    
    #check types
    #if a functional bed has been specified, convert to an array
    bed_array = zeros(grid.nx,grid.nx) #initialize a bed array
    try
    bed_array = get_bed_elevation(bed_elevation, grid)
    catch
    throw(ArgumentError("bed elevation must be of type function or array"))
    end
            #check the size of the bed
    #(Base.size(bed_array) = (grid.nx, grid.ny)) || throw(ArgumentError("Bed and grid sizes must be identical"))
    
    #Check initial conditions, and revert to default values if not
    initial_conditions = check_initial_conditions(initial_conditions, params, grid)

    ## Parameter fields checks 
    # if drag_coefficient passed as a scalar, replace drag_coefficient parameters with matrix of this value
    if isdefined(sliding_law, :drag_coefficient)
        if isa(sliding_law.drag_coefficient, Number) 
            sliding_law = @set sliding_law.drag_coefficient = sliding_law.drag_coefficient*ones(grid.nx,grid.ny)
        end
        #check size compatibility of resulting drag C
        (size(sliding_law.drag_coefficient)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input drag_coefficient must match grid size (i.e. $(grid.nx) x $(grid.ny))"))
    end

    #if coulomb_coefficient passed as a scalar, replace coulomb_coefficient parameters with matrix of this value
    if isdefined(sliding_law, :coulomb_coefficient)
        if isa(sliding_law.coulomb_coefficient, Number)
            sliding_law = @set sliding_law.coulomb_coefficient = sliding_law.coulomb_coefficient * ones(grid.nx, grid.ny)
        end
        #check size compatibility of resulting drag C
        (size(sliding_law.coulomb_coefficient)==(grid.nx, grid.ny)) || throw(DimensionMismatch("Size of input coulomb_coefficient must match grid size (i.e. $(grid.nx) x $(grid.ny))"))
    end

    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.accumulation_rate, Number) 
        params = @set params.accumulation_rate = params.accumulation_rate*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting accumulation rate
    (size(params.accumulation_rate)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input accumulation must match grid size (i.e. $(grid.nx) x $(grid.ny))"))

    #if accumulation is passed as a scalar, replace accumulation parameters with matrix of this value
    if isa(params.glen_a_ref, Number) 
        params = @set params.glen_a_ref = params.glen_a_ref*ones(grid.nx,grid.ny)
    end
    #check size compatibility of resulting glen a ref
    (size(params.glen_a_ref)==(grid.nx,grid.ny)) || throw(DimensionMismatch("Size of input glen_a_ref must match grid size (i.e. $(grid.nx) x $(grid.ny))"))




    #Setup the fields 
    fields = setup_fields(grid, initial_conditions, solver_params, params, bed_array)

    #Use type constructor to build initial state with no extra physics
    model=Model(grid,params,solver_params,initial_conditions,fields,shelf_melt_rate,parallel_spec,sliding_law,basal_hydrology,thermo_dynamics)

    return model
end

include("model_utilities.jl")
include("update_state.jl")
include("update_velocities.jl")