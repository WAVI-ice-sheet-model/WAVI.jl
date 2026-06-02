export BinfileAccumulationRate

struct BinfileAccumulationRate <: AbstractSurfaceMassBalance
    input_filename::String    #specify accumulation rate filename
end

function BinfileAccumulationRate(;
                        input_filename = nothing)

    #input file exceptions
    ~(input_filename === nothing) || throw(ArgumentError("You must pass an input filename"))
    isfile(input_filename) || throw(ArgumentError("Did not find the specified binary file"))
    
    return BinfileAccumulationRate(input_filename)
end


function update_accumulation_rate!(surface_mass_balance::BinfileAccumulationRate, model::AbstractModel, clock::Clock)
    @unpack input_filename = surface_mass_balance
    @unpack accumulation = model.fields.gh
    
    file_size = stat(input_filename).size #bytes in input file must match matrix returned
    accumulation_rate_from_file = zeros(grid.nx,grid.ny)
    (file_size == sizeof(accumulation_rate_from_file)) || throw(DimensionMismatch("Size of input file incompatible with specified nx, ny"))
    try 
        read!(input_filename, accumulation_rate_from_file)
    catch
        Error("Input file read error")
    end

   # accumulation_rate .= ntoh.(accumulation_rate)
    accumulation .= accumulation_rate_from_file
    return nothing
end

