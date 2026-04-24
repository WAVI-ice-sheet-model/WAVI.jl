struct BinfileMeltRate <: AbstractMeltRate
    input_filename::String    #specify melt filename
end

function BinfileMeltRate(;
                        input_filename = nothing)

    #input file exceptions
    ~(input_filename === nothing) || throw(ArgumentError("You must pass an input filename"))
    isfile(input_filename) || throw(ArgumentError("Did not find the specified binary file"))
    
    return BinfileMeltRate(input_filename)
end


function update_shelf_melt_rate!(shelf_melt_model::BinfileMeltRate, fields, grid, clock)
    @unpack input_filename = melt_model
    @unpack shelf_basal_melt,grounded_fraction = fields.gh
    
    file_size = stat(input_filename).size #bytes in input file must match matrix returned
    shelf_melt_rate = zeros(grid.nx,grid.ny)
    (file_size == sizeof(shelf_melt_rate)) || throw(DimensionMismatch("Size of input file incompatible with specified nx, ny"))
    try 
        read!(input_filename, shelf_melt_rate)
    catch
        Error("Input file read error")
    end

   # shelf_melt_rate .= ntoh.(shelf_melt_rate)
    shelf_basal_melt[:] = shelf_melt_rate[:] .* ((1 .- grounded_fraction))
    return nothing
end

