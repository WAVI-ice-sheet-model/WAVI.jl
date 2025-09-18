# TODO: Utility functions from full WAVI - highlights module issue (would have to import everything)
function get_bed_elevation(bed_elevation::F, grid) where (F <: Function)
    bed_array = bed_elevation.(grid.xxh, grid.yyh)
    return bed_array
end

function get_bed_elevation(bed_elevation::Array{T,2}, grid) where (T <: Real)
    bed_array = bed_elevation
    return bed_array
end
