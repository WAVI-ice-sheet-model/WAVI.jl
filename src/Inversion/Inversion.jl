struct Inversion{M,IP,D,IO}
    model::M
    data::D
    inversion_params::IP
    inversion_output::IO
end

"""
        Inversion(;
            model = Model(),
            data = data,
            inversion_params = InversionParams(),
            inversion_output = InversionOutput(),
            )

Construct a WAVI.jl inversion object.

Keyword arguments
=================
    - `model`: model
    - `data`: data to be used in the inversion
    - `inversion_params`: inversion parameters
    - `inversion_output`: output from the inversion

"""
function Inversion(;
    model = Model(),
    data = DataFields(),
    inversion_params = InversionParams(),
    inversion_output = InversionOutput())

    if !all(data.ghdata.accumulation_rate .== model.params.accumulation_rate)
    println("WARNING: a different accumulation_rate is being used in model and inversion.data!")
    else
    #println("The same accumulation_rate is being used in model and inversion.data!")
    end

    inversion=Inversion(model,data,inversion_params,inversion_output)

    return inversion
end

include("inversion_utilities.jl")
