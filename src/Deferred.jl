module Deferred

export Collector, clear!, collect!, register_field!

using WAVI: AbstractModel

abstract type AbstractCollector end
abstract type AbstractProperty end

function get_property(item::AbstractProperty, model::AbstractModel) end

mutable struct Collector <: AbstractCollector
    items::Dict{String, AbstractProperty}
    data::Union{Nothing, Dict{String, Any}}
    
    function Collector()
        new(Dict{String, AbstractProperty}(), Dict{String, Any}())
    end
end

function clear!(collector::Collector)
    collector.data = nothing
end

function collect!(collector::Collector, model::AbstractModel)
    data = Dict{String, Any}()
    
    for property_name in keys(collector.items)
        data[property_name] = get_property(collector.items[property_name], model)
    end
    collector.data = data
end

function register_item!(collector::Collector, property::AbstractProperty)
    collector.items[property.name] = property
end

# Specific implementation of a AbstractProperty for fields, allowing deferred collection of nested fields
struct FieldExtractor{T} <: AbstractProperty
    name::String
    field_accessor::Function
    path::Vector{Symbol}
end

function field_extractor(name::String, accessor::Function, path::Vector{Symbol})
    return FieldExtractor{Any}(name, accessor, path)
end

function get_property(item::FieldExtractor, model::AbstractModel)
    return item.field_accessor(model)
end

function register_field!(collector::Collector, name::String, path::Vector{Symbol})
    accessor = function(model)
        result = model
        for field in path
            result = getproperty(result, field)
        end
        return result
    end
    extractor = field_extractor(name, accessor, path)
    register_item!(collector, extractor)
end
register_field!(collector, name, path::String) = register_field!(collector, name, Symbol.(split(path, ".")))
register_field!(collector, name::Symbol, path::Vector{Symbol}) = register_field!(collector, string(name), path)

end