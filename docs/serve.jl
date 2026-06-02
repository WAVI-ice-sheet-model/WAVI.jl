#!/usr/bin/env julia

using Pkg

# Use docs project
Pkg.activate(@__DIR__)

# Add WAVI source path
Pkg.develop(path = joinpath(@__DIR__, "..")) # This adds local WAVI to docs/Project.toml
#insert!(LOAD_PATH, 1, joinpath(@__DIR__, "..")) # Add local WAVI git clone - prioritise over a normal install

Pkg.instantiate()

if "--live" in ARGS
    try
        @info "Starting LiveServer with auto-reload..."
        using LiveServer
        # Monitor WAVI source code for changes as well
        servedocs(include_dirs=["src/"])
    catch e
        @error "LiveServer failed to start: $e"
    end
else
    try
        include("make.jl")
        @info "Documentation building completed."
    catch e
        @error "Error building documentation."
    end
end

