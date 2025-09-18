using Revise
using WAVI
using LiveServer

Revise.revise()
include("make.jl")
servedocs(
    foldername=".",
    include_dirs=["./src/"]
)
