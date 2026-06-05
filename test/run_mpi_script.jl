# Shared helper: spawn an MPI test script via MPI.jl's mpiexec() (see MPI.jl docs "Writing MPI tests").
function run_mpi_script(script::String; nprocs::Int = 2)
    script_path = abspath(joinpath(@__DIR__, "mpi_tests", script))
    project = Base.active_project()
    cmd = `$(mpiexec()) -n $nprocs $(Base.julia_cmd()) --project=$project $script_path`
    return run(ignorestatus(cmd)), cmd
end
