# Running on HPC Clusters

WAVI.jl is designed to be highly scalable and can be run efficiently on High-Performance Computing (HPC) clusters. Below are examples of how to submit jobs using the SLURM workload manager for serial, multi-threaded, and MPI configurations.

These examples are based on the `MISMIP_PLUS` driver provided in the `example_drivers/MISMIP_PLUS` directory for running on the British Antarctic Survey's HPC.

## Serial Execution

For smaller test runs or non-parallel configurations, you can run WAVI.jl in serial. 

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.serial.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```

## Multi-threaded Execution

To take advantage of shared-memory parallelism on a single node, you can run Julia with multiple threads. Set the `--cpus-per-task` in your SLURM script and pass this value to Julia using the `-t` flag.

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.threaded.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```

## MPI Execution

For large-scale simulations spanning multiple nodes or making use of distributed memory, WAVI.jl supports MPI. You will need to load the appropriate MPI module (e.g., `mpich`) and use `mpiexecjl` to launch the Julia script. 

*Note: Ensure that you have `MPI.jl` correctly configured on your cluster to work with the system's MPI binaries.*

```@eval
using Markdown
using WAVI
path = joinpath(pkgdir(WAVI), "example_drivers", "MISMIP_PLUS", "submit-driver.mpi.sh")
content = read(path, String)
Markdown.MD(Markdown.Code("bash", content))
```
