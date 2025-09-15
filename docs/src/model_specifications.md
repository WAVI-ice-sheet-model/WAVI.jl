# WAVI Model Specifications

## Overview

Drives the execution
Describe setup: model/outputs/time -> simulation -> specifications

A late addition to WAVI is the idea of "Specs", which guide the setup of the model itself. These can be seen to operate underneath `Simulation`, which brings together `Model` and `Clock` to apply a lifecycle for WAVI as an full model. 

Semantically, we can think of the components this way: 

* `Model` brings together all the elements of compute and data into a singular operable entity
* `Spec` determines how the model actually operates at the computational level
* `Simulation` determines how the model will be used to produce results, incorporating and controlling it's development over time irrespective of the underlying `Spec`

## Types

### Model Specifications

There are three types of specifications currently available in WAVI: 

* `BasicSpec` - this is the original single-core, single-thread implementation. You will run into limits with bigger domains.
* `ThreadedSpec` - this is an implementation of schwarz decomposition for solving subdomains within individual threads. It won't get you around the memory limits, but it will speed up compute.
* `MPISpec` - this is the currently experimental fully-distributed specification for WAVI. The model domain is split apart and computed individually.

### Implementation Considerations

Depending on specifications, some considerations need to be taken into account:

* Examples in the docs
* How you use `Model` interface changes in a distributed setting, to a fully functional interface

## Setting up

### Driver files

### MPI specific execution

MPI.jl / MPIPreferences


`mpiexecjl --project -np <num> julia <path to driver> <driver args...>`