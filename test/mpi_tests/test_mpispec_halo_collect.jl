# Unit tests for MPISpec halo exchange and global field collection (run under MPI with 2 ranks).
using Test
using MPI
using WAVI

# Halo exchange, thickness sync, and global collect on a 2-rank MPISpec split.
@testset "MPISpec halo/global collect unit tests" begin
    MPI.Initialized() || MPI.Init()
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)

    @test nprocs == 2

    """
        build_model(; nx, ny, halo, pou)

    Build a small MPISpec model split across MPI ranks (2×1 layout) for halo tests.
    """
    function build_model(; nx = 12, ny = 8, halo = 1, pou = false)
        grid = Grid(nx = nx, ny = ny, nσ = 3, dx = 1.0, dy = 1.0)
        spec = MPISpec(nprocs, 1, halo, grid; pou = pou, niterations = 2)
        bed = -500.0 .* ones(nx, ny)
        model = Model(
            grid = grid,
            bed_elevation = bed,
            params = Params(accumulation_rate = 0.1),
            solver_params = SolverParams(maxiter_picard = 1),
            initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(nx, ny)),
            spec = spec,
        )
        return model
    end

    # After exchange, the left halo should hold the left neighbour's core thickness (rank 1 on a 2×1 split).
    @testset "halo_exchange! updates rank-neighbour h halos" begin
        model = build_model(halo = 1, pou = false)
        th, rh, bh, lh = WAVI.Specs.get_halos(model.spec)
        h = model.fields.gh.h
        h .= -99.0
        h[(1 + lh):(end - rh), (1 + th):(end - bh)] .= rank + 1.0

        WAVI.Specs.halo_exchange!(model; fields = [:h])

        if model.spec.left > -1 && lh > 0
            @test all(h[1:lh, :] .== (model.spec.left + 1.0))
        end
    end

    # With partition of unity, thickness sync must not overwrite assembled u and v halos.
    @testset "mpi_sync_halos_after_thickness! preserves u/v when pou=true" begin
        model = build_model(halo = 2, pou = true)
        u_before = copy(model.fields.gu.u)
        v_before = copy(model.fields.gv.v)

        WAVI.Specs.mpi_sync_halos_after_thickness!(model)

        @test model.fields.gu.u == u_before
        @test model.fields.gv.v == v_before
    end

    # Global h on rank 0 should contain core values from both MPI ranks after collect.
    @testset "collect_mpi_field! assembles both rank patches" begin
        model = build_model(halo = 1, pou = false)
        th, rh, bh, lh = WAVI.Specs.get_halos(model.spec)
        model.fields.gh.h .= -1.0
        model.fields.gh.h[(1 + lh):(end - rh), (1 + th):(end - bh)] .= rank + 10.0

        gathered = WAVI.Specs.collect_mpi_field!(model, [:global_fields, :gh, :h])
        if rank == 0
            vals = Set(vec(gathered))
            @test 10.0 in vals
            @test 11.0 in vals
        end
    end

    # Nested ThreadedSpec local solve should complete a velocity update with finite fields.
    @testset "MPISpec with nested ThreadedSpec local solve" begin
        grid = Grid(nx = 12, ny = 8, nσ = 3, dx = 1.0, dy = 1.0)
        local_spec = ThreadedSpec(ngridsx = 2, ngridsy = 1, overlap = 1, niterations = 1)
        spec = MPISpec(nprocs, 1, 2, grid; pou = true, niterations = 1, local_spec = local_spec)
        model = Model(
            grid = grid,
            bed_elevation = -500.0 .* ones(grid.nx, grid.ny),
            params = Params(accumulation_rate = 0.1),
            solver_params = SolverParams(maxiter_picard = 1),
            initial_conditions = InitialConditions(initial_thickness = 100.0 .* ones(grid.nx, grid.ny)),
            spec = spec,
        )
        @test model.spec.local_spec isa ThreadedSpec
        WAVI.Processes.update_model_velocities!(model)
        @test all(isfinite.(model.fields.gu.u))
        @test all(isfinite.(model.fields.gv.v))
    end

    # 3D halo exchange: each vertical level of a 3D field should be updated by the exchange.
    @testset "halo_exchange! updates halos of 3D fields (all levels)" begin
        model = build_model(halo = 1, pou = false)
        th, rh, bh, lh = WAVI.Specs.get_halos(model.spec)

        # g3d.η is a genuine (nx × ny × nσ) field.
        # Fill the whole field with a known background value, then write a unique
        # value into each rank's core for each vertical level so we can verify
        # that the halo exchange copies the correct data from each neighbour.
        η = model.fields.g3d.η
        nσ = size(η, 3)
        η .= -99.0
        for k in 1:nσ
            η[(1 + lh):(end - rh), (1 + th):(end - bh), k] .= (rank + 1.0) * (k + 10.0)
        end

        WAVI.Specs.halo_exchange!(model; fields = [:η])

        # Rank that has a left neighbour: its left halo should hold that neighbour's core value.
        if model.spec.left > -1 && lh > 0
            for k in 1:nσ
                expected = (model.spec.left + 1.0) * (k + 10.0)
                @test all(η[1:lh, :, k] .≈ expected)
            end
        end
    end

    # 3D global collect: collect_mpi_field! must assemble all vertical levels onto root.
    @testset "collect_mpi_field! assembles both rank patches for a 3D field" begin
        model = build_model(halo = 1, pou = false)
        th, rh, bh, lh = WAVI.Specs.get_halos(model.spec)

        η = model.fields.g3d.η
        nσ = size(η, 3)
        η .= -1.0
        for k in 1:nσ
            η[(1 + lh):(end - rh), (1 + th):(end - bh), k] .= (rank + 1.0) * (k + 10.0)
        end

        gathered = WAVI.Specs.collect_mpi_field!(model, [:global_fields, :g3d, :η])
        if rank == 0
            # Every (rank, level) sentinel must appear somewhere in the global slice.
            for k in 1:nσ
                @test any(gathered[:, :, k] .≈ (1.0) * (k + 10.0))  # rank-0 contribution
                @test any(gathered[:, :, k] .≈ (2.0) * (k + 10.0))  # rank-1 contribution
            end
        end
    end

end

MPI.Finalize()
