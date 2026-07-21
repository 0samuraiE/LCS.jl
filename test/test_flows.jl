using Test

using LCS: LCS, Utils
using Parallel
using Random
using ReferenceTests

import LCS.Flows as F
import Topologies as T
import KernelAbstractions as KA

backend = KA.CPU()
grid = LCS.Grid((16, 16, 16), 3)
ts = LCS.TimeStep(0.3, 0.01)
params = F.FlowParams(1, 1, 1)

t = 0.0
dt = 0.01

topo = T.Topology()

fields = F.Fields(backend, grid, topo)
scratches = F.Scratches(backend, grid, topo)
halobuf = T.HaloBuffer(backend, topo, LCS.FP, grid.dims, grid.halo_size, LCS.get_max_comm_arrays())

REFERENCE = joinpath(@__DIR__, "references/flows.jld2")

@testset "flows" begin
    @testset "init" begin
        @testset "ideal-flow" begin
            let fields = deepcopy(fields)
                init = F.IdealFlow()
                F.init!(init; fields.Us, fields, grid, backend, topo)
                @test_reference REFERENCE Dict(
                    "flows/init/ideal-flow/U" => LCS.unhalo(fields.Us[1], grid),
                    "flows/init/ideal-flow/V" => LCS.unhalo(fields.Us[2], grid),
                    "flows/init/ideal-flow/W" => LCS.unhalo(fields.Us[3], grid),
                )
            end
        end

        @testset "random-flow" begin
            let fields = deepcopy(fields)
                init = F.RandomFlow(0)
                F.init!(init; fields.Us, fields, grid, topo)
                @test_reference REFERENCE Dict(
                    "flows/init/random-flow/U" => LCS.unhalo(fields.Us[1], grid),
                    "flows/init/random-flow/V" => LCS.unhalo(fields.Us[2], grid),
                    "flows/init/random-flow/W" => LCS.unhalo(fields.Us[3], grid),
                )
            end
        end
    end

    @testset rng = Random.Xoshiro(0) "rhs" begin
        let fields = deepcopy(fields)
            Random.rand!.(fields.Us)

            F.rhs!(; fields.Us, fields, scratches, grid, params, backend, topo)
            @test_reference REFERENCE Dict(
                "flows/rhs/dU" => LCS.unhalo(fields.dUdts[1], grid),
                "flows/rhs/dV" => LCS.unhalo(fields.dUdts[2], grid),
                "flows/rhs/dW" => LCS.unhalo(fields.dUdts[3], grid),
            )
        end
    end

    @testset rng = Random.Xoshiro(0) "update" begin
        let fields = deepcopy(fields)
            Random.rand!.(fields.Us)
            Random.rand!.(fields.dUdts)

            F.update!(; Us=fields.Us2, Us0=fields.Us, fields, dt, grid, backend, topo)
            @test_reference REFERENCE Dict(
                "flows/update/U2" => LCS.unhalo(fields.Us2[1], grid),
                "flows/update/V2" => LCS.unhalo(fields.Us2[2], grid),
                "flows/update/W2" => LCS.unhalo(fields.Us2[3], grid),
            )
        end
    end

    @testset "couple" begin
        @testset rng = Random.Xoshiro(0) for time_blocking in (true, false), overlap in (true, false)
            let fields = deepcopy(fields)
                Random.rand!.(fields.Us)
                Random.rand!(fields.P)

                hsmac = F.HSMAC(time_blocking, overlap, Parallel.RedBlackFast(), 1.7, 1E-3, 6, 6)
                F.couple!(hsmac; fields.Us, fields, scratches, halobuf, dt, grid, backend, topo)

                case = "hsmac"
                case *= time_blocking ? "-time-blocking" : ""
                case *= overlap ? "-overlap" : ""
                @test_reference REFERENCE Dict(
                    "flows/couple/$(case)/U" => LCS.unhalo(fields.Us[1], grid),
                    "flows/couple/$(case)/V" => LCS.unhalo(fields.Us[2], grid),
                    "flows/couple/$(case)/W" => LCS.unhalo(fields.Us[3], grid),
                    "flows/couple/$(case)/P" => LCS.unhalo(fields.P, grid),
                )
            end
        end
    end

    @testset "forcing" begin
        @testset "linear" begin
            @testset rng = Random.Xoshiro(0) "constant-power" begin
                let fields = deepcopy(fields)
                    Random.rand!.(fields.Us)

                    forcing = F.ConstantPowerLF(1.0)
                    F.forcing!(forcing; fields.Us, U2=scratches[1], halobuf, t, dt, grid, backend, topo)
                    @test_reference REFERENCE Dict(
                        "flows/forcing/constant-power-lf/U" => LCS.unhalo(fields.Us[1], grid),
                        "flows/forcing/constant-power-lf/V" => LCS.unhalo(fields.Us[2], grid),
                        "flows/forcing/constant-power-lf/W" => LCS.unhalo(fields.Us[3], grid),
                    )
                end
            end
        end

        @testset "reduced-communication" begin
            @testset "boxmean" begin
                let
                    dims = (8, 8, 8)
                    coarse_dims = (2, 2, 2)
                    filter_dims = div.(dims, coarse_dims)
                    halo_size = 3

                    A = zeros(Float64, dims .+ 2 * halo_size)
                    A_coarse_ref = reshape(1:prod(coarse_dims), coarse_dims...)
                    @views A[map(n -> (1:n) .+ halo_size, dims)...] .= repeat(A_coarse_ref; inner=filter_dims)
                    A_coarse = zeros(eltype(A), coarse_dims...)
                    F.boxmean_l!(backend, A_coarse, A, halo_size, filter_dims)
                    @test A_coarse == A_coarse_ref
                end
            end

            @testset "helpers" begin
                let
                    forcing = F.EnergyPreserveRCF(1.0, 1, 2, 4)

                    mocktopo = T.Mock.mocktopology(; linear_rank=1, proc_dims=(4, 2, 1))
                    @test F.coarse_dims_g(forcing) == (4, 4, 4)
                    @test F.coarse_dims_l(forcing, mocktopo) == (1, 2, 4)
                    @test F.filter_size(forcing, grid) == 4
                    @test F.phase_shift(forcing, grid) == 2
                end
            end

            @testset rng = Random.Xoshiro(0) "energy-preserve" begin
                let fields = deepcopy(fields)
                    Random.rand!.(fields.Us)

                    forcing = F.EnergyPreserveRCF(1.0, 1, 2, 4)
                    forcingbuf = F.RCFBuffer(forcing, backend, topo)
                    Power = scratches[1]
                    F.forcing!(forcing; fields.Us, Power, forcingbuf, dt, grid, backend, topo)
                    @test_reference REFERENCE Dict(
                        "flows/forcing/energy-preserve-rcf/U" => LCS.unhalo(fields.Us[1], grid),
                        "flows/forcing/energy-preserve-rcf/V" => LCS.unhalo(fields.Us[2], grid),
                        "flows/forcing/energy-preserve-rcf/W" => LCS.unhalo(fields.Us[3], grid),
                    )
                end
            end
        end
    end

    @testset rng = Random.Xoshiro(0) "stat" begin
        let fields = deepcopy(fields)
            Random.rand!.(fields.Us)

            ilbuf = F.IntegralLengthBuffer(backend, grid, topo)
            statbuf = F.StatBuffer(ilbuf)
            stat = F.stat(; fields.Us, scratches, halobuf, statbuf, grid, params, backend, topo)
            @test_reference REFERENCE Dict("flows/stat/stat" => stat)

            Q = F.Q(; fields.Us, scratches, halobuf, grid, backend, topo)
            @test_reference REFERENCE Dict("flows/stat/Q" => Q)
        end
    end

    @testset "timestep" begin
        let fields = deepcopy(fields)
            Random.rand!.(fields.Us)

            dt = F.timestep(; fields.Us, scratches, halobuf, grid, ts, params, backend, topo)
            @test_reference REFERENCE Dict("flows/timestep/dt" => dt)
        end
    end
end
