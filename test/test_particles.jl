using Test

using LCS
using Random
using ReferenceTests

import KernelAbstractions as KA
import LCS.Flows as F
import LCS.Particles as P
import Topologies as T

REFERENCE = joinpath(@__DIR__, "references/particles.jld2")

backend = KA.CPU()
topo = T.Topology()

grid = LCS.Grid((16, 16, 16), 1)

fparams = F.FlowParams(1.0, 1.0, 1.0)
pparams = P.ParticleParams(9.8, 1.0)
population = P.Population(10, 0.0, 0.0)
nvalid = P.nvalid_l(population, topo)

dt = 0.01

cell = P.Cell(4)
stat = P.Stat(P.RDF(1E-3, 4, 1E-2, 1))

fields = F.Fields(backend, grid, topo)
props = P.Properties(backend, population, topo)
cibuf = P.CellIndexBuffer(backend, grid, population, cell, topo)

@testset "init" begin
    @testset "helper" begin
        nvalid_l_all = [3, 4, 5, 6]
        @test P.id_offset(nvalid_l_all, 0) == 0
        @test P.id_offset(nvalid_l_all, 1) == 3
        @test P.id_offset(nvalid_l_all, 2) == 7
        @test P.id_offset(nvalid_l_all, 3) == 12
    end

    @testset "id" begin
        let props = deepcopy(props)
            P.init_id!(P.GenerateId(); props.ids, nvalid, backend, topo)
            @test_reference REFERENCE Dict("particles/init/id/ids" => view(props.ids, 1:nvalid))
        end
    end

    @testset rng = Random.Xoshiro(0) "position" begin
        let props = deepcopy(props)
            P.init_position!(P.RandomPosition(); props.xss, nvalid, backend, topo)
            @test_reference REFERENCE Dict(
                "particles/init/position/x" => view(props.xss[1], 1:nvalid),
                "particles/init/position/y" => view(props.xss[2], 1:nvalid),
                "particles/init/position/z" => view(props.xss[3], 1:nvalid),
            )
        end
    end

    @testset "velocity" begin
        let props = deepcopy(props)
            P.init_velocity!(P.RestVelocity(); props.uss, nvalid, backend)
            @test_reference REFERENCE Dict(
                "particles/init/velocity/u" => view(props.uss[1], 1:nvalid),
                "particles/init/velocity/v" => view(props.uss[2], 1:nvalid),
                "particles/init/velocity/w" => view(props.uss[3], 1:nvalid),
            )
        end
    end

    @testset "size" begin
        let props = deepcopy(props)
            P.init_size!(P.ConstSize(1.0); props.diams, nvalid, fparams, backend)
            @test_reference REFERENCE Dict("particles/init/size/diams" => view(props.diams, 1:nvalid))
        end
    end
end

@testset rng = Random.Xoshiro(0) "makeindex" begin
    let props = deepcopy(props), cibuf = deepcopy(cibuf)
        P.init_id!(P.GenerateId(); props.ids, nvalid, backend, topo)
        P.init_position!(P.RandomPosition(); props.xss, nvalid, backend, topo)

        P.makeindex!(; props.xss, props, cibuf, nvalid, grid, cell, backend, topo)
        @test sum(cibuf.stops - cibuf.starts .+ 1) == nvalid

        @test_reference REFERENCE Dict(
            "particles/makeindex/hashes" => view(cibuf.hashes, 1:nvalid),
            "particles/makeindex/starts" => cibuf.starts,
            "particles/makeindex/stops" => cibuf.stops,
        )
    end
end

@testset rng = Random.Xoshiro(0) "motion!" begin
    @test P.coeff(P.LinearDrag(), 1.0) == 1.0
    @test P.coeff(P.NonlinearDrag(), 1.0) == 1.15

    @testset "linear-drag" begin
        let fields = deepcopy(fields), props = deepcopy(props)
            Random.rand!.(fields.Us)
            Random.rand!.(props.xss)
            Random.rand!.(props.uss)

            drag = P.LinearDrag()
            P.motion!(; fields.Us, props.xss, props.uss, props, nvalid, grid, fparams, pparams, drag, backend, topo)
            @test_reference REFERENCE Dict(
                "particles/motion/linear-drag/du" => view(props.dudtss[1], 1:nvalid),
                "particles/motion/linear-drag/dv" => view(props.dudtss[2], 1:nvalid),
                "particles/motion/linear-drag/dw" => view(props.dudtss[3], 1:nvalid),
            )
        end
    end

    @testset "nonlinear-drag" begin
        let fields = deepcopy(fields), props = deepcopy(props)
            Random.rand!.(fields.Us)
            Random.rand!.(props.xss)
            Random.rand!.(props.uss)

            drag = P.NonlinearDrag()
            P.motion!(; fields.Us, props.xss, props.uss, props, nvalid, grid, fparams, pparams, drag, backend, topo)
            @test_reference REFERENCE Dict(
                "particles/motion/nonlinear-drag/du" => view(props.dudtss[1], 1:nvalid),
                "particles/motion/nonlinear-drag/dv" => view(props.dudtss[2], 1:nvalid),
                "particles/motion/nonlinear-drag/dw" => view(props.dudtss[3], 1:nvalid),
            )
        end
    end
end

@testset rng = Random.Xoshiro(0) "update!" begin
    let props = deepcopy(props)
        Random.rand!.(props.xss)
        Random.rand!.(props.uss)
        Random.rand!.(props.dudtss)

        P.update!(; xss=props.xss2, xss0=props.xss, uss=props.uss2, uss0=props.uss, props, dt, nvalid, backend)
        @test_reference REFERENCE Dict(
            "particles/update/x2" => view(props.xss2[1], 1:nvalid),
            "particles/update/y2" => view(props.xss2[2], 1:nvalid),
            "particles/update/z2" => view(props.xss2[3], 1:nvalid),
            "particles/update/u2" => view(props.uss2[1], 1:nvalid),
            "particles/update/v2" => view(props.uss2[2], 1:nvalid),
            "particles/update/w2" => view(props.uss2[3], 1:nvalid),
        )
    end
end

@testset "timestep" begin
    let props = deepcopy(props)
        Random.rand!(props.diams)

        ts = LCS.TimeStep(0.3, 0.01)
        dt = P.timestep(props.diams, ts, fparams, pparams, topo)
        @test_reference REFERENCE Dict("particles/timestep/dt" => dt)
    end
end

@testset rng = Random.Xoshiro(0) "rdf" begin
    let props = deepcopy(props), cibuf = deepcopy(cibuf)
        Random.rand!.(props.xss)

        P.makeindex!(; xss=props.xss, props, cibuf, nvalid, grid, cell, backend, topo)
        rdf = P.rdf_assume_ci(; diam=1.0, props.xss, props, cibuf, nvalid, grid, population, cell, topo, stat, backend)
        @test_reference REFERENCE Dict(
            "particles/rdf/gr_contact" => rdf.gr_contact,
            "particles/rdf/edges" => rdf.edges,
            "particles/rdf/gr" => rdf.gr,
            "particles/rdf/npairs" => rdf.npairs,
        )
    end
end

@testset rng = Random.Xoshiro(0) "density" begin
    let props = deepcopy(props), cibuf = deepcopy(cibuf)
        Random.rand!.(props.xss)

        densities = P.density(; props.xss, props, cibuf, nvalid, grid, cell, backend, topo)
        @test_reference REFERENCE Dict("particles/density/densities" => densities)
    end
end
