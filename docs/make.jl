using Documenter
using LCS
using PolySerde
using Topologies
using Parallel
using Offsets
using PPrint
using ReferenceTests
using Schema
using FourierTools

makedocs(;
    modules=[
        LCS,
        LCS.Flows,
        LCS.Particles,
        LCS.LCSIO,
        FourierTools,
        Offsets,
        Parallel,
        PolySerde,
        PPrint,
        ReferenceTests,
        Schema,
        Topologies,
    ],
    sitename="LCS.jl",
    authors="Taketo Tominaga",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://0samuraiE.github.io/LCS.jl",
        assets=String[],
        size_threshold=nothing,
    ),
    pages=[
        "Home" => "index.md",
        "Guides" => ["Configuration & Environment" => "guide/configuration.md"],
        "API Reference" => [
            "LCS Core" => "api/core.md",
            "Flows" => "api/flows.md",
            "Particles" => "api/particles.md",
            "LCSIO" => "api/lcsio.md",
        ],
        "Libraries" => [
            "FourierTools" => "lib/fouriertools.md",
            "Offsets" => "lib/offsets.md",
            "Parallel" => "lib/parallel.md",
            "PolySerde" => "lib/polyserde.md",
            "PPrint" => "lib/pprint.md",
            "ReferenceTests" => "lib/referencetests.md",
            "Schema" => "lib/schema.md",
            "Topologies" => "lib/topologies.md",
        ],
    ],
    checkdocs=:none,
    remotes=nothing,
    doctest=false,
    warnonly=[:docs_block, :cross_references, :missing_docs],
)
