using PkgBenchmark
using Dates
using LCS: LCS

output = get(ARGS, 1, "benchmark.md")

results = benchmarkpkg(LCS, BenchmarkConfig(; juliacmd=Base.julia_cmd(), env=Dict("JULIA_NUM_THREADS" => "auto")))

export_markdown(output, results)
