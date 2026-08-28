dir = @__DIR__
files = readdir(dir)
excludes = Set(["common.jl", "runall.jl"])

files = filter(files) do f
    endswith(f, ".jl") && !(f in excludes)
end

for f in files
    m = gensym()
    ex = :(
        module $m
        include(joinpath($dir, $f))
        end
    )
    Core.eval(Main, ex)
end
