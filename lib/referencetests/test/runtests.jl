using Test

using JLD2
using ReferenceTests

@testset "ReferenceTests" begin
    @testset "test_reference" begin
        let
            file = tempname()
            try
                jldopen(file, "w") do f
                    f["a"] = 1
                    f["b"] = [1, 2, 3]
                    f["c"] = "test"
                end

                @test ReferenceTests.test_reference(file, Dict("a" => 1))
                @test ReferenceTests.test_reference(file, Dict("a" => 1, "b" => [1, 2, 3]))
                @test ReferenceTests.test_reference(file, Dict("a" => 1, "b" => [1, 2, 3], "c" => "test"))

                @test !ReferenceTests.test_reference(file, Dict("a" => 2))
                @test !ReferenceTests.test_reference(file, Dict("b" => [1, 2]))

                @test_throws ArgumentError ReferenceTests.test_reference(file, Dict("nonexistent" => 1))

                @test ReferenceTests.test_reference(file, Dict("a" => 1.0); by=isapprox)

                @test_logs (:warn, r"mismatch for key \"a\"") !ReferenceTests.test_reference(file, Dict("a" => 2))

                @test_logs (:warn, r"mismatch for key \"b\"") !ReferenceTests.test_reference(
                    file, Dict("b" => [4, 5, 6])
                )

                @test_logs (:warn, r"mismatch for key \"b\"") !ReferenceTests.test_reference(file, Dict("b" => [1, 2]))

                @test_logs(
                    (:warn, r"mismatch for key"),
                    (:warn, r"mismatch for key"),
                    !ReferenceTests.test_reference(file, Dict("a" => 2, "b" => [4, 5, 6]))
                )
            finally
                isfile(file) && rm(file)
            end
        end
    end

    @testset "update_reference" begin
        let
            file = tempname()
            try
                jldopen(file, "w") do f
                    f["a"] = 1
                    f["b"] = 2
                end

                ReferenceTests.update_reference(file, Dict("a" => 10, "c" => 3))

                jldopen(file, "r") do f
                    @test f["a"] == 10
                    @test f["b"] == 2
                    @test f["c"] == 3
                end
            finally
                isfile(file) && rm(file)
            end
        end
    end

    @testset "jldrepack" begin
        let
            file = tempname()
            try
                jldopen(file, "w") do f
                    f["a"] = 1
                    f["b/c"] = 2
                    delete!(f, "a")
                    f["a"] = 3
                    delete!(f, "b")
                    f["b/c"] = 4
                end

                size1 = filesize(file)

                ReferenceTests.jldrepack(file)

                @test isfile(file)

                jldopen(file, "r") do f
                    @test f["a"] == 3
                    @test f["b/c"] == 4
                end

                size2 = filesize(file)
                @test size2 < size1
            finally
                isfile(file) && rm(file)
            end
        end
    end

    @testset "@test_reference" begin
        let
            file = tempname()
            try
                jldopen(file, "w") do f
                    g = JLD2.Group(f, "group")
                    g["x"] = 1
                    g["y"] = 2
                end

                withenv(ReferenceTests.ENV_JULIA_REFERENCETESTS_UPDATE => "false") do
                    @test_reference file Dict("group/x" => 1, "group/y" => 2)
                end

                withenv(ReferenceTests.ENV_JULIA_REFERENCETESTS_UPDATE => "true") do
                    @test_logs (:info, r"Updated reference file") (@test_reference file Dict("group/x" => 100))
                end

                jldopen(file, "r") do f
                    @test f["group/x"] == 100
                    @test f["group/y"] == 2
                end
            finally
                isfile(file) && rm(file)
            end
        end
    end
end
