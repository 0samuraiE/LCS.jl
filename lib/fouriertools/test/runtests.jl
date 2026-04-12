using Test

using FFTW
using StaticArrays

import FourierTools as FT

@test FT.realprod(1 + 2im, 3 + 4im) == real((1 + 2im) * (3 + 4im))

@testset "fftindex" begin
    # Test positive indices
    @test FT.fftindex(0, 8) == 1
    @test FT.fftindex(3, 8) == 4
    @test FT.fftindex(7, 8) == 8

    # Test negative indices
    @test FT.fftindex(-1, 8) == 8
    @test FT.fftindex(-3, 8) == 6
    @test FT.fftindex(-7, 8) == 2
end

@testset "fftcoeff" begin
    @testset "1D" begin
        x = randn(8)
        F_fft = fft(x)
        F_rfft = rfft(x)

        for k in -3:4
            @test FT.fftcoeff(F_fft, (k,)) ≈ FT.rfftcoeff(F_rfft, (k,))
        end
    end

    @testset "2D" begin
        x = randn(8, 6)
        F_fft = fft(x)
        F_rfft = rfft(x)

        for k2 in -3:2, k1 in -3:4
            @test FT.fftcoeff(F_fft, (k1, k2)) ≈ FT.rfftcoeff(F_rfft, (k1, k2))
        end
    end

    @testset "3D" begin
        x = randn(4, 6, 8)
        F_fft = fft(x)
        F_rfft = rfft(x)

        for k3 in -3:4, k2 in -2:3, k1 in -1:2
            @test FT.fftcoeff(F_fft, (k1, k2, k3)) ≈ FT.rfftcoeff(F_rfft, (k1, k2, k3))
        end
    end

    @testset "conjugate symmetry" begin
        x = randn(8, 8)
        F_rfft = rfft(x)

        @test FT.rfftcoeff(F_rfft, (-1, 2)) == conj(FT.rfftcoeff(F_rfft, (1, -2)))
        @test FT.rfftcoeff(F_rfft, (-2, -1)) == conj(FT.rfftcoeff(F_rfft, (2, 1)))
    end
end

@testset "canon" begin
    @test FT.canon((1, 2, 3)) == (1, 2, 3)
    @test FT.canon((-1, -2, -3)) == (1, 2, 3)
    @test FT.canon((0, 1, 2)) == (0, 1, 2)
    @test FT.canon((0, -1, -2)) == (0, 1, 2)
    @test FT.canon((0, 0, 1)) == (0, 0, 1)
    @test FT.canon((0, 0, -1)) == (0, 0, 1)
end

@testset "kset" begin
    @test Set(FT.kset(1)) == Set([#
        (0, 1, 1),
        (0, 0, 1),
        (1, 0, -1),
        (0, 1, -1),
        (0, 1, 0),
        (1, -1, 0),
        (1, 0, 0),
        (1, 1, 0),
        (1, 0, 1),
    ])
end

@testset "kcount" begin
    @test FT.kcount(1) == 18
    @test FT.kcount(2) == 62
end

@testset "enespe" begin
    U_hat = fill(one(ComplexF64) * 16^3, 9, 16, 16)
    @test FT.enespe(U_hat, U_hat, U_hat, 1) == FT.kcount(1) * 3 / 2
    @test FT.enespe(U_hat, U_hat, U_hat, 2) == FT.kcount(2) * 3 / 2
end

@testset "sidft3" begin
    n = 8
    kmin = 2
    kmax = 3

    A = rand(n, n, n)
    p = FFTW.plan_rfft(A)
    A_hat = rfft(A)

    B_hat = FT.crop(A_hat, kmax)

    A2 = fill!(similar(A), 0)
    for k in axes(A2, 3), j in axes(A2, 2), i in axes(A2, 1)
        exp1k = cispi(2 * (k - 1) / n)
        exp1j = cispi(2 * (j - 1) / n)
        exp1i = cispi(2 * (i - 1) / n)

        a = FT.sirdft(Val(kmin), Val(kmax), B_hat, exp1i, exp1j, exp1k)

        A2[i, j, k] = a
    end

    A3_hat = fill!(similar(A_hat), 0)
    for kz in (-kmax):kmax, ky in (-kmax):kmax, kx in 0:kmax
        k = (kx, ky, kz)
        if kmin - 0.5 <= sqrt(sum(abs2, k)) < kmax + 0.5
            A3_hat[FT.fftindex.((kx, ky, kz), size(A3_hat))...] = A_hat[FT.fftindex.((kx, ky, kz), size(A_hat))...]
        end
    end

    A3 = irfft(A3_hat, size(A, 1))
    @test A2 ≈ A3
end
