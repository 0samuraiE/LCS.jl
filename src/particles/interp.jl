#=
  3-|-----------|-----------|-----------|
    |           |           |           |
    |         (1,3)      (2,3)        (3,3)
   ==>         ==>        ==>          ==>
    |           |           |           |
    |           |           |           |
  2-|-----------|-----------|-----------|
    |           |           |           |
    |         (1,2)       (2,2)       (2,3)
   ==>         ==>         ==>         ==>
    |           |           |           |
    |           |           |           |
  1-|-----------|(1.4,0.9)--|-----------|
    |           |    *      |           |
    |        @(1,1)  .    (2,1)       (3,1)
   ==>         ==>....r2   ==>         ==>
    |           |  r1       |           |
    |           |           |           |
  0-|-----------|-----------|-----------|
    0           1           2           3

(x, y) = (1.4, 0.9)
(i, j) = (1, 1) == ([x], [y + 0.5])
(r1, r2) = (0.4, 0.4) == (x - i, y - j + 0.5)

  3-|-----------|-----------|-----------|
    |           |           |           |
    |         (1,3)      (2,3)        (3,3)
   ==>         ==>        ==>          ==>
    |           |           |           |
    |           |           |           |
  2-|-----------|-----------|-----------|
    |           |           |           |
    |         (1,2)       (2,2)       (2,3)
   ==>         ==>         ==>         ==>
    |           |(1.4,1.1)  |           |
    |           |    *      |           |
  1-|-----------|----.------|-----------|
    |           |    .r2    |           |
    |        @(1,1)  .    (2,1)       (3,1)
   ==>         ==>....     ==>         ==>
    |           |  r1       |           |
    |           |           |           |
  0-|-----------|-----------|-----------|
    0          1          2          3

(x, y) = (1.4, 1.1)
(i, j) = (1, 1) == ([x], [y + 0.5])
(r1, r2) = (0.3, 0.6) == (x - i, y - j + 0.5)
=#

struct Mapping{Int<:Integer,Float<:AbstractFloat}
    i  :: Int
    j  :: Int
    k  :: Int
    r1 :: Float
    r2 :: Float
    r3 :: Float
end

Base.@propagate_inbounds function interp(U::Field, halo_size::Integer, m::Mapping)
    (; i, j, k, r1, r2, r3) = m

    @offsetviews halo_size begin
        c0 = U[i, j, k]
        c1 = U[i + 1, j, k]
        c2 = U[i + 1, j + 1, k]
        c3 = U[i, j + 1, k]
        c4 = U[i, j, k + 1]
        c5 = U[i + 1, j, k + 1]
        c6 = U[i + 1, j + 1, k + 1]
        c7 = U[i, j + 1, k + 1]
    end

    _interp_trilinear(c0, c1, c2, c3, c4, c5, c6, c7, r1, r2, r3)
end

Base.@propagate_inbounds interp(U::Field, grid::LCS.Grid, m::Mapping) = interp(U, grid.halo_size, m)

@inline function _interp_trilinear(c0, c1, c2, c3, c4, c5, c6, c7, r1, r2, r3)
    (1 - r1) * (1 - r2) * (1 - r3) * c0 +
    r1 * (1 - r2) * (1 - r3) * c1 +
    r1 * r2 * (1 - r3) * c2 +
    (1 - r1) * r2 * (1 - r3) * c3 +
    (1 - r1) * (1 - r2) * r3 * c4 +
    r1 * (1 - r2) * r3 * c5 +
    r1 * r2 * r3 * c6 +
    (1 - r1) * r2 * r3 * c7
end

@inline function xmapping(x, y, z, dx, dy, dz)
    i = Utils.unsafe_floor(x / dx)
    j = Utils.unsafe_floor(y / dy + 0.5)
    k = Utils.unsafe_floor(z / dz + 0.5)

    r1 = x / dx - i
    r2 = y / dy - j + 0.5
    r3 = z / dz - k + 0.5

    Mapping(i, j, k, r1, r2, r3)
end

@inline function ymapping(x, y, z, dx, dy, dz)
    j = Utils.unsafe_floor(y / dy)
    k = Utils.unsafe_floor(z / dz + 0.5)
    i = Utils.unsafe_floor(x / dx + 0.5)

    r2 = y / dy - j
    r3 = z / dz - k + 0.5
    r1 = x / dx - i + 0.5

    Mapping(i, j, k, r1, r2, r3)
end

@inline function zmapping(x, y, z, dx, dy, dz)
    k = Utils.unsafe_floor(z / dz)
    i = Utils.unsafe_floor(x / dx + 0.5)
    j = Utils.unsafe_floor(y / dy + 0.5)

    r3 = z / dz - k
    r2 = y / dy - j + 0.5
    r1 = x / dx - i + 0.5

    Mapping(i, j, k, r1, r2, r3)
end
