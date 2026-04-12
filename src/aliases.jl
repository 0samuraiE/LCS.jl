module Aliases
export Field, Property, Tuple3, Integer3, Real3

const Field = AbstractArray{<:Real,3}
const Property{T} = AbstractVector{T}
const Tuple3{T} = NTuple{3,T}
const Integer3 = Tuple3{<:Integer}
const Real3 = Tuple3{<:Real}
end
