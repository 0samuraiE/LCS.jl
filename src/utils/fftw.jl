using AbstractFFTs: Plan
using LinearAlgebra: mul!, ldiv!

rfft!(A_hat, p::Plan, A) = mul!(A_hat, p, A)
irfft!(A, p::Plan, A_hat) = ldiv!(A, p, A_hat)
