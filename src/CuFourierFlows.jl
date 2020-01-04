using .CuArrays

# Discard `effort` argument for CuArrays
plan_flows_fft(a::CuArray, effort) = plan_fft(a)
plan_flows_rfft(a::CuArray, effort) = plan_rfft(a)

OneDGrid(dev::GPU, args...; kwargs...) = OneDGrid(args...; ArrayType=CuArray, kwargs...)
TwoDGrid(dev::GPU, args...; kwargs...) = TwoDGrid(args...; ArrayType=CuArray, kwargs...)

function Base.zeros(::GPU, T, dims)
  a = CuArray{T}(undef, dims...)
  a .= 0
  return a
end

ArrayType(::GPU, T, dim) = CuArray{T, dim}
ArrayType(::GPU) = CuArray

supersize(a::CuArray) = size(a)

getetdcoeffs(dt, L::CuArray; kwargs...) = 
    (CuArray(ζ) for ζ in getetdcoeffs(dt, Array(L); kwargs...))

makefilter(K::CuArray; kwargs...) = CuArray(makefilter(Array(K); kwargs...))

function makefilter(g::AbstractGrid{Tg, <:CuArray}, T, sz; kwargs...) where Tg
    CuArray(ones(T, sz)) .* makefilter(g; realvars=sz[1]==g.nkr, kwargs...)
end

import Base.exp, Base.sin, Base.cos
using CUDAnative

exp(x::Complex{Float32}) = CUDAnative.exp(x)
exp(x::Complex{Float64}) = CUDAnative.exp(x)
exp(x::Float32) = CUDAnative.exp(x)
exp(x::Float64) = CUDAnative.exp(x)

sin(x::Float32) = CUDAnative.sin(x)
sin(x::Float64) = CUDAnative.sin(x)

cos(x::Float32) = CUDAnative.cos(x)
cos(x::Float64) = CUDAnative.cos(x)
