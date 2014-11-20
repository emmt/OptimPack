include("OptimPack.jl")

function rosenbrock_init!{T<:Real}(x0::Array{T,1})
  x0[1:2:end] = -1.2
  x0[2:2:end] =  1.0
end

function rosenbrock_fg!{T<:Real}(x::Array{T,1}, gx::Array{T,1})
  const c1::T = 1
  const c2::T = 2
  const c10::T = 10
  const c200::T = 200
  x1 = x[1:2:end]
  x2 = x[2:2:end]
  t1 = c1 - x1
  t2 = c10*(x2 - x1.*x1)
  g2 = c200*(x2 - x1.*x1)
  gx[1:2:end] = -c2*(x1.*g2 + t1)
  gx[2:2:end] = g2
  return sum(t1.*t1) + sum(t2.*t2)
end

function rosenbrock_test(n::Integer=20, m::Integer=3; single::Bool=false)
  T = (single ? Float32 : Float64)
  x0 = Array(T, n)
  rosenbrock_init!(x0)
  lbfgs(rosenbrock_fg!, x0, m, verb=true)
end

#space = OptimPack.OptimPackShapedVectorSpace(Float64, 10)

x0 = Array(Float64, 20)
rosenbrock_init!(x0)
x1 = OptimPack.nlcg(rosenbrock_fg!, x0, verb=true)
println(x1)
x2 = OptimPack.vmlm(rosenbrock_fg!, x0, verb=true)
println(x2)

x0 = Array(Float32, 20)
rosenbrock_init!(x0)
x1 = OptimPack.nlcg(rosenbrock_fg!, x0, verb=true)
println(x1)
x2 = OptimPack.vmlm(rosenbrock_fg!, x0, verb=true)
println(x2)

