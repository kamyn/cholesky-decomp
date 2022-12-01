include("cholesky_decomposition.jl")
using LinearAlgebra;
using BenchmarkTools;
using Plots;
using Main.CholeskyDecomp;

ns = [10, 20, 30, 40, 50, 100, 150, 200, 300, 400, 500, 1000, 1500, 2000]
medians = zeros(lastindex(ns), 1)
for idx in 1:lastindex(ns)
    local n = ns[idx]
    local A = 100*rand(n, n)
    A = A'*A # symmetric positive definite matrix
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = 1
    t = @benchmark CholeskyDecomp.decomposition($A)
    medians[idx] = time(median(t)) / 1000000000
end

plot(ns, medians, seriestype=:scatter, label=false)
xlabel!("n")
ylabel!("time, seconds")