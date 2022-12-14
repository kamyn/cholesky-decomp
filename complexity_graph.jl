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
    BenchmarkTools.DEFAULT_PARAMETERS.seconds = n/10
    local t = @benchmark CholeskyDecomp.decomposition($A)
    medians[idx] = time(median(t)) / 1000000000
end

nrate = 0.
for i in 1:lastindex(ns)-1
    for j in i+1:lastindex(ns)
        global nrate += medians[j] / medians[i] * 2^3/(ns[j]/ns[i])^3
    end
end
nrate = 2*nrate / ( (lastindex(ns)-1)*lastindex(ns) ) 

plot(ns, medians, seriestype=:scatter, label=false)
xlabel!("n")
ylabel!("time, seconds")

println("Increase time from 2n input data: ", nrate)