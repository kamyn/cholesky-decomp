module CholeskyDecomp
    export decomposition
    using LinearAlgebra;

    # Cholesky decomposition for positive definite symmetric matrix A
    function decomposition(A::Matrix{Float64})
        local n = size(A)[1]
        L = zeros(n, n)
        L[1,1] = sqrt(A[1,1])
        for j in 2:n
            L[j,1] = A[j,1]/L[1,1]
        end
        for i in 2:n
            L[i,i] = sqrt(A[i,i] - sum(x -> x^2, L[i,1:i-1]))
            for j in i+1:n
                L[j,i] = (A[j,i] - sum([L[i,p] * L[j,p] for p in 1:i-1])) / L[i,i]
            end
        end
        return L
    end
end