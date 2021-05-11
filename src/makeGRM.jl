###### Julia implementation of VanRaden GRM

using LinearAlgebra
using PositiveFactorizations
using Statistics;

function GRM(M)

	M = Matrix(M);
	#M = 1.0 .+ Matrix(M); # change to 0,1,2
	P = mean(M, dims=1) ./2 ## freq
	varsum = 2 * sum(P .* (1.0 .- P))
	A = (M .- (2 .* P)) ;
	K = (A * A') .* (1/varsum);
	K = 0.99*K + 0.01*I
	return(K)

end


function GRMwted(M, D)

	M = Matrix(M);
	#M = 1.0 .+ Matrix(M); # change to 0,1,2
	P = mean(M, dims=1) ./2 ## freq
	varsum = 2 * sum(P .* (1.0 .- P))
	A = (M .- (2 .* P)) ;
	K = (A * D * A') .* (1/varsum);
	K = 0.99*K + 0.01*I
	return(K)

end

function GRMinv(M)

        M = Matrix(M);
        #M = 1.0 .+ Matrix(M);
        P = mean(M, dims=1) ./2
        varsum = 2 * sum(P .* (1.0 .- P))
        A = (M .- (2 .* P));
        K = (A * A') .* (1/varsum);
        K = 0.99*K + 0.01*I
	K = cholesky(Positive, K);
	K = inv(K);
        return(K)

end

function GRMwtedinv(M, D)

        M = Matrix(M);
        #M = 1.0 .+ Matrix(M);
        P = mean(M, dims=1) ./2
        varsum = 2 * sum(P .* (1.0 .- P))
        A = (M .- (2 .* P));
        K = (A * D * A') .* (1/varsum);
        K = 0.99*K + 0.01*I
	K = cholesky(Positive, K);
	K = inv(K);
        return(K)

end
