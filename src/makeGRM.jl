###### Julia implementation of VanRaden GRM

using LinearAlgebra
using PositiveFactorizations
using Statistics;

function GRM(M)

	M = Matrix(M);
	M = 1.0 .+ Matrix(M);
	freq = mean(M, dims=1) ./2
	P = 2 .* freq
	varsum = sum(P .* (1.0 .- freq))
	A = (M .- P) ./ varsum;
	A = A * A';
	A = 0.95*A + 0.05*I
	return(A)

end


function GRMinv(M)

        M = Matrix(M);
        M = 1.0 .+ Matrix(M);
        freq = mean(M, dims=1) ./2
        P = 2 .* freq
        varsum = sum(P .* (1.0 .- freq))
        A = (M .- P) ./ varsum;
        A = A * A';
        A = 0.95*A + 0.05*I
	A = cholesky(Positive, A);
	A = inv(A);
        return(A)

end

