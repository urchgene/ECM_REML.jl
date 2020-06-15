###################################################################
##### This code runs RKHS and gets a cholesky of it for fitting 
##### MT-RKHS model...
##### Developer: Uche Godfrey Okeke
###################################################################


###### Julia implementation of RKHS

using LinearAlgebra
using PositiveFactorizations;
 

##### SqEuclidean calculation #####
function SqEuclid(M)
           X = M; Y= M; K = size(M,1); N = size(M,1);
           D = sum(X.^2, dims=2) * ones(1,K) + ones(N,1)*sum( Y.^2, dims=2 )' - 2 .* X*Y'
           #D = sum(X.^2, 2) * ones(1,K) + ones(N,1)*sum(Y.^2, 2)' - 2 .* X*Y'
end
##################################

function RKHS(M, w)

	M = Matrix(M);
	M = 1.0 .+ Matrix(M);
	w = Float64.(w)
	w = 1/(2*(w*w))
	K = SqEuclid(M);
	A = exp.(-w*K)
	A = 0.95*A + (1-0.95)*I  ## bending
	A = A + 0.001*I
	return(A)

end


function RKHSinv(M, w)

        M = Matrix(M);
        M = 1.0 .+ Matrix(M);
        w = Float64.(w)
        w = 1/(2*(w*w))
        K = SqEuclid(M);
        A = exp.(-w*K)
        A = 0.95*A + (1-0.95)*I  ## bending
        A = A + 0.001*I
	A = cholesky(Positive, A);
	A = inv(A);
        return(A)

end

