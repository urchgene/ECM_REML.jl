########################################################################
## EMMREML (univariate) translated to Julia for big data computing
## Translated by Uche Godfrey Okeke
## EMMREML was originally written in R by Deniz Akdemir and Uche G. Okeke
## This is a fast and big Data version with autodifferentiation in Julia...
## Large problem version using SuperLU
## use for Julia 1.6 and above
######################################################################### 

using Optim;
using SuperLU, SparseArrays;
using ForwardDiff, PositiveFactorizations;
using LinearAlgebra, DataFrames;

function emmreml_LU(y, X, Z, K, linenames)

q = size(X,2);
n = size(y,1);
spI = one(rand(n,n));

#### Hat matrix via SVD
U, D, V = svd(X);
S = spI - U*U';
ZK = Z*K;
offset = 0.000001;
ZKZt = ZK*Z';
ZKZtandoffset = ZKZt + (offset * I);
SZKZtSandoffset = (S * ZKZtandoffset)*S;

### Change to use SVD or eigen decomposition here...
U, D, V = svd(SZKZtSandoffset);
Ur = U[:, :1:(n - q)];
lambda = D[1:(n - q)] .- offset;
eta = Ur'y;


function minimfunc(delta)
  (n - q) * log.(sum(eta.^2 ./(lambda .+ delta))) + sum(log.(lambda .+ delta))
end

nvar = 1
lower = ([0.00000000001])
upper = ([Inf])
od = OnceDifferentiable(vars -> minimfunc(vars[1]), ones(nvar); autodiff=:forward);
inner_optimizer = LBFGS()
optimout = optimize(od, lower, upper, ones(nvar), Fminbox(inner_optimizer), Optim.Options(show_trace=true))

deltahat = Optim.minimizer(optimout);
deltahat = reshape(deltahat)[1];

### use superlu to get inv not pinv
Hinvhat = lu(sparse(ZKZt + (deltahat * spI)));
Hinvhat = inv(Hinvhat);
XtHinvhat = X'Hinvhat;

#### Do cholesky solve for betahat
F = lu(sparse(XtHinvhat * X));
betahat = F \ XtHinvhat * y;
ehat = y .- (X * betahat);
Hinvhatehat = Hinvhat * ehat;
sigmausqhat = sum(eta.^2 ./(lambda .+ deltahat))/(n - q);
Vinv = (1/sigmausqhat) * Hinvhat;
sigmaesqhat = deltahat * sigmausqhat;
uhat = ZK'Hinvhatehat;
df = n - q;
loglik = -0.5 * (Optim.minimum(optimout) + df + df * log.(2 * pi/df));


### also use LU solve here....
F = lu(sparse(X'Vinv * X));
jjj = F \ X'Vinv;
P = Vinv - Vinv * X * jjj;
varuhat = sigmausqhat.^2 * ZK'P * ZK;
PEVuhat = sigmausqhat * K - varuhat;

varbetahat = lu(sparse(X'Vinv * X));
varbetahat = inv(X'Vinv * X);

uhat = DataFrame(Lines=linenames, Uhat=uhat);
Vu = sigmausqhat; Ve = sigmaesqhat;  
varuhat = diag(varuhat); varbetahat = diag(varbetahat); PEVuhat = diag(PEVuhat);
h2 = Vu ./(Vu + Ve); rel = 1 .- (PEVuhat ./(Vu * diag(K)));

m11 =  
  Dict(
  :Vu => Vu,
  :Ve => Ve,
  :betahat => betahat,
  :uhat => uhat,
  :varuhat => varuhat,
  :varbetahat => varbetahat,
  :PEVuhat => PEVuhat,
  :loglik => loglik,
  :h2 => h2,
  :rel => rel
)

return(m11)
  
end
