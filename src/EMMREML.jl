module EMMREML

using Optim;
using Statistics;
using ForwardDiff, PositiveFactorizations;
using LinearAlgebra, DataFrames;

include("emmremlJulia.jl")
include("emmremlMultivariate_Varcomp.jl")
include("makeGRM.jl")
include("makeRKHS.jl")

export emmreml, emmremlMultivariate
export GRM, GRMinv, RKHS, RKHSinv, SqEuclid
export GRMwted, GRMwtedinv, GRMiter

end # module
