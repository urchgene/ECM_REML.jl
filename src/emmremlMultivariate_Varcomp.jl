
#################################################
### EMMREML-MULTIVARIATE for Julia
### Modified and translated by Uche Godfrey Okeke
### Suitable for Big data purposes
### This is specifically for estimating Vg and Ve
##################################################


using LinearAlgebra;
#using RCall;
using Statistics;

#BLAS.set_num_threads(20)
#my_BLAS_set_num_threads(n) =  ccall((:openblas_set_num_threads64_, Base.libblas_name), Cvoid, (Int32,), n)
#my_BLAS_set_num_threads(20)

function emmremlMultivariate(Y, X, Z, K)

## call with regular Y, X, K but transpose Z ###

	Z = Z; X = X'; Y = Y'; tolpar = 1e-06; tolparinv = 1e-06;


    function ECM1(ytl, xtl, Vgt, Vet, Bt, deltal)
        Vlt = deltal * Vgt + Vet
        invVlt = inv(Vlt + tolparinv * I)
        gtl = deltal * Vgt * invVlt * (ytl - Bt * xtl)
        Sigmalt = deltal * Vgt - deltal * Vgt * invVlt * (deltal * Vgt)
        return(Dict(:Vlt => Vlt, :gtl => gtl, :Sigmalt => Sigmalt))
    end

    function wrapperECM1(l)
        ytl = Yt[:, l]
        xtl = Xt[:, l]
        deltal = eigZKZt.values[l]
        return(ECM1(ytl, xtl, Vgt, Vet, Bt, deltal))
    end


    function Vgfunc(l, outfromECM1)
        Vgl = outfromECM1[l][:gtl] * outfromECM1[l][:gtl]'
        return((1/n) * (1 ./ eigZKZt.values[l]) * (Vgl + outfromECM1[l][:Sigmalt]))
    end


    function Vefunc(l, outfromECM1)
        etl = Yt[:, l] - Bt * Xt[:, l] - outfromECM1[l][:gtl]
        return((1/n) * (etl*etl' + outfromECM1[l][:Sigmalt]))
    end

    if (!ismissing(Y))
        N = size(K,1)
        KZt = K * Z'
        ZKZt = Z * KZt
        ZKZt  = ZKZt + 0.0001*I;
        eigZKZt = eigen(ZKZt)
        n = size(ZKZt, 1)
        d = size(Y, 1)
        Yt = Y * eigZKZt.vectors
        Xt = X * eigZKZt.vectors
        Vgt = cov(Y')/2
        Vet = cov(Y')/2
        XttinvXtXtt = Xt' * pinv(Xt*Xt')
        Bt = Yt * XttinvXtXtt
        Vetm1 = Vet

        j = 1;
        while true
            println("iteration ...", j)
            outfromECM1 = [wrapperECM1(l) for l in 1:n]
            Vetm1 = Vet
            hh = [outfromECM1[i][:gtl] for i in 1:n]; 
            Gt = reduce(hcat, hh);
            Bt = (Yt - Gt) * XttinvXtXtt
            listVgts = [Vgfunc(l, outfromECM1) for l in 1:n]
            Vgt = reduce(+, listVgts)
            listVets = [Vefunc(l, outfromECM1) for l in 1:n]
            Vet = reduce(+, listVets)
            convnum = abs(sum(diag(Vet - Vetm1)))/abs(sum(diag(Vetm1)))

            if convnum < tolpar
                break
            end
            j += 1

        end

    end

    h2 = diag(Vgt ./ (Vgt + Vet))


    ### Write ouput in R for pretty compatible format
    
    m11 =  Dict(
			:Vg => Vgt,
			:Ve => Vet,
			:h2 => h2)

    return(m11)

end 



