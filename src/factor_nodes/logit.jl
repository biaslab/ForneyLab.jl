export Logit

"""
Description:

    Logit mapping between a real variable in1 ∈ R and a binary variable out ∈ {0, 1}.

    f(out, in1, xi) =  Ber(out | σ(in1))
                    >= exp(in1*out) σ(xi) exp[-(in1 + xi)/2 - λ(xi)(in1^2 - xi^2)], where
               σ(x) =  1/(1 + exp(-x))
               λ(x) =  (σ(x) - 1/2)/(2*x)

Interfaces:

    1. out (binary)
    2. in1 (real)
    3. xi (auxiliary variable)

Construction:

    Logit(out, in1, xi)
"""
mutable struct Logit <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}

    function Logit(out, in1, xi; id=generateId(Logit))
        @ensureVariables(out, in1, xi)
        self = new(id, Vector{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:xi] = self.interfaces[3] = associate!(Interface(self), xi)

        return self
    end
end

slug(::Type{Logit}) = "σ"

logisticSigmoid(x::Float64) = 1/(1 + exp(-x))
logLogisticSigmoid(x::Float64) = -log(1 + exp(-x))
logisticLambda(x::Float64) = (logisticSigmoid(x) - 0.5)/(2*x)

# Average energy functionals
function averageEnergy(::Type{Logit}, marg_out::ProbabilityDistribution{Univariate}, marg_in1::ProbabilityDistribution{Univariate}, marg_xi::ProbabilityDistribution{Univariate})
    xi_hat = unsafeMode(marg_xi)

    logisticLambda(xi_hat)*(unsafeMean(marg_in1)^2 + unsafeCov(marg_in1) - xi_hat^2) +
    0.5*(unsafeMean(marg_in1) + xi_hat) -
    logLogisticSigmoid(xi_hat) -
    unsafeMean(marg_in1)*unsafeMean(marg_out)
end