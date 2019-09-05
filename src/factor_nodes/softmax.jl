export Softmax

"""
Description:

    Softmax mapping between a real variable in1 ∈ R^d and discrete variable out ∈ {0, 1}^d.

    f(out, in1, xi, a) = Cat(out_j | exp(in1_j)/Σ_k exp(in1_k)),
    where log(Σ_k exp(x_k)) is upper-bounded according to (Bouchard, 2007).

Interfaces:

    1. out (discrete)
    2. in1 (real)
    3. xi  (auxiliary variable)
    4. a   (auxiliary variable)

Construction:

    Softmax(out, in1, xi, a)
"""
mutable struct Softmax <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}

    function Softmax(out, in1, xi, a; id=generateId(Softmax))
        @ensureVariables(out, in1, xi, a)
        self = new(id, Vector{Interface}(undef, 4), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:xi] = self.interfaces[3] = associate!(Interface(self), xi)
        self.i[:a] = self.interfaces[4] = associate!(Interface(self), a)

        return self
    end
end

slug(::Type{Softmax}) = "σ"

function unsafeBoundMean(marg_in1::ProbabilityDistribution{Multivariate}, marg_xi::ProbabilityDistribution{Multivariate}, marg_a::ProbabilityDistribution{Univariate})
    xi_hat = unsafeMode(marg_xi)
    a_hat = unsafeMode(marg_a)
    x_bar = unsafeMean(marg_in1)
    v_x = unsafeVar(marg_in1)
    d = length(x_bar)

    s = 0.5*(x_bar .- a_hat - xi_hat) + 
        logisticLambda.(xi_hat).*(v_x + (x_bar .- a_hat).^2 - xi_hat) +
        log.(1 .+ exp.(xi_hat))
    
    return a_hat + sum(s)
end

# Average energy functionals
function averageEnergy(::Type{Softmax}, marg_out::ProbabilityDistribution, marg_in1::ProbabilityDistribution{Multivariate}, marg_xi::ProbabilityDistribution{Multivariate}, marg_a::ProbabilityDistribution{Univariate})
    return unsafeBoundMean(marg_in1, marg_xi, marg_a) - unsafeMean(marg_out)'*unsafeMean(marg_in1)
end