export Contingency

"""
Description:

    Contingency factor node

    The contingency distribution is a multivariate generalization of
    the categorical distribution. As a bivariate distribution, the
    contingency distribution defines the joint probability
    over two unit vectors. The parameter p encodes a contingency matrix
    that specifies the probability of co-occurrence.

    out1 ∈ {0, 1}^d1 where Σ_j out1_j = 1
    out2 ∈ {0, 1}^d2 where Σ_k out2_k = 1
    p ∈ [0, 1]^{d1 × d2}, where Σ_jk p_jk = 1

    f(out1, out2, p) = Con(out1, out2 | p)
                     = Π_jk p_jk^{out1_j * out2_k}

    A Contingency distribution over more than two variables requires
    higher-order tensors as parameters; these are not implemented in ForneyLab.

Interfaces:

    1. out1
    2. out2
    3. p

Construction:

    Contingency(id=:some_id)
"""
mutable struct Contingency <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Contingency(out1, out2, p; id=generateId(Contingency))
        @ensureVariables(out1, out2, p)
        self = new(id, Array{Interface}(undef, 3), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out1] = self.interfaces[1] = associate!(Interface(self), out1)
        self.i[:out2] = self.interfaces[2] = associate!(Interface(self), out2)
        self.i[:p] = self.interfaces[3] = associate!(Interface(self), p)

        return self
    end
end

slug(::Type{Contingency}) = "Con"

format(dist::ProbabilityDistribution{Multivariate, Contingency}) = "$(slug(Contingency))(p=$(format(dist.params[:p])))"

ProbabilityDistribution(::Type{Multivariate}, ::Type{Contingency}; p=1/9*ones(3,3)) = ProbabilityDistribution{Multivariate, Contingency}(Dict(:p=>p))
ProbabilityDistribution(::Type{Contingency}; p=1/9*ones(3,3)) = ProbabilityDistribution{Multivariate, Contingency}(Dict(:p=>p))

dims(dist::ProbabilityDistribution{Multivariate, Contingency}) = length(size(dist.params[:p]))

vague(::Type{Contingency}, n_factors::Tuple{Int64, Int64}=(3,3)) = ProbabilityDistribution(Multivariate, Contingency, p=(1/prod(n_factors))*ones(n_factors))

isProper(dist::ProbabilityDistribution{Multivariate, Contingency}) = (abs(sum(dist.params[:p])-1.) < 1e-6)

# Entropy functional
function differentialEntropy(dist::ProbabilityDistribution{Multivariate, Contingency})
    -sum(dist.params[:p].*log.(clamp.(dist.params[:p], tiny, Inf)))
end