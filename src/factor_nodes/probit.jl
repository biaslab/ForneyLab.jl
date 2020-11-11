export Probit

"""
Description:

    Constrains a continuous, real-valued variable in1 ∈ R with a binary (boolean) variable out ∈ {0, 1} through a probit link function.

    f(out, in1) = Ber(out | Φ(in1))

Interfaces:

    1. out (binary)
    2. in1 (real)

Construction:

    Probit(out, in1, id=:some_id)
"""
mutable struct Probit <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Probit(out, in1; id=generateId(Probit))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{Probit}) = "Φ"

function breakerParameters(interface::Interface, partner_interface::Interface, partner_node::Probit)
    (partner_interface == partner_node.i[:in1]) || error("Breaker initialization requested for non-breaker interface: $(interface)")

    return (Message{GaussianMeanVariance, Univariate}, ()) # Univariate only
end

# Average energy functional
function averageEnergy(::Type{Probit}, marg_out::ProbabilityDistribution{Univariate, Bernoulli}, marg_in1::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    (marg_in1_m, marg_in1_v) = unsafeMeanCov(marg_in1)
    h = (x -> log(0.5*erf(x) + 0.5 + tiny)) # Add `tiny` for numeric stability

    (1 - marg_out.params[:p])*gaussianQuadrature(h, m=-marg_in1_m, v=marg_in1_v) +
    marg_out.params[:p]*gaussianQuadrature(h, m=marg_in1_m, v=marg_in1_v)
end

function averageEnergy(::Type{Probit}, marg_out::ProbabilityDistribution{Univariate, PointMass}, marg_in1::ProbabilityDistribution{Univariate, F}) where F<:Gaussian
    p = marg_out.params[:m]
    isnan(p) && (p = 0.5)
    (0.0 <= p <= 1.0) || error("Binary input $p must be between 0 and 1")

    return averageEnergy(Probit, ProbabilityDistribution(Univariate, Bernoulli, p=p), marg_in1)
end

function averageEnergy(::Type{Probit}, marg_out::ProbabilityDistribution{Univariate,PointMass}, marg_in1::ProbabilityDistribution{Univariate,SampleList})
    samples = marg_in1.params[:s]
    weights = marg_in1.params[:w]
    p = marg_out.params[:m]

    -p*sum(weights.*normlogcdf.(samples)) -
    (1 - p)*sum(weights.*normlogcdf.(-samples))
end