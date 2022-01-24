export Probit

using StatsFuns: normcdf

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

    return (Message{Gaussian{Moments}, Univariate}, ()) # Univariate only
end

# Average energy functional
function averageEnergy(::Type{Probit}, marg_out::Distribution{Univariate, Bernoulli}, marg_in1::Distribution{Univariate, F}) where F<:Gaussian

    # extract parameters
    (m, v) = unsafeMeanCov(marg_in1)
    p = unsafeMean(marg_out)
    h(x)  = -p*log(normcdf(x)) - (1-p)*log(normcdf(-x))

    gaussianQuadrature(h, m=m, v=v)
end

function averageEnergy(::Type{Probit}, marg_out::Distribution{Univariate, PointMass}, marg_in1::Distribution{Univariate, F}) where F<:Gaussian
    p = marg_out.params[:m]
    isnan(p) && (p = 0.5)
    (0.0 <= p <= 1.0) || error("Binary input $p must be between 0 and 1")

    return averageEnergy(Probit, Distribution(Univariate, Bernoulli, p=p), marg_in1)
end

function averageEnergy(::Type{Probit}, marg_out::Distribution{Univariate,PointMass}, marg_in1::Distribution{Univariate,SampleList})
    samples = marg_in1.params[:s]
    weights = marg_in1.params[:w]
    p = marg_out.params[:m]

    -p*sum(weights.*normlogcdf.(samples)) -
    (1 - p)*sum(weights.*normlogcdf.(-samples))
end