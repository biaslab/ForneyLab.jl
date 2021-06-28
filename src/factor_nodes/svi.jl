export SVI, Svi

"""
Description:

    Stochastic Variational Inference node allows VMP to be scaled to large datasets
    by generating mirror of global parameters inside plates which is not often used
    in FFGs.

    Maps a location to a scale parameter by exponentiation

    f(out,in1) = δ(out - in1)

Interfaces:

    1. out
    2. in1

Construction:

    Svi(out, in1, q, opt, M, N, id=:some_id)

    where q: initial variational factor for global variable,
          opt: an optimizer from Flux.Optimise family together with ForgetDelayDescent defined in helpers.jl
          M: batch size
          N: dataset size
"""

mutable struct SVI <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}
    q::ProbabilityDistribution
    opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent}
    batch_size::Int
    dataset_size::Int

    function SVI(out, in1, q, opt, batch_size, dataset_size; id=generateId(SVI))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Int,Interface}(), q, opt, batch_size, dataset_size)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{SVI}) = "svi"

function Svi(out::Variable, in1::Variable, q::ProbabilityDistribution, opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent}, batch_size::Int, dataset_size::Int)
    SVI(out, in1, q, opt, batch_size, dataset_size)
end
