export SVI, Svi

"""
Description:

    Paper: https://www.jmlr.org/papers/volume14/hoffman13a/hoffman13a.pdf

    Stochastic Variational Inference node allows VMP to be scaled to large datasets
    by generating mirror of global parameters inside plates which is not often used
    in FFGs.

    Maps a global variable to a local mirror variable by

    f(out,in1) = δ(out - in1)

Interfaces:

    1. out
    2. in1

Construction:

    Svi(out, in1, opt=opt, q=q, batch_size=M, dataset_size=N, id=:some_id)

    where q: initial variational factor for global variable,
          opt: an optimizer from Flux.Optimise family together with ForgetDelayDescent defined in helpers.jl
          batch_size: batch size
          dataset_size: dataset size
"""

mutable struct SVI <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}
    q::ProbabilityDistribution
    q_memory::ProbabilityDistribution
    opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent}
    batch_size::Int
    dataset_size::Int

    function SVI(out, in1, opt, q, batch_size, dataset_size; id=generateId(SVI))
        @ensureVariables(out, in1)
        self = new(id, Array{Interface}(undef, 2), Dict{Int,Interface}(), q, deepcopy(q), opt, batch_size, dataset_size)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)

        return self
    end
end

slug(::Type{SVI}) = "svi"

function Svi(out::Variable, in1::Variable;
             opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent},
             q::ProbabilityDistribution,
             batch_size::Int, dataset_size::Int)
    SVI(out, in1, opt, q, batch_size, dataset_size)
end
