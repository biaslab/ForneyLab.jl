export CVI, Cvi

"""
Description:

    Conjugate Computation Variational Inference node allows VMP to be applied to nonconjugate factor pairs.

    Maps a variable through

    f(out,in1) = δ(out - g(in1))

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

mutable struct CVI <: DeltaFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    g::Function # Vector function that expresses the output as a function of the inputs
    #opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent, Vector{Any}}
    opt::Any
    num_iterations::Union{Int,Vector{Int}}
    #num_samples::Union{Int,Vector{Int}}
    num_samples::Int
    #q::Union{ProbabilityDistribution,Vector{ProbabilityDistribution}}
    q::Vector{<:ProbabilityDistribution}
    q_memory::Vector{<:ProbabilityDistribution}
    infer_memory::Int
    proper_message::Bool
    online_inference::Union{Bool,Vector{Bool}}
    batch_size::Int
    dataset_size::Int

    function CVI(id::Symbol, g::Function,
                    #opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent, Vector{Any}},
                    opt::Any,
                    num_iterations::Union{Int,Vector{Int}}, num_samples::Union{Int,Vector{Int}}, q::Vector{<:ProbabilityDistribution},
                    infer_memory::Int, proper_message::Bool, online_inference::Union{Bool,Vector{Bool}},
                    batch_size::Int, dataset_size::Int, out::Variable, args::Vararg)
        @ensureVariables(out)
        n_args = length(args)
        for i=1:n_args
            @ensureVariables(args[i])
        end
        self = new(id, Array{Interface}(undef, n_args+1), Dict{Int,Interface}(), g, opt, num_iterations, num_samples,
                    q, deepcopy(q), infer_memory, proper_message, online_inference, batch_size, dataset_size)
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        for k = 1:n_args
            self.i[:in*k] = self.interfaces[k+1] = associate!(Interface(self), args[k])
        end

        return self
    end

end

slug(::Type{CVI}) = "cvi"

function Cvi(out::Variable, args::Vararg; g::Function,
    #opt::Union{Descent, Momentum, Nesterov, RMSProp, ADAM, ForgetDelayDescent, Vector{Any}},
             opt::Any,
                num_samples::Union{Int,Vector{Int}}, num_iterations::Union{Int,Vector{Int}},
                q=[ProbabilityDistribution(Univariate,GaussianMeanVariance,m=0,v=1)], infer_memory=0,
                proper_message=false, online_inference=false,
                batch_size=1, dataset_size=1, id=ForneyLab.generateId(CVI))
    CVI(id, g, opt, num_iterations, num_samples, q, infer_memory, proper_message,
        online_inference, batch_size, dataset_size, out, args...)
end
