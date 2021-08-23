export CVI, Cvi

"""
Description:

    Paper: https://arxiv.org/pdf/1703.04265.pdf

    Conjugate Computation Variational Inference node allows VMP to be applied to nonconjugate factor pairs.

    Maps a variable through

    f(out,in1,in2,in3,...) = δ(out - g(in1,in2,in3,...))

Interfaces:

    1. out
    2. in1,in2,in3,...

Construction:

    Cvi(out, in1,in2,in3,..., g=g, opt=opt,
        num_iterations=[I1,I2,I3,...], num_samples=I,
        q=[q1,q2,q3,...], proper_message=B, online_inference=[B1,B2,B3,...]
        batch_size=M, dataset_size=N, id=:some_id)

    where q: initial variational factors vector,
          opt: an optimizer from Flux.Optimise family together with ForgetDelayDescent defined in helpers.jl
          num_iterations: Vector of integers, num iterations in optimization per variable
          num_samples: # of samples to be employed in estimations
          proper_message: flag to indicate if the message is needed to be a proper dist.
          online_inference: vector of bools show if online inference rules are applied per variable
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

function Cvi(out::Variable, args::Vararg; g::Function, opt::Any,
             num_samples::Union{Int,Vector{Int}}, num_iterations::Union{Int,Vector{Int}},
             q=[ProbabilityDistribution(Univariate,GaussianMeanVariance,m=0,v=1)], infer_memory=0,
             proper_message=false, online_inference=false,
             batch_size=1, dataset_size=1, id=ForneyLab.generateId(CVI))
    CVI(id, g, opt, num_iterations, num_samples, q, infer_memory, proper_message,
        online_inference, batch_size, dataset_size, out, args...)
end
