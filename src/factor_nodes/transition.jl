export Transition

"""
Description:

    The transition node models a transition between discrete
    random variables, with node function

    f(out, in1, a) = Cat(out | a*in1)

Interfaces:

    1. out
    2. in1
    3. a

Construction:

    Transition(out, in1, a, id=:some_id)
"""
mutable struct Transition <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function Transition(out, in1, a; id=generateId(Transition))
        @ensureVariables(out, in1, a)
        self = new(id, Array{Interface}(undef, 3), Dict{Int,Interface}())
        addNode!(currentGraph(), self)
        self.i[:out] = self.interfaces[1] = associate!(Interface(self), out)
        self.i[:in1] = self.interfaces[2] = associate!(Interface(self), in1)
        self.i[:a] = self.interfaces[3] = associate!(Interface(self), a)

        return self
    end
end

slug(::Type{Transition}) = "T"

# Average energy functional
function averageEnergy(::Type{Transition}, marg_out::ProbabilityDistribution, marg_in1::ProbabilityDistribution, marg_a::ProbabilityDistribution{MatrixVariate})
    -unsafeMeanVector(marg_out)'*unsafeLogMean(marg_a)*unsafeMeanVector(marg_in1)
end

function averageEnergy(::Type{Transition}, marg_out_in1::ProbabilityDistribution{Multivariate, Contingency}, marg_a::ProbabilityDistribution{MatrixVariate})
    -tr(marg_out_in1.params[:p]'*unsafeLogMean(marg_a))
end
