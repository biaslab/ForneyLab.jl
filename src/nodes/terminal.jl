export TerminalNode, PriorNode, DataNode

"""
Description:

    Sends out a predefined message.

        out
    [T]----->

    out = T.value

Interfaces:

    1 i[:out]

Construction:

    TerminalNode(Gaussian(), id=:my_node)
"""
type TerminalNode <: Node
    id::Symbol
    value::ProbabilityDistribution
    interfaces::Vector{Interface}
    i::Dict{Symbol,Interface}

    function TerminalNode(value::ProbabilityDistribution; id=generateNodeId(TerminalNode))
        self = new(id, deepcopy(value), Vector{Interface}(1), Dict{Symbol,Interface}())
        addNode!(currentGraph(), self)

        self.i[:out] = self.interfaces[1] = Interface(self)

        return self
    end
end

TerminalNode(num::Number; id=generateNodeId(TerminalNode)) = TerminalNode(convert(Delta, num), id=id)

TerminalNode(bool::Bool; id=generateNodeId(TerminalNode)) = TerminalNode(convert(Delta, bool), id=id)

TerminalNode{T<:Number}(vect::Vector{T}; id=generateNodeId(TerminalNode)) = TerminalNode(convert(MvDelta, vect), id=id)

TerminalNode{T<:Number}(mat::AbstractMatrix{T}; id=generateNodeId(TerminalNode)) = TerminalNode(convert(MatrixDelta, mat), id=id)

TerminalNode(; id=generateNodeId(TerminalNode)) = TerminalNode(vague(Gaussian), id=id)

# Aliases for more overview during graph construction
typealias PriorNode TerminalNode
typealias DataNode TerminalNode

isDeterministic(node::TerminalNode) = (typeof(node.value) <: AbstractDelta) ? true : false # Edge case for deterministicness

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::TerminalNode) = (node.interfaces[1].partner==nothing) ? node.interfaces[1] : error("No free interface on $(typeof(node)) $(node.id)")

collectAllOutboundTypes(rule::Function, call_signature::Vector, node::ForneyLab.TerminalNode) = DataType[typeof(node.value)]


############################################
# Sumproduct rules
############################################

"""
TerminalNode

    [T]--->
        -->
"""
function sumProductRule!(   node::TerminalNode,
                            outbound_interface_index::Type{Val{1}},
                            outbound_dist::Any,
                            msg_out::Any)

    # Fill the fields of outbound_dist with node.value
    return injectParameters!(outbound_dist, node.value)
end


############################
# Average energy functionals
############################

"""
Compute average energy
"""
averageEnergy(::Type{TerminalNode}, value::Gaussian, marg_out::Univariate) = averageEnergy(GaussianNode, Delta(unsafeMean(value)), Delta(1.0/unsafeCov(value)), marg_out)
averageEnergy(::Type{TerminalNode}, value::MvGaussian, marg_out::Multivariate) = averageEnergy(GaussianNode, MvDelta(unsafeMean(value)), MatrixDelta(cholinv(unsafeCov(value))), marg_out)

function averageEnergy(::Type{TerminalNode}, value::Beta, marg_out::Univariate)
    lbeta(value.a, value.b) -
    (value.a - 1)*unsafeLogMean(marg_out) -
    (value.b - 1)*unsafeMirroredLogMean(marg_out)
end