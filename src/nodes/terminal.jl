export TerminalNode, PriorNode

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

typealias PriorNode TerminalNode # For more overview during graph construction

isDeterministic(node::TerminalNode) = (typeof(node.value) <: AbstractDelta) ? true : false # Edge case for deterministicness

# Implement firstFreeInterface since EqualityNode is symmetrical in its interfaces
firstFreeInterface(node::TerminalNode) = (node.interfaces[1].partner==nothing) ? node.interfaces[1] : error("No free interface on $(typeof(node)) $(node.id)")

collectAllOutboundTypes(rule::Function, call_signature::Vector, node::TerminalNode) = DataType[typeof(node.value)]


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
Compute average energy as U[q] = -E_q[log f(x)]
"""
function U(node::TerminalNode, dist_y::Gaussian)
    (typeof(node.value) <: Gaussian) || error("Average energy for TerminalNode only defined for Gaussian value")

    ensureParameters!(node.value, (:m, :W))
    ensureParameters!(dist_y, (:m, :V))

    mean = node.value.m
    prec = node.value.W

    return  (1/2)*log(2*pi) -
            0.5*log(prec) +
            0.5*(prec)*( dist_y.V + (dist_y.m - mean)^2)
end

function U{dims}(node::TerminalNode, dist_y::MvGaussian{dims})
    (typeof(node.value) <: MvGaussian) || error("Average energy for TerminalNode only defined for Gaussian value")

    ensureParameters!(node.value, (:m, :W))
    ensureParameters!(dist_y, (:m, :V))

    mean = node.value.m
    prec = node.value.W

    return  (dims/2)*log(2*pi) -
            0.5*log(det(prec)) +
            0.5*trace( prec*( dist_y.V + (dist_y.m - mean)*(dist_y.m - mean)' ) )
end