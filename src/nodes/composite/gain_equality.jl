############################################
# GainEqualityCompositeNode
############################################
# Description:
#   Composite gain-equality node.
#
#        _________
#    in1 |       | in2
#   -----|->[=]<-|-----
#        |   |   |
#        |   v   |
#        |  [A]  |
#        |___|___|
#            | out
#            v
#
#   out = A*in1 = A*in2
#   Example:
#       GainEqualityCompositeNode([1.0]; name="my_node")
#   Gain A is fixed and defined through the constructor.
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       Message{GaussianDistribution}
#   2. in2:
#       Message{GaussianDistribution}
#   3. out:
#       Message{GaussianDistribution}
#
############################################

export GainEqualityCompositeNode

type GainEqualityCompositeNode <: CompositeNode
    # Basic node properties.
    # use_composite_update_rules is a flag for indicating use of
    # shortcut update rules. When set to true (default), the outbound
    # message is calculated via the provided update rules.
    # When set to false, messages are passed through the internal
    # graph of the composite node.
    A::Array
    use_composite_update_rules::Bool
    name::ASCIIString
    interfaces::Array{Interface,1}
    # Pointers to internal nodes
    equality_node::EqualityNode
    fixed_gain_node::FixedGainNode
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface

    function GainEqualityCompositeNode(A::Union(Array{Float64},Float64)=1.0, use_composite_update_rules::Bool=true; name=unnamedStr(), args...)
        if typeof(A)==Float64
            A = fill!(Array(Float64,1,1),A)
        elseif use_composite_update_rules
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            # In case we don't use composite update rules, A is passed to the internal FixedGainNode.
            A = ensureMatrix(deepcopy(A))
        end
        self = new(A, use_composite_update_rules, name, Array(Interface, 3))

        # Define the internals of the composite node
        self.equality_node = EqualityNode(name="$(name)_internal_equality")
        self.fixed_gain_node = FixedGainNode(A, name="$(name)_internal_gain")
        Edge(self.equality_node.interfaces[2], self.fixed_gain_node.in1, GaussianDistribution, add_to_graph=false) # Internal edge

        named_handle_list = [:in1, :in2, :out]
        for i = 1:length(named_handle_list)
            self.interfaces[i] = Interface(self) # Initialize the composite node interfaces belonging to the composite node itself.
            setfield!(self, named_handle_list[i], self.interfaces[i]) # Init named interface handles
        end

        # Initialize the interfaces as references to the internal node interfaces.
        self.interfaces[1].child = self.equality_node.interfaces[1]
        self.interfaces[2].child = self.equality_node.interfaces[3]
        self.interfaces[3].child = self.fixed_gain_node.out

        # Set internal message passing schedules
        self.in1.internal_schedule = convert_to_schedule([self.fixed_gain_node.in1, self.equality_node.interfaces[1]])
        self.in2.internal_schedule = convert_to_schedule([self.fixed_gain_node.in1, self.equality_node.interfaces[3]])
        self.out.internal_schedule = convert_to_schedule([self.equality_node.interfaces[2], self.fixed_gain_node.out])

        return self
    end
end

isDeterministic(::GainEqualityCompositeNode) = true

############################################
# GaussianDistribution methods
############################################

# Rule set for backward propagation ({in2,out}-->in1 or {in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x - V_x * A' * inv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, V_x::Array{T, 2}, m_y::Array{T, 1}, V_y::Array{T, 2}) = m_x + V_x * A' * inv(V_y + A * V_x * A') * (m_y - A * m_x)

function sumProduct!(node::GainEqualityCompositeNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Nothing)
    # Forward message (towards out)
    # We don't have a shortcut rule for this one, so we always use the internal nodes to calculate the outbound msg
    return node.interfaces[outbound_interface_id].message = executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
end

function sumProduct!(node::GainEqualityCompositeNode,
                            outbound_interface_id::Int,
                            msg_in1::Message{GaussianDistribution},
                            msg_in2::Nothing,
                            msg_out::Message{GaussianDistribution})
    # Backward message (towards in2)
    return applyBackwardRule!(node, outbound_interface_id, msg_in1, msg_out)
end

function sumProduct!(node::GainEqualityCompositeNode,
                            outbound_interface_id::Int,
                            msg_in1::Nothing,
                            msg_in2::Message{GaussianDistribution},
                            msg_out::Message{GaussianDistribution})
    # Backward message (towards in1)
    return applyBackwardRule!(node, outbound_interface_id, msg_in2, msg_out)
end

function applyBackwardRule!(node::GainEqualityCompositeNode,
                            outbound_interface_id::Int,
                            msg_in::Message{GaussianDistribution},
                            msg_out::Message{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.
    # Backward message (towards in1 or in2)
    if !node.use_composite_update_rules
        node.interfaces[outbound_interface_id].message = executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
    else
        dist_result = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).payload
        dist_3 = msg_out.payload
        dist_in = msg_in.payload

        # Select parameterization
        # Order is from least to most computationally intensive
        if dist_3.xi != nothing && dist_3.W != nothing && dist_in.xi != nothing && dist_in.W != nothing
            dist_result.m = nothing
            dist_result.V = nothing
            dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
            dist_result.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
        elseif dist_3.m != nothing && dist_3.V != nothing && dist_in.m != nothing && dist_in.V != nothing
            dist_result.m = backwardGainEqualityMRule(node.A, dist_in.m, dist_in.V, dist_3.m, dist_3.V)
            dist_result.V = backwardGainEqualityVRule(node.A, dist_in.V, dist_3.V)
            dist_result.W = nothing
            dist_result.xi = nothing
        elseif dist_3.m != nothing && dist_3.W != nothing && dist_in.m != nothing && dist_in.W != nothing
            dist_result.m = backwardGainEqualityMRule(node.A, dist_in.m, inv(dist_in.W), dist_3.m, inv(dist_3.W))
            dist_result.V = nothing
            dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
            dist_result.xi = nothing
        else
            # Fallback: convert inbound messages to (xi,W) parametrization and then use efficient rules
            ensureXiWParametrization!(dist_in)
            ensureXiWParametrization!(dist_3)
            dist_result.m = nothing
            dist_result.V = nothing
            dist_result.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
            dist_result.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
        end
    end
    # Return the outbound message
    return node.interfaces[outbound_interface_id].message
end