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

    function GainEqualityCompositeNode(A::Array=[1.0], use_composite_update_rules::Bool=true; name="unnamed")
        if use_composite_update_rules
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            # In case we don't use composite update rules, A is passed to the internal FixedGainNode.
            A = ensureMatrix(deepcopy(A))
        end
        self = new(A, use_composite_update_rules, name, Array(Interface, 3))

        # Define the internals of the composite node
        self.equality_node = EqualityNode(3, name="$(name)_internal_equality")
        self.fixed_gain_node = FixedGainNode(A, name="$(name)_internal_gain")
        Edge(self.equality_node.interfaces[2], self.fixed_gain_node.in1) # Internal edge

        # Initialize the composite node interfaces belonging to the composite node itself.
        self.interfaces[1] = Interface(self)
        self.interfaces[2] = Interface(self)
        self.interfaces[3] = Interface(self)
        # Initialize the interfaces as references to the internal node interfaces.
        self.interfaces[1].child = self.equality_node.interfaces[1]
        self.interfaces[2].child = self.equality_node.interfaces[3]
        self.interfaces[3].child = self.fixed_gain_node.out
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
        self.out = self.interfaces[3]
        # Set internal message passing schedules
        self.in1.internal_schedule = [self.fixed_gain_node.in1, self.equality_node.interfaces[1]]
        self.in2.internal_schedule = [self.fixed_gain_node.in1, self.equality_node.interfaces[3]]
        self.out.internal_schedule = [self.equality_node.interfaces[2], self.fixed_gain_node.out]

        return self
    end
end

############################################
# GaussianDistribution methods
############################################

# Rule set for backward propagation ({in2,out}-->in1 or {in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x - V_x * A' * inv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, V_x::Array{T, 2}, m_y::Array{T, 1}, V_y::Array{T, 2}) = m_x + V_x * A' * inv(V_y + A * V_x * A') * (m_y - A * m_x)

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GainEqualityCompositeNode,
                            inbound_messages_value_types::Type{GaussianDistribution})
    # Calculate an outbound message based on the inbound messages and the node function.
    # This function is not exported, and is only meant for internal use.

    if !node.use_composite_update_rules
        executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
    else
        if outbound_interface_id == 3
            # Forward message
            # We don't have a shortcut rule for this one, so we use the internal nodes to calculate the outbound msg
            executeSchedule(node.interfaces[outbound_interface_id].internal_schedule)
        elseif outbound_interface_id == 1 || outbound_interface_id == 2
            dist_out = getOrCreateMessage(node.interfaces[outbound_interface_id], GaussianDistribution).value

            # Backward messages
            dist_3 = node.interfaces[3].partner.message.value
            dist_in = node.interfaces[outbound_interface_id == 1 ? 2 : 1].partner.message.value # the other input interface

            # Select parameterization
            # Order is from least to most computationally intensive
            if dist_3.xi != nothing && dist_3.W != nothing && dist_in.xi != nothing && dist_in.W != nothing
                dist_out.m = nothing
                dist_out.V = nothing
                dist_out.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
                dist_out.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
            elseif dist_3.m != nothing && dist_3.V != nothing && dist_in.m != nothing && dist_in.V != nothing
                dist_out.m = backwardGainEqualityMRule(node.A, dist_in.m, dist_in.V, dist_3.m, dist_3.V)
                dist_out.V = backwardGainEqualityVRule(node.A, dist_in.V, dist_3.V)
                dist_out.W = nothing
                dist_out.xi = nothing
            elseif dist_3.m != nothing && dist_3.W != nothing && dist_in.m != nothing && dist_in.W != nothing
                # TODO: Not very efficient!
                dist_out.m = backwardGainEqualityMRule(node.A, dist_in.m, inv(dist_in.W), dist_3.m, inv(dist_3.W))
                dist_out.V = nothing
                dist_out.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
                dist_out.xi = nothing
            else
                # Fallback: convert inbound messages to (xi,W) parametrization and then use efficient rules
                ensureXiWParametrization!(dist_in)
                ensureXiWParametrization!(dist_3)
                dist_out.m = nothing
                dist_out.V = nothing
                dist_out.W = backwardGainEqualityWRule(node.A, dist_in.W, dist_3.W)
                dist_out.xi = backwardGainEqualityXiRule(node.A, dist_in.xi, dist_3.xi)
            end
        else
            error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
        end
    end
    # Return the outbound message
    return node.interfaces[outbound_interface_id].message
end

