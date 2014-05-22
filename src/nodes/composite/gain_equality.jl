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
#       GaussianMessage
#   2. in2:
#       GaussianMessage
#   3. out:
#       GaussianMessage
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
    parent::Union(CompositeNode, Nothing)
    interfaces::Array{Interface,1}
    # Pointers to internal nodes
    equality_node::EqualityNode
    fixed_gain_node::FixedGainNode
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface

    function GainEqualityCompositeNode(A::Array, use_composite_update_rules::Bool, parent::Union(CompositeNode, Nothing)=nothing; args...)
        (name = getArgumentValue(args, :name))!=false || (name = "unnamed")

        if use_composite_update_rules
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            # In case we don't use composite update rules, A is passed to the internal FixedGainNode.
            A = ensureMatrix(deepcopy(A))
        end
        self = new(A, use_composite_update_rules, name, parent, Array(Interface, 3))

        # Define the internals of the composite node
        self.equality_node = EqualityNode(3, self, name="$(name)_internal_equality")
        self.fixed_gain_node = FixedGainNode(A, self, name="$(name)_internal_gain")
        Edge(self.equality_node.interfaces[2], self.fixed_gain_node.in1) # Internal edge

        if use_composite_update_rules
            # Initialize the composite node interfaces belonging to the composite node itself.
            self.interfaces[1] = Interface(self)
            self.interfaces[2] = Interface(self)
            self.interfaces[3] = Interface(self)
        else
            # Initialize the interfaces as references to the internal node interfaces.
            self.interfaces[1] = self.equality_node.interfaces[1]
            self.interfaces[2] = self.equality_node.interfaces[3]
            self.interfaces[3] = self.fixed_gain_node.out
        end
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
        self.out = self.interfaces[3]
        return self
    end
end
GainEqualityCompositeNode(A::Array; args...) = GainEqualityCompositeNode(A, true; args...)
GainEqualityCompositeNode(; args...) = GainEqualityCompositeNode([1.0], true; args...)

############################################
# GaussianMessage methods
############################################

# Rule set for backward propagation ({in2,out}-->in1 or {in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y
backwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x - V_x * A' * inv(V_y + A * V_x * A') * A * V_x
backwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, V_x::Array{T, 2}, m_y::Array{T, 1}, V_y::Array{T, 2}) = m_x + V_x * A' * inv(V_y + A * V_x * A') * (m_y - A * m_x)

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GainEqualityCompositeNode,
                            inbound_messages::Array{GaussianMessage, 1})
    # Calculate an outbound message based on the inbound_messages array and the node function.
    # This function is not exported, and is only meant for internal use.
    # inbound_messages is indexed with the interface ids of the node.
    # inbound_messages[outbound_interface_id] should be #undef to indicate that the inbound message on this interface is not relevant.

    if isdefined(inbound_messages, outbound_interface_id)
        warn("The inbound message on the outbound interface is not undefined ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
    end
    if !node.use_composite_update_rules
        error("You can't call updateNodeMessage!() on a composite node ($(typeof(node)) $(node.name) interface $(outbound_interface_id)) that uses its internal graph to pass messages. Use calculateMessage!(node.interface) instead.")
    end

    if outbound_interface_id == 3
        # Forward message
        # We don't have a shortcut rule for this one, so we use the internal nodes to calculate the outbound msg
        # msg_i = inbound message on interface i, msg_out is always the calculated outbound message.
        msg_1 = inbound_messages[1]
        msg_2 = inbound_messages[2]
        msg_out = GaussianMessage()

        # First pass through the internal equality node
        inbound_messages = Array(GaussianMessage, 3)
        inbound_messages[1] = msg_1
        inbound_messages[3] = msg_2
        msg_temp = updateNodeMessage!(2, node.equality_node, inbound_messages)

        # Then go forward through the internal gain node
        inbound_messages = Array(GaussianMessage, 2)
        inbound_messages[1] = msg_temp
        msg_out = updateNodeMessage!(2, node.fixed_gain_node, inbound_messages)
    elseif outbound_interface_id == 1 || outbound_interface_id == 2
        # Backward messages
        msg_out = GaussianMessage()
        msg_3 = inbound_messages[3]
        msg_in = inbound_messages[outbound_interface_id == 1 ? 2 : 1] # the other input interface

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_3.xi != nothing && msg_3.W != nothing && msg_in.xi != nothing && msg_in.W != nothing
            msg_out.m = nothing
            msg_out.V = nothing
            msg_out.W = backwardGainEqualityWRule(node.A, msg_in.W, msg_3.W)
            msg_out.xi = backwardGainEqualityXiRule(node.A, msg_in.xi, msg_3.xi)
        elseif msg_3.m != nothing && msg_3.V != nothing && msg_in.m != nothing && msg_in.V != nothing
            msg_out.m = backwardGainEqualityMRule(node.A, msg_in.m, msg_in.V, msg_3.m, msg_3.V)
            msg_out.V = backwardGainEqualityVRule(node.A, msg_in.V, msg_3.V)
            msg_out.W = nothing
            msg_out.xi = nothing
        elseif msg_3.m != nothing && msg_3.W != nothing && msg_in.m != nothing && msg_in.W != nothing
            # TODO: Not very efficient!
            msg_out.m = backwardGainEqualityMRule(node.A, msg_in.m, inv(msg_in.W), msg_3.m, inv(msg_3.W))
            msg_out.V = nothing
            msg_out.W = backwardGainEqualityWRule(node.A, msg_in.W, msg_3.W)
            msg_out.xi = nothing
        else
            # Fallback: convert inbound messages to (xi,W) parametrization and then use efficient rules
            ensureXiWParametrization!(msg_in)
            ensureXiWParametrization!(msg_3)
            msg_out.m = nothing
            msg_out.V = nothing
            msg_out.W = backwardGainEqualityWRule(node.A, msg_in.W, msg_3.W)
            msg_out.xi = backwardGainEqualityXiRule(node.A, msg_in.xi, msg_3.xi)
        end
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end
    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end

