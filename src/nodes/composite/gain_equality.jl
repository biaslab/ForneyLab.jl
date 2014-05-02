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

type GainEqualityCompositeNode <: Node
    # Basic node properties.
    # Use_composite_update_rule is a flag for indicating use of shortcut update rules. When set to true, the internals are ignored and the updated
    # outgoing message is calculated via the provided update rules. When set to false, messages are passed through the internals of the composite node.
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
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function GainEqualityCompositeNode(A::Array, use_composite_update_rules::Bool; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end

        if use_composite_update_rules
            # Build node, use composite rules without building internal graph.
            # This initializes the composite node interfaces belonging to the composite node itself.
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            self = new(ensureMatrix(deepcopy(A)), true, name, Array(Interface, 3))
            # Try to precompute inv(A)
            try
                self.A_inv = inv(self.A)
            catch
                warn("The specified multiplier for ", string(typeof(self)), " ", self.name, " is not invertible. This might cause problems. Please check if this is what you really want.")
            end
            # Initialize interfaces
            self.interfaces[1] = Interface(self)
            self.interfaces[2] = Interface(self)
            self.interfaces[3] = Interface(self)
        else
            # Build internal graph, ignore composite rules.
            # This initializes the interfaces as references to the internal node interfaces.
            self = new(A, false, name, Array(Interface, 3))
            # Define the internals of the composite node
            self.equality_node = EqualityNode()
            self.fixed_gain_node = FixedGainNode(A)
            Edge(self.equality_node.interfaces[2], self.fixed_gain_node.in1) # Internal edge
            # Set pointers of the composite node interfaces to the internal nodes
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

# Rule set for forward propagation
forwardGainEqualityWRule{T<:Number}(A_inv::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = A_inv' * (W_x + W_y) * A_inv
forwardGainEqualityXiRule{T<:Number}(A_inv::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = A_inv' * (xi_x + xi_y)
forwardGainEqualityVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = A * V_x * pinv(V_x + V_y) * V_y * A'
forwardGainEqualityMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, W_x::Array{T, 2}, m_y::Array{T, 1}, W_y::Array{T, 2}) = A * pinv(W_x + W_y) * (W_x * m_x + W_y * m_y)

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
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

    # Use the shortcut update rules
    if outbound_interface_id == 3
        # Forward message; msg_i = inbound message on interface i, msg_out is always the calculated outbound message.
        msg_1 = inbound_messages[1]
        msg_2 = inbound_messages[2]
        msg_out = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.xi != nothing && msg_1.W != nothing && msg_2.xi != nothing && msg_2.W != nothing && isRoundedPosDef(node.A) && isRoundedPosDef(msg_1.W + msg_2.W) && isdefined(node, :A_inv)
            msg_out.m = nothing
            msg_out.V = nothing
            msg_out.W = forwardGainEqualityWRule(node.A_inv, msg_1.W, msg_2.W)
            msg_out.xi = forwardGainEqualityXiRule(node.A_inv, msg_1.xi, msg_2.xi)
        elseif msg_1.m != nothing && msg_1.W != nothing && msg_2.m != nothing && msg_2.W != nothing && isdefined(node, :A_inv)
            msg_out.m = forwardGainEqualityMRule(node.A, msg_1.m, msg_1.W, msg_2.m, msg_2.W)
            msg_out.V = nothing
            msg_out.W = forwardGainEqualityWRule(node.A_inv, msg_1.W, msg_2.W)
            msg_out.xi = nothing
        elseif msg_1.m != nothing && msg_1.V != nothing && msg_2.m != nothing && msg_2.V != nothing
            # TODO: Not very efficient!
            msg_out.m = forwardGainEqualityMRule(node.A, msg_1.m, inv(msg_1.V), msg_2.m, inv(msg_2.V))
            msg_out.V = forwardGainEqualityVRule(node.A, msg_1.V, msg_2.V)
            msg_out.W = nothing
            msg_out.xi = nothing
        else
            # Alternative parameterization not caught by the above rules
            warn("Parameterization unknown; using xi, W parameterization instead. ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
            ensureXiWParametrization!(msg_1)
            ensureXiWParametrization!(msg_2)
            # Check preconditions
            if isRoundedPosDef(node.A) && isRoundedPosDef(msg_1.W + msg_2.W)
                msg_out.m = nothing
                msg_out.V = nothing
                msg_out.W = forwardGainEqualityWRule(node.A_inv, msg_1.W, msg_2.W)
                msg_out.xi = forwardGainEqualityXiRule(node.A_inv, msg_1.xi, msg_2.xi)
            else
                error("Cannot calculate outbound message on ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
            end
        end
        # Set the outbound message
        return node.interfaces[outbound_interface_id].message = msg_out
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
            # Alternative parameterization not caught by the above rules
            warn("Parameterization unknown; using xi, W parameterization instead ($(typeof(node)) $(node.name) interface $(outbound_interface_id))")
            ensureXiWParametrization!(msg_in)
            ensureXiWParametrization!(msg_3)
            msg_out.m = nothing
            msg_out.V = nothing
            msg_out.W = backwardGainEqualityWRule(node.A, msg_in.W, msg_3.W)
            msg_out.xi = backwardGainEqualityXiRule(node.A, msg_in.xi, msg_3.xi)
        end
        # Set the outbound message
        return node.interfaces[outbound_interface_id].message = msg_out
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end

end

