############################################
# GainAdditionCompositeNode
############################################
# Description:
#   Composite gain-addition node.
#
#            | in1
#            |
#        ____|____
#        |   v   |
#        |  [A]  |
#        |   |   |
#    in2 |   v   | out
#   -----|->[+]--|---->
#        |_______|
#
#
#   out = A*in1 + in2
#   Example:
#       GainAdditionCompositeNode([1.0]; name="my_node")
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
export GainAdditionCompositeNode

type GainAdditionCompositeNode <: CompositeNode
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
    addition_node::AdditionNode
    fixed_gain_node::FixedGainNode
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface

    function GainAdditionCompositeNode(A::Array=[1.0], use_composite_update_rules::Bool=true, parent::Union(CompositeNode, Nothing)=nothing; args...)
        (name = getArgumentValue(args, :name))!=false || (name = "unnamed")

        if use_composite_update_rules
            # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
            # In case we don't use composite update rules, A is passed to the internal FixedGainNode.
            A = ensureMatrix(deepcopy(A))
        end
        self = new(A, use_composite_update_rules, name, parent, Array(Interface, 3))

        # Define the internals of the composite node
        self.addition_node = AdditionNode(self, name="$(name)_internal_addition")
        self.fixed_gain_node = FixedGainNode(A, self, name="$(name)_internal_gain")
        Edge(self.fixed_gain_node.out, self.addition_node.in1) # Internal edge

        if use_composite_update_rules
            # Initialize the composite node interfaces belonging to the composite node itself.
            self.interfaces[1] = Interface(self)
            self.interfaces[2] = Interface(self)
            self.interfaces[3] = Interface(self)
        else
            # Initialize the interfaces as references to the internal node interfaces.
            self.interfaces[1] = self.fixed_gain_node.in1
            self.interfaces[2] = self.addition_node.in2
            self.interfaces[3] = self.addition_node.out
        end
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.in2 = self.interfaces[2]
        self.out = self.interfaces[3]
        return self
    end
end
GainAdditionCompositeNode(A::Array, parent::CompositeNode; args...) = GainAdditionCompositeNode(A, true, parent; args...)
GainAdditionCompositeNode(parent::CompositeNode; args...) = GainAdditionCompositeNode([1.0], true, parent; args...)

############################################
# GaussianMessage methods
############################################

# Rule set for forward propagation ({in1,in2}-->out)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainAdditionMRule{T<:Number}(A::Array{T, 2}, m_x::Array{T, 1}, m_y::Array{T, 1}) = m_x + A*m_y
forwardGainAdditionVRule{T<:Number}(A::Array{T, 2}, V_x::Array{T, 2}, V_y::Array{T, 2}) = V_x + A*V_y*A'
forwardGainAdditionWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x - W_x * A * inv(W_y+A'*W_x*A) * A' * W_x
forwardGainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}, W_x::Array{T, 2}, W_y::Array{T, 2}) = xi_x + W_x*A*inv(W_y+A'*W_x*A)*(xi_y-A'*xi_x)

# Rule set for backward propagation ({in1,out}-->in2)
# From: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardIn2GainAdditionMRule{T<:Number}(A::Array{T, 2}, m_y::Array{T, 1}, m_z::Array{T, 1}) = m_z - A*m_y
backwardIn2GainAdditionVRule{T<:Number}(A::Array{T, 2}, V_y::Array{T, 2}, V_z::Array{T, 2}) = V_z + A*V_y*A'
backwardIn2GainAdditionWRule{T<:Number}(A::Array{T, 2}, W_y::Array{T, 2}, W_z::Array{T, 2}) = W_z - W_z * A * inv(W_y+A'*W_z*A) * A' * W_z
backwardIn2GainAdditionXiRule{T<:Number}(A::Array{T, 2}, xi_y::Array{T, 1}, xi_z::Array{T, 1}, W_y::Array{T, 2}, W_z::Array{T, 2}) = xi_z - W_z*A*inv(W_y+A'*W_z*A)*(xi_y+A'*xi_z)

function updateNodeMessage!(outbound_interface_id::Int,
                            node::GainAdditionCompositeNode,
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
        # Forward message towards "out" interface
        # msg_i = inbound message on interface i, msg_out is always the calculated outbound message.
        msg_1 = inbound_messages[1]
        msg_2 = inbound_messages[2]
        msg_out = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.m != nothing && msg_1.V != nothing && msg_2.m != nothing && msg_2.V != nothing
            msg_out.m  = forwardGainAdditionMRule(node.A, msg_2.m, msg_1.m)
            msg_out.V  = forwardGainAdditionVRule(node.A, msg_2.V, msg_1.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        elseif msg_1.m != nothing && msg_1.W != nothing && msg_2.m != nothing && msg_2.W != nothing
            msg_out.m  = forwardGainAdditionMRule(node.A, msg_2.m, msg_1.m)
            msg_out.V  = nothing
            msg_out.W  = forwardGainAdditionWRule(node.A, msg_2.W, msg_1.W)
            msg_out.xi = nothing
        elseif msg_1.xi != nothing && msg_1.W != nothing && msg_2.xi != nothing && msg_2.W != nothing
            msg_out.m  = nothing
            msg_out.V  = nothing
            msg_out.W  = forwardGainAdditionWRule(node.A, msg_2.W, msg_1.W)
            msg_out.xi = forwardGainAdditionXiRule(node.A, msg_2.xi, msg_1.xi, msg_2.W, msg_1.W)
        elseif (msg_1.m != nothing && msg_1.V != nothing) || (msg_2.m != nothing && msg_2.V != nothing)
            # Fallback: at least one inbound msg is in (m,V) parametrization
            # Convert the other one to (m,V)
            ensureMVParametrization!(msg_1)
            ensureMVParametrization!(msg_2)
            msg_out.m  = forwardGainAdditionMRule(node.A, msg_2.m, msg_1.m)
            msg_out.V  = forwardGainAdditionVRule(node.A, msg_2.V, msg_1.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        elseif (msg_1.m != nothing && msg_1.W != nothing) || (msg_2.m != nothing && msg_2.W != nothing)
            # Fallback: at least one inbound msg is in (m,W) parametrization
            # Convert the other one to (m,W)
            ensureMWParametrization!(msg_1)
            ensureMWParametrization!(msg_2)
            msg_out.m  = forwardGainAdditionMRule(node.A, msg_2.m, msg_1.m)
            msg_out.V  = nothing
            msg_out.W  = forwardGainAdditionWRule(node.A, msg_2.W, msg_1.W)
            msg_out.xi = nothing
        else
            # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
            ensureMVParametrization!(msg_1)
            ensureMVParametrization!(msg_2)
            msg_out.m  = forwardGainAdditionMRule(node.A, msg_2.m, msg_1.m)
            msg_out.V  = forwardGainAdditionVRule(node.A, msg_2.V, msg_1.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        end
    elseif outbound_interface_id == 2
        # Backward message towards "in2" interface
        msg_1 = inbound_messages[1]
        msg_3 = inbound_messages[3]
        msg_out = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.m != nothing && msg_1.V != nothing && msg_3.m != nothing && msg_3.V != nothing
            msg_out.m  = backwardIn2GainAdditionMRule(node.A, msg_1.m, msg_3.m)
            msg_out.V  = backwardIn2GainAdditionVRule(node.A, msg_1.V, msg_3.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        elseif msg_1.m != nothing && msg_1.W != nothing && msg_3.m != nothing && msg_3.W != nothing
            msg_out.m  = backwardIn2GainAdditionMRule(node.A, msg_1.m, msg_3.m)
            msg_out.V  = nothing
            msg_out.W  = backwardIn2GainAdditionWRule(node.A, msg_1.W, msg_3.W)
            msg_out.xi = nothing
        elseif msg_1.xi != nothing && msg_1.W != nothing && msg_3.xi != nothing && msg_3.W != nothing
            msg_out.m  = nothing
            msg_out.V  = nothing
            msg_out.W  = backwardIn2GainAdditionWRule(node.A, msg_1.W, msg_3.W)
            msg_out.xi = backwardIn2GainAdditionXiRule(node.A, msg_1.xi, msg_3.xi, msg_1.W, msg_3.W)
        elseif (msg_1.m != nothing && msg_1.V != nothing) || (msg_3.m != nothing && msg_3.V != nothing)
            # Fallback: at least one inbound msg is in (m,V) parametrization
            # Convert the other one to (m,V)
            ensureMVParametrization!(msg_1)
            ensureMVParametrization!(msg_3)
            msg_out.m  = backwardIn2GainAdditionMRule(node.A, msg_1.m, msg_3.m)
            msg_out.V  = backwardIn2GainAdditionVRule(node.A, msg_1.V, msg_3.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        elseif (msg_1.m != nothing && msg_1.W != nothing) || (msg_3.m != nothing && msg_3.W != nothing)
            # Fallback: at least one inbound msg is in (m,W) parametrization
            # Convert the other one to (m,W)
            ensureMWParametrization!(msg_1)
            ensureMWParametrization!(msg_3)
            msg_out.m  = backwardIn2GainAdditionMRule(node.A, msg_1.m, msg_3.m)
            msg_out.V  = nothing
            msg_out.W  = backwardIn2GainAdditionWRule(node.A, msg_1.W, msg_3.W)
            msg_out.xi = nothing
        else
            # Fallback: if all else fails, convert everything to (m,V) and then use efficient rule
            ensureMVParametrization!(msg_1)
            ensureMVParametrization!(msg_3)
            msg_out.m  = backwardIn2GainAdditionMRule(node.A, msg_1.m, msg_3.m)
            msg_out.V  = backwardIn2GainAdditionVRule(node.A, msg_1.V, msg_3.V)
            msg_out.W  = nothing
            msg_out.xi = nothing
        end
    elseif outbound_interface_id == 1
        # Backward message towards "in1" interface
        # We don't have a shortcut rule for this one, so we use the internal nodes to calculate the outbound msg
        msg_2 = inbound_messages[2]
        msg_3 = inbound_messages[3]
        msg_out = GaussianMessage()

        # First pass through the internal addition node
        inbound_messages = Array(GaussianMessage, 3)
        inbound_messages[2] = msg_2
        inbound_messages[3] = msg_3
        msg_temp = updateNodeMessage!(1, node.addition_node, inbound_messages)

        # Then go backwards through the internal gain node
        inbound_messages = Array(GaussianMessage, 2)
        inbound_messages[2] = msg_temp
        msg_out = updateNodeMessage!(1, node.fixed_gain_node, inbound_messages)
    else
        error("Invalid outbound interface id $(outbound_interface_id), on $(typeof(node)) $(node.name).")
    end
    # Set the outbound message
    return node.interfaces[outbound_interface_id].message = msg_out
end
