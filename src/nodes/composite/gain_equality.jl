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

# Rule set for backward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
backwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
backwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y

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
        # Backward message; msg_i = inbound message on interface i, msg_out is always the calculated outbound message.
        msg_1 = inbound_messages[1]
        msg_2 = inbound_messages[2]
        msg_out = GaussianMessage()

        # Select parameterization
        # Order is from least to most computationally intensive
        if msg_1.xi != nothing && msg_1.W != nothing && msg_2.xi != nothing && msg_2.W != nothing
            msg_out.m = nothing
            msg_out.V = nothing
            msg_out.W = backwardGainEqualityWRule(node.A, msg_1.W, msg_2.W)
            msg_out.xi = backwardGainEqualityXiRule(node.A, msg_1.xi, msg_2.xi)
        else
            # TODO: other parametrizations
            error("No other parametrizations are implemented yet")
        end
        # Set the outbound message
        return node.interfaces[outbound_interface_id].message = msg_out
    else
        # TODO: other interfaces; 3 and 1 are the same by symmetry
        error("No other interfaces are implemented yet")
    end

end

