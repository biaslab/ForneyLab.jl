############################################
# GainEqualityCompositeNode
############################################
# Description:
#   Composite gain-equality node.
#
#        _________
#    in1 |       | out2
#   -----|->[=]--|------>
#        |   |   |
#        |   v   |
#        |  [A]  |
#        |___|___|
#            | out1
#            v 
#
#   in1 = A*out1 = out2
#   Example:
#       GainEqualityCompositeNode([1.0]; name="my_node")
#   Gain A is fixed and defined through the constructor.
#
# Interface ids, (names) and supported message types:
#   1. in1:
#       GaussianMessage
#   2. out1:
#       GaussianMessage
#   3. out2:
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
    out1::Interface
    out2::Interface
    A_inv::Array{Float64, 2} # holds pre-computed inv(A) if possible

    function GainEqualityCompositeNode(A::Array, use_composite_update_rules::Bool; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end

        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        A_copy = ensureMatrix(deepcopy(A))
        # Build composite node
        self = new(A_copy, use_composite_update_rules, name, Array(Interface, 3))
        # Define the internals of the composite node
        self.equality_node = EqualityNode()
        self.fixed_gain_node = FixedGainNode(A_copy)
        Edge(self.equality_node.interfaces[2], self.fixed_gain_node.in1) # Internal edge

        # Set pointers of the composite node interfaces to the internal nodes
        self.interfaces[1] = self.equality_node.interfaces[1]
        self.interfaces[2] = self.fixed_gain_node.out
        self.interfaces[3] = self.equality_node.interfaces[3]
        # Init named interface handles
        self.in1 = self.interfaces[1]
        self.out1 = self.interfaces[2]
        self.out2 = self.interfaces[3]
        # Pointer to precomputed inv(A), needed for quick calculate rules
        self.A_inv = self.fixed_gain_node.A_inv
        return self
    end
end
GainEqualityCompositeNode(A::Array; args...) = GainEqualityCompositeNode(A, true; args...)
GainEqualityCompositeNode(use_composite_update_rules::Bool; args...) = GainEqualityCompositeNode([1.0], use_composite_update_rules; args...)
GainEqualityCompositeNode(; args...) = GainEqualityCompositeNode([1.0], true; args...)

############################################
# GaussianMessage methods
############################################

# Rule set for forward propagation, from: Korl (2005), "A Factor graph approach to signal modelling, system identification and filtering", Table 4.1
forwardGainEqualityWRule{T<:Number}(A::Array{T, 2}, W_x::Array{T, 2}, W_y::Array{T, 2}) = W_x + A' * W_y * A
forwardGainEqualityXiRule{T<:Number}(A::Array{T, 2}, xi_x::Array{T, 1}, xi_y::Array{T, 1}) = xi_x + A' * xi_y

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

    if node.use_composite_update_rules
        # Use the shortcut update rules
        if outbound_interface_id == 3
            # Forward message
            msg_in = inbound_messages[1]
            msg_in1 = deepcopy(msg_in)
            msg_out = inbound_messages[2]
            msg_out1 = deepcopy(msg_out)
            msg_out2 = GaussianMessage()

            # Select parameterization
            # Order is from least to most computationally intensive
            if msg_in1.xi != nothing && msg_in1.W != nothing && msg_out1.xi != nothing && msg_out1.W != nothing
                msg_out2.m = nothing
                msg_out2.V = nothing
                msg_out2.W = forwardGainEqualityWRule(node.A, msg_in1.W, msg_out1.W)
                msg_out2.xi = forwardGainEqualityXiRule(node.A, msg_in1.xi, msg_out1.xi)
            else
                # TODO: other parametrizations
                error("No other parametrizations are implemented yet")
            end
        else
            # TODO: other interfaces; 3 and 1 are the samy by symmetry
            error("No other interfaces are implemented yet")
        end

        # Set the outbound message
        return node.interfaces[outbound_interface_id].message = msg_out2
    else
        # Pass the message through the internals

        # TODO: this approach needs the graph to be defined, bacause it calls calculateMessage!(). That's not the idea behind the updateNodeMessage! function?
        warn("Explicit passing in composite nodes is still a little shabby")
               
        return node.interfaces[outbound_interface_id].message = calculateMessage!(node.interfaces[outbound_interface_id]) # Just call calculateMessage on the interface, this interface is a pointer to the subnode interface
    end
end

