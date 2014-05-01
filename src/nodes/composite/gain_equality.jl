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
    #use_composite_update_rules::Bool
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

    function GainEqualityCompositeNode(A::Array; args...)
        name = "#undef"
        for (key, val) in args
            if key==:name
                name=val
            end
        end

        # Deepcopy A to avoid an unexpected change of the input argument A. Ensure that A is a matrix.
        A_copy = ensureMatrix(deepcopy(A))
        # Build composite node
        self = new(A_copy, name, Array(Interface, 3))
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
GainEqualityCompositeNode(; args...) = GainEqualityCompositeNode([1.0]; args...)

