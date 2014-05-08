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

type GainAdditionCompositeNode <: Node
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
    addition_node::EqualityNode
    fixed_gain_node::FixedGainNode
    # Helper fields filled by constructor
    in1::Interface
    in2::Interface
    out::Interface

    function GainAdditionCompositeNode(A::Array, use_composite_update_rules::Bool; args...)
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
            Edge(self.fixed_gain_node.out, self.addition_node.in1) # Internal edge
            # Set pointers of the composite node interfaces to the internal nodes
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
GainAdditionCompositeNode(A::Array; args...) = GainAdditionCompositeNode(A, true; args...)
GainAdditionCompositeNode(; args...) = GainAdditionCompositeNode([1.0], true; args...)

############################################
# GaussianMessage methods
############################################

# TODO