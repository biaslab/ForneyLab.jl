# This file contains ntgration helper functions for constructing and validating context graphs

#############
# Mocks
#############

type MockNode <: Node
	# MockNode is a node with an arbitrary function, that when created
	# initiates a valid message on its only interface 'out'
	interfaces::Array{Interface, 1}
	out::Interface
    name::ASCIIString
	function MockNode(num_interfaces::Int=1)
        self = new(Array(Interface, num_interfaces))
        for interface_id = 1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        self.out = self.interfaces[1]
        self.name = "#undefined"
        return(self)
	end
end
function MockNode(message::Message, num_interfaces::Int=1)
    self = MockNode(num_interfaces)
    for interface in self.interfaces
        interface.message = message
        interface.message_valid = true
    end
    return(self)
end

#############
# Backgrounds
#############

function initializePairOfMockNodes()
    # Two unconnected MockNodes
    #
    # [M]--|
    #
    # [M]--|
    
    return MockNode(), MockNode()
end

function initializePairOfNodes(; A=[1.0], msg_gain_1=GeneralMessage(2.0), msg_gain_2=GeneralMessage(3.0), msg_const=GeneralMessage(1.0))
	# Helper function for initializing an unconnected pair of nodes
	#     
	# |--[A]--|
	#
	# |--[C]

    node1 = FixedGainNode(A)
    node1.interfaces[1].message = msg_gain_1
    node1.interfaces[2].message = msg_gain_2
    node2 = ConstantNode(msg_const)
    node2.interfaces[1].message = msg_const
    return node1, node2
end

function initializeFixedGainNode()    
	# Helper function for initializing a fixed gain node
	#     
	# |--[A]--|

	node = FixedGainNode()
    node.interfaces[1].message = GaussianMessage()
    node.interfaces[2].message = GaussianMessage()
    return node
end

function initializeChainOfNodes()
	# Chain of three nodes
	#
	#  1     2     3
	# [C]-->[A]-->[B]-|

    node1 = ConstantNode(GeneralMessage(3.0))
    node2 = FixedGainNode([2.0])
    node3 = FixedGainNode([2.0])
    Edge(node1.out, node2.in1)
    Edge(node2.out, node3.in1)
    return (node1, node2, node3)
end

function initializeLoopyGraph(; A=[2.0], B=[0.5], noise_m=[1.0], noise_V=[0.1])
    # Set up a loopy graph
    #    (driver)
    #   -->[A]---
    #   |       |
    #   |      [+]<-[N]
    #   |       |
    #   ---[B]<--
    #  (inhibitor)

    driver      = FixedGainNode(A, name="driver")
    inhibitor   = FixedGainNode(B, name="inhibitor")
    noise       = ConstantNode(GaussianMessage(m=noise_m, V=noise_V), name="noise")
    add         = AdditionNode(name="adder")
    Edge(add.out, inhibitor.in1)
    Edge(inhibitor.out, driver.in1)
    Edge(driver.out, add.in1)
    Edge(noise.out, add.in2)
    return (driver, inhibitor, noise, add)
end

function initializeTreeGraph()
	# Set up some tree graph
    #
    #          (c2)
    #           |
    #           v
    # (c1)---->[+]---->[=]----->
    #                   ^    y
    #                   |
    #                  (c3)
    #
    c1 = ConstantNode(GaussianMessage())
    c2 = ConstantNode(GaussianMessage())
    c3 = ConstantNode(GaussianMessage(m=[-2.], V=[3.]))
    add = AdditionNode()
    equ = EqualityNode()
    # Edges from left to right
    Edge(c1.out, add.in1)
    Edge(c2.out, add.in2)
    Edge(add.out, equ.interfaces[1])
    Edge(c3.out, equ.interfaces[2])
    Edge(c3.out, equ.interfaces[2])
    return (c1, c2, c3, add, equ)
end

function initializeAdditionNode(msgs::Array{Any})
	# Set up an addition node and prepare the messages
	# A MockNode is connected for each argument message.
	#
	# [M]-->[+]<--[M]	
	#        |        

	add_node = AdditionNode()
	interface_count = 1
	for msg=msgs
		if msg!=nothing
			Edge(MockNode(msg).out, add_node.interfaces[interface_count])
		end
		interface_count += 1
	end
	return add_node
end

function initializeEqualityNode(msgs::Array{Any})
	# Set up an equality node and prepare the messages
	# A MockNode is connected for each argument message
	#
	# [M]-->[=]<--[M] (as many incoming edges as length(msgs))	
	#        |        

	eq_node = EqualityNode(length(msgs))
	interface_count = 1
	for msg=msgs
		if msg!=nothing
			Edge(MockNode(msg).out, eq_node.interfaces[interface_count])
		end
		interface_count += 1
	end
	return eq_node
end

function initializeFixedGainNode(A::Array, msgs::Array{Any})
	# Set up a fixed gain node and prepare the messages
	# A MockNode is connected for each argument message
	#
	# [M]-->[A]--|	
	#                

	fg_node = FixedGainNode(A)
	interface_count = 1
	for msg=msgs
		if msg!=nothing
			Edge(MockNode(msg).out, fg_node.interfaces[interface_count])
		end
		interface_count += 1
	end
	return fg_node
end

function initializeGaussianNode(msgs::Array{Any})
    # Set up a Gaussian node and prepare the messages
    # A MockNode is connected for each argument message
    #
    #         [M] (mean)
    #          |
    #          v  out
    #   [M]-->[N]----->
    #  (prec.)
    #                

    g_node = GaussianNode()
    interface_count = 1
    for msg=msgs
        if msg!=nothing
            Edge(MockNode(msg).out, g_node.interfaces[interface_count])
        end
        interface_count += 1
    end
    return g_node
end

function initializeConstantAndGainAddNode()
    # Initialize some nodes
    #
    #    node
    #    [N]--| 
    #       out
    #
    #       c_node
    #    ------------
    #    |          |
    # |--| |--[+]-| |--|
    # in2| in2 |    |
    #    |    ...   |

    c_node = GainAdditionCompositeNode([1.0], false)
    node = ConstantNode()
    return(c_node, node)
end

function initializeGainAdditionCompositeNode(A::Array, use_composite_update_rules::Bool, msgs::Array{Any})
	# Set up a gain addition node and prepare the messages
	# A MockNode is connected for each argument message
	#
    #           [M]
    #            | in1
    #            |
    #        ____|____
    #        |   v   |
    #        |  [A]  |
    #        |   |   |
    #    in2 |   v   | out
    #[M]-----|->[+]--|---->
    #        |_______|

	gac_node = GainAdditionCompositeNode(A, use_composite_update_rules)
	interface_count = 1
	for msg=msgs
		if msg!=nothing
			Edge(MockNode(msg).out, gac_node.interfaces[interface_count])
		end
		interface_count += 1
	end
	return gac_node
end

function initializeConstantAndGainEqNode()
    # Initialize some nodes
    #
    #    node
    #    [N]--| 
    #       out
    #
    #       c_node
    #    ------------
    #    |          |
    # |--| |--[=]-| |--|
    # in1| in1 |    |
    #    |    ...   |

    c_node = GainEqualityCompositeNode([1.0], false)
    node = ConstantNode()
    return(c_node, node)
end

function initializeGainEqualityCompositeNode(A::Array, use_composite_update_rules::Bool, msgs::Array{Any})
	# Set up a gain equality node and prepare the messages
	# A MockNode is connected for each argument message
	#
    #         _________
    #     in1 |       | in2
    # [M]-----|->[=]<-|-----[M]
    #         |   |   |
    #         |   v   |
    #         |  [A]  |
    #         |___|___|
    #             | out
    #             v

	gec_node = GainEqualityCompositeNode(A, use_composite_update_rules)
	interface_count = 1
	for msg=msgs
		if msg!=nothing
			Edge(MockNode(msg).out, gec_node.interfaces[interface_count])
		end
		interface_count += 1
	end
	return gec_node
end

#############
# Validations
#############

function testInterfaceConnections(node1::FixedGainNode, node2::ConstantNode)
	# Helper function for node comparison

    # Check that nodes are properly connected
    @fact node1.interfaces[1].message.value => 2.0
    @fact node2.interfaces[1].message.value => 1.0
    @fact node1.interfaces[1].partner.message.value => 1.0
    @fact node2.interfaces[1].partner.message.value => 2.0
    # Check that pointers are initiatized correctly
    @fact node1.out.message.value => 3.0
    @fact node2.out.message.value => 1.0
    @fact node1.in1.partner.message.value => 1.0
    @fact node2.out.partner.message.value => 2.0
end
