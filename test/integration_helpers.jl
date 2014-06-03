# This file contains ntgration helper functions for constructing and validating context graphs

#############
# Backgrounds
#############

function initializePairOfNodes(; A=[1.0], msg_gain_1::Message=GeneralMessage(2.0), msg_gain_2=GeneralMessage(3.0), msg_const=GeneralMessage(1.0))
	# Helper function for initializing an unconnected pair of nodes
	#     
	# |--[A]--|
	#
	# |--[C]

    node1 = FixedGainNode(A)
    node1.interfaces[1].message = msg_gain_1
    node1.interfaces[2].message = msg_gain_2
    node2 = ConstantNode()
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
