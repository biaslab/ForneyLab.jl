# This file contains intgration helper functions for constructing and validating context graphs

#############
# Mocks
#############

type MockNode <: Node
    # MockNode is a node with an arbitrary function, that when created
    # initiates a message on all its interfaces.
    # Interface 1 is named "out"
    interfaces::Array{Interface, 1}
    out::Interface
    name::ASCIIString
    function MockNode(num_interfaces::Int=1)
        self = new(Array(Interface, num_interfaces))
        for interface_id = 1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        self.out = self.interfaces[1]
        self.name = ForneyLab.unnamedStr()
        return(self)
    end
end
function MockNode(message::Message, num_interfaces::Int=1)
    self = MockNode(num_interfaces)
    for interface in self.interfaces
        interface.message = message
    end
    return(self)
end
ForneyLab.isDeterministic(::MockNode) = false # Edge case, same as terminal node

#############
# Backgrounds
#############

function initializePairOfMockNodes()
    # Two unconnected MockNodes
    #
    # [M]--|
    #
    # [M]--|

    FactorGraph()
    return MockNode(), MockNode()
end

function initializePairOfNodes(; A=[1.0], msg_gain_1=Message(2.0), msg_gain_2=Message(3.0), msg_terminal=Message(1.0))
    # Helper function for initializing an unconnected pair of nodes
    #     
    # |--[A]--|
    #
    # |--[C]

    FactorGraph()
    node1 = FixedGainNode(A)
    node1.interfaces[1].message = msg_gain_1
    node1.interfaces[2].message = msg_gain_2
    node2 = TerminalNode(msg_terminal.payload)
    node2.interfaces[1].message = msg_terminal
    return node1, node2
end

function initializeChainOfNodes()
    # Chain of three nodes
    #
    #  1     2     3
    # [C]-->[A]-->[B]-->[M]

    FactorGraph()
    node1 = TerminalNode(reshape([3.0], 1, 1))
    node2 = FixedGainNode([2.0])
    node3 = FixedGainNode([2.0])
    Edge(node1.out, node2.in1, Array{Float64, 2})
    Edge(node2.out, node3.in1, Array{Float64, 2})
    Edge(node3.out, MockNode().out, Array{Float64, 2})
    return (node1, node2, node3)
end

function initializeLoopyGraph(; A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
    # Set up a loopy graph
    #    (driver)
    #   -->[A]---
    #   |       |
    #   |      [+]<-[N]
    #   |       |
    #   ---[B]<--
    #  (inhibitor)

    FactorGraph()
    driver      = FixedGainNode(A, name="driver")
    inhibitor   = FixedGainNode(B, name="inhibitor")
    noise       = TerminalNode(GaussianDistribution(m=noise_m, V=noise_V), name="noise")
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

    FactorGraph()
    c1 = TerminalNode(GaussianDistribution())
    c2 = TerminalNode(GaussianDistribution())
    c3 = TerminalNode(GaussianDistribution(m=-2.0, V=3.0))
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

function initializeFactoringGraph()
    # Set up a graph to test factorize function
    #             [T]    
    #              |     -----
    #              v     v   |
    # [T]-->[A]-->[N]-->[+] [N]
    #                    |   ^
    #                    -----

    FactorGraph()
    t1 = TerminalNode(name="t1")
    a1 = FixedGainNode(name="a1")
    g1 = GaussianNode(form="moment", name="g1")
    t2 = TerminalNode(name="t2")
    add1 = AdditionNode(name="add1")
    g2 = GaussianNode(name="g2", V=0.1)
    Edge(t1.out, a1.in1)
    Edge(a1.out, g1.mean)
    Edge(t2.out, g1.variance, InverseGammaDistribution)
    Edge(g1.out, add1.in1)
    Edge(add1.out, g2.mean)
    Edge(g2.out, add1.in2)
    return (t1, a1, g1, t2, add1, g2)
end

function initializeFactoringGraphWithoutLoop()
    # Set up a graph to test factorize function
    #             [T]    
    #              |     
    #              v     
    # [T]-->[A]-->[N]-->[T]
    #                    
    #                    

    FactorGraph()
    t1 = TerminalNode(name="t1")
    a1 = FixedGainNode(name="a1")
    g1 = GaussianNode(form="moment", name="g1")
    t2 = TerminalNode(name="t2")
    t3 = TerminalNode(name="t3")
    Edge(t1.out, a1.in1)
    Edge(a1.out, g1.mean)
    Edge(t2.out, g1.variance, InverseGammaDistribution)
    Edge(g1.out, t3.out)
    return (t1, a1, g1, t2, t3)
end


function initializeAdditionNode(msgs::Array{Any})
    # Set up an addition node and prepare the messages
    # A MockNode is connected for each argument message.
    #
    # [M]-->[+]<--[M]   
    #        |        

    FactorGraph()
    add_node = AdditionNode()
    interface_count = 1
    for msg in msgs
        if msg != nothing
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

    FactorGraph()
    eq_node = EqualityNode(length(msgs))
    interface_count = 1
    for msg in msgs
        if msg!=nothing
            Edge(MockNode(msg).out, eq_node.interfaces[interface_count])
        end
        interface_count += 1
    end
    return eq_node
end

function initializeTerminalAndGainAddNode()
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

    FactorGraph()
    c_node = GainAdditionCompositeNode([1.0], false)
    node = TerminalNode()
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

    FactorGraph()
    gac_node = GainAdditionCompositeNode(A, use_composite_update_rules)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().out, gac_node.interfaces[interface_count])
        else
            Edge(MockNode(msg).out, gac_node.interfaces[interface_count])
        end
        interface_count += 1
    end
    return gac_node
end

function initializeTerminalAndGainEqNode()
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

    FactorGraph()
    c_node = GainEqualityCompositeNode([1.0], false)
    node = TerminalNode()
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

    FactorGraph()
    gec_node = GainEqualityCompositeNode(A, use_composite_update_rules)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().out, gec_node.interfaces[interface_count])
        else
            Edge(MockNode(msg).out, gec_node.interfaces[interface_count])
        end
        interface_count += 1
    end
    return gec_node
end

function initializeGaussianNode(; y_type::DataType=Float64)
    # Initialize a Gaussian node
    #
    #    mean   precision
    #  [M]-->[N]<--[M]
    #         |
    #         v y
    #        [M]

    graph = FactorGraph()
    node = GaussianNode(form="precision")
    edges = Array(Edge, 3)
    edges[1] = Edge(MockNode().out, node.mean, GaussianDistribution)
    edges[1].tail.message = Message(GaussianDistribution())
    edges[1].head.message = Message(GaussianDistribution())
    edges[2] = Edge(MockNode().out, node.precision, GammaDistribution)
    edges[2].tail.message = Message(GammaDistribution())
    edges[2].head.message = Message(GammaDistribution())
    edges[3] = Edge(node.out, MockNode().out, GaussianDistribution)
    edges[3].tail.message = Message(GaussianDistribution())
    edges[3].head.message = Message(uninformative(y_type))

    # Set messages and marginals
    return (node, edges)
end

function initializeGaussianNodeChain(y::Array{Float64, 1})
    # Set up a chain of Gaussian nodes for mean-precision estimation
    #
    #     [gam_0]-------->[=]---------->[=]---->    -->[gam_N]
    #                      |             |     etc...
    #     [m_0]-->[=]---------->[=]------------>    -->[m_N]
    #          q(m)| q(gam)|     |       |
    #              -->[N]<--     -->[N]<--
    #                  | q(y)        |
    #                  v             v
    #                [y_1]         [y_2]

    FactorGraph()
    # Initial settings
    n_samples = length(y) # Number of observed samples

    # Pre-assign arrays for later reference
    g_nodes = Array(GaussianNode, n_samples)
    m_eq_nodes = Array(EqualityNode, n_samples)
    gam_eq_nodes = Array(EqualityNode, n_samples)
    y_nodes = Array(TerminalNode, n_samples)
    q_gam_edges = Array(Edge, n_samples)
    q_m_edges = Array(Edge, n_samples)
    q_y_edges = Array(Edge, n_samples)

    # Build graph
    for section=1:n_samples
        g_node = GaussianNode(; name="g_node_$(section)", form="precision") # Variational flag set to true, so updateNodeMessage knows what formula to use
        m_eq_node = EqualityNode(; name="m_eq_$(section)") # Equality node chain for mean
        gam_eq_node = EqualityNode(; name="s_eq_$(section)") # Equality node chain for variance
        y_node = TerminalNode(y[section], name="c_obs_$(section)") # Observed y values are stored in terminal node
        g_nodes[section] = g_node
        m_eq_nodes[section] = m_eq_node
        gam_eq_nodes[section] = gam_eq_node
        y_nodes[section] = y_node
        q_y_edges[section] = Edge(g_node.out, y_node.out, Float64)
        q_m_edges[section] = Edge(m_eq_node.interfaces[3], g_node.mean, GaussianDistribution)
        q_gam_edges[section] = Edge(gam_eq_node.interfaces[3], g_node.precision, GammaDistribution)
        if section > 1 # Connect sections
            Edge(m_eq_nodes[section-1].interfaces[2], m_eq_nodes[section].interfaces[1], GaussianDistribution)
            Edge(gam_eq_nodes[section-1].interfaces[2], gam_eq_nodes[section].interfaces[1], GammaDistribution)
        end
    end
    # Attach beginning and end nodes
    m_0 = TerminalNode(GaussianDistribution(m=0.0, V=100.0)) # Prior
    gam_0 = TerminalNode(GammaDistribution(a=0.01, b=0.01)) # Prior
    m_N = TerminalNode(uninformative(GaussianDistribution)) # Neutral 'one' message
    gam_N = TerminalNode(uninformative(GammaDistribution)) # Neutral 'one' message
    Edge(m_0.out, m_eq_nodes[1].interfaces[1])
    Edge(gam_0.out, gam_eq_nodes[1].interfaces[1], GammaDistribution)
    Edge(m_eq_nodes[end].interfaces[2], m_N.out)
    Edge(gam_eq_nodes[end].interfaces[2], gam_N.out, GammaDistribution)

    return (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges)
end

function initializeGaussianNodeChainForSvmp(y::Array{Float64, 1})
    # Set up a chain of Gaussian nodes for joint mean-precision estimation
    # through structured variational message passing
    #
    #     [gam_0]-------->[=]-->[1] gam_N
    #                      |
    #     [m_0]-->[=]---------->[1] m_N
    #              |       |
    #     q(m,gam) -->[N]<--
    # - - - - - - - - -|- - - - - - - - - -
    #                  |q(y)      
    #                  v       
    #                [y_1]

    # Batch estimation with multiple samples will intruduce cycles in the subgraph.
    # Therefore we implement a forward algorithm that uses forward estimation only.

    FactorGraph()

    g_node = GaussianNode(; name="g_node", form="precision")
    m_eq_node = EqualityNode(; name="m_eq_node") # Equality node chain for mean
    gam_eq_node = EqualityNode(; name="gam_eq_node") # Equality node chain for variance
    y_node = TerminalNode(GaussianDistribution(), name="y_node") # Observed y values are stored in terminal node
    y_edge = Edge(g_node.out, y_node.out, GaussianDistribution)
    m_edge = Edge(m_eq_node.interfaces[3], g_node.mean, GaussianDistribution)
    gam_edge = Edge(gam_eq_node.interfaces[3], g_node.precision, GammaDistribution)

    # Attach beginning and end nodes
    m_0_node = TerminalNode(GaussianDistribution(m=0.0, V=100.0)) # Prior
    gam_0_node = TerminalNode(GammaDistribution(a=0.01, b=0.01)) # Prior
    m_N_node = TerminalNode(uninformative(GaussianDistribution)) # Neutral 'one' message
    gam_N_node = TerminalNode(uninformative(GammaDistribution)) # Neutral 'one' message
    m_0_eq_edge = Edge(m_0_node.out, m_eq_node.interfaces[1], GaussianDistribution)
    gam_0_eq_edge = Edge(gam_0_node.out, gam_eq_node.interfaces[1], GammaDistribution)
    m_N_eq_edge = Edge(m_eq_node.interfaces[2], m_N_node.out, GaussianDistribution)
    gam_N_eq_edge = Edge(gam_eq_node.interfaces[2], gam_N_node.out, GammaDistribution)

    return (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge)
end


#############
# Validations
#############

function testInterfaceConnections(node1::FixedGainNode, node2::TerminalNode)
    # Helper function for node comparison

    # Check that nodes are properly connected
    @fact node1.interfaces[1].message.payload => 2.0
    @fact node2.interfaces[1].message.payload => 1.0
    @fact node1.interfaces[1].partner.message.payload => 1.0
    @fact node2.interfaces[1].partner.message.payload => 2.0
    # Check that pointers are initiatized correctly
    @fact node1.out.message.payload => 3.0
    @fact node2.out.message.payload => 1.0
    @fact node1.in1.partner.message.payload => 1.0
    @fact node2.out.partner.message.payload => 2.0
end

function validateOutboundMessage(node::Node, outbound_interface_id::Int, inbound_messages::Array, correct_outbound_value)
    msg = ForneyLab.updateNodeMessage!(node, outbound_interface_id, inbound_messages...)
    @fact node.interfaces[outbound_interface_id].message => msg
    @fact node.interfaces[outbound_interface_id].message.payload => correct_outbound_value

    return node.interfaces[outbound_interface_id].message
end
