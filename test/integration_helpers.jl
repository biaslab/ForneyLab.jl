# This file contains ntgration helper functions for constructing and validating context graphs

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
        self.name = "#undefined"
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

function initializePairOfNodes(; A=[1.0], msg_gain_1=Message(2.0), msg_gain_2=Message(3.0), msg_const=Message(1.0))
    # Helper function for initializing an unconnected pair of nodes
    #     
    # |--[A]--|
    #
    # |--[C]

    node1 = FixedGainNode(A)
    node1.interfaces[1].message = msg_gain_1
    node1.interfaces[2].message = msg_gain_2
    node2 = ConstantNode(msg_const.value)
    node2.interfaces[1].message = msg_const
    return node1, node2
end

function initializeFixedGainNode()    
    # Helper function for initializing a fixed gain node
    #     
    # |--[A]--|

    node = FixedGainNode()
    node.interfaces[1].message = Message(GaussianDistribution())
    node.interfaces[2].message = Message(GaussianDistribution())
    return node
end

function initializeChainOfNodes()
    # Chain of three nodes
    #
    #  1     2     3
    # [C]-->[A]-->[B]-|

    node1 = ConstantNode(3.0)
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
    noise       = ConstantNode(GaussianDistribution(m=noise_m, V=noise_V), name="noise")
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
    c1 = ConstantNode(GaussianDistribution())
    c2 = ConstantNode(GaussianDistribution())
    c3 = ConstantNode(GaussianDistribution(m=[-2.], V=[3.]))
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

function initializeVariationalGaussianNode(msgs::Array{Any})
    # Set up a Gaussian node and prepare the messages/marginals
    # A MockNode is connected for each argument message
    #
    #         [M] (mean)
    #          |
    #          v  out
    #   [M]-->[N]----->
    #  (prec.)
    #                

    g_node = GaussianNode(true)
    edge = Edge(MockNode(Message(GaussianDistribution())).out, g_node.in1)
    if msgs[1] != nothing
        edge.marginal = msgs[1].value
    end
    edge = Edge(MockNode(Message(GammaDistribution())).out, g_node.in2)
    if msgs[2] != nothing
        edge.marginal = msgs[2].value
    end
    Edge(g_node.out, MockNode(msgs[3]).out)
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

function initializeGaussianNodeChain()
    # Set up a chain of Gaussian nodes for joint mean-precision estimation
    #
    # [gam_prior]-------->[=]---------->[=]--->    -->[C_gam]
    #                      |             |   etc...
    # [m_prior]-->[=]---------->[=]--------->      -->[C_m]
    #          q(m)| q(gam)|     |       |
    #              -->[N]<--     -->[N]<--
    #                  |             |
    #                  v y_1         v y_2

    # Initial settings
    # Fixed observations drawn from N(5.0, 2.0)
    y_observations = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
    n_samples = length(y_observations) # Number of observed samples

    # Pre-assign arrays for later reference
    g_nodes = Array(GaussianNode, n_samples)
    m_eq_nodes = Array(EqualityNode, n_samples)
    gam_eq_nodes = Array(EqualityNode, n_samples)
    obs_nodes = Array(ConstantNode, n_samples)
    q_gam_edges = Array(Edge, n_samples)
    q_m_edges = Array(Edge, n_samples)
    y_edges = Array(Edge, n_samples)

    # Build graph
    for section=1:n_samples
        g_node = GaussianNode(true; name="g_node_$(section)") # Variational flag set to true, so updateNodeMessage knows what formula to use
        m_eq_node = EqualityNode(; name="m_eq_$(section)") # Equality node chain for mean
        gam_eq_node = EqualityNode(; name="s_eq_$(section)") # Equality node chain for variance
        obs_node = ConstantNode(y_observations[section], name="c_obs_$(section)") # Observed y values are stored as constant node values
        g_nodes[section] = g_node
        m_eq_nodes[section] = m_eq_node
        gam_eq_nodes[section] = gam_eq_node
        obs_nodes[section] = obs_node
        y_edges[section] = Edge(g_node.out, obs_node.out)
        q_m_edges[section] = Edge(m_eq_node.interfaces[3], g_node.in1)
        q_gam_edges[section] = Edge(gam_eq_node.interfaces[3], g_node.in2)
        # Preset uninformative ('one') messages
        setMarginal!(q_gam_edges[section], uninformative(GammaDistribution))
        setMarginal!(q_m_edges[section], uninformative(GaussianDistribution))
        if section > 1 # Connect sections
            Edge(m_eq_nodes[section-1].interfaces[2], m_eq_nodes[section].interfaces[1])
            Edge(gam_eq_nodes[section-1].interfaces[2], gam_eq_nodes[section].interfaces[1])
        end
    end
    # Attach beginning and end nodes
    m_prior = ConstantNode(GaussianDistribution(m=[0.0], V=[100.0])) # Prior
    s_prior = ConstantNode(GammaDistribution(a=0.01, b=0.01)) # Prior
    c_m = ConstantNode(uninformative(GaussianDistribution)) # Neutral 'one' message
    c_gam = ConstantNode(uninformative(GammaDistribution)) # Neutral 'one' message
    Edge(m_prior.out, m_eq_nodes[1].interfaces[1])
    Edge(s_prior.out, gam_eq_nodes[1].interfaces[1])
    Edge(m_eq_nodes[end].interfaces[2], c_m.out)
    Edge(gam_eq_nodes[end].interfaces[2], c_gam.out)

    return(g_nodes, obs_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges)
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
