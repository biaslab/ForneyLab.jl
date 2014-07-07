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
    edge = Edge(MockNode(GaussianMessage()).out, g_node.in1)
    if msgs[1] != nothing
        edge.marginal = msgs[1]
    end
    edge = Edge(MockNode(GammaMessage()).out, g_node.in2)
    if msgs[2] != nothing
        edge.marginal = msgs[2]
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
        g_node = GaussianNode(true) # Variational flag set to true, so updateNodeMessage knows what formula to use
        m_eq_node = EqualityNode() # Equality node chain for mean
        gam_eq_node = EqualityNode() # Equality node chain for variance
        obs_node = ConstantNode(GeneralMessage(y_observations[section])) # Observed y values are stored as constant node values
        g_nodes[section] = g_node
        m_eq_nodes[section] = m_eq_node
        gam_eq_nodes[section] = gam_eq_node
        obs_nodes[section] = obs_node
        y_edges[section] = Edge(g_node.out, obs_node.out)
        q_m_edges[section] = Edge(m_eq_node.interfaces[3], g_node.in1)
        q_gam_edges[section] = Edge(gam_eq_node.interfaces[3], g_node.in2)
        # Preset uninformative ('one') messages
        setMarginal!(q_gam_edges[section], uninformative(GammaMessage))
        setMarginal!(q_m_edges[section], uninformative(GaussianMessage))
        if section > 1 # Connect sections
            Edge(m_eq_nodes[section-1].interfaces[2], m_eq_nodes[section].interfaces[1])
            Edge(gam_eq_nodes[section-1].interfaces[2], gam_eq_nodes[section].interfaces[1])
        end
    end
    # Attach beginning and end nodes
    m_prior = ConstantNode(GaussianMessage(m=[0.0], V=[100.0])) # Prior
    s_prior = ConstantNode(GammaMessage(a=0.01, b=0.01)) # Prior
    c_m = ConstantNode(uninformative(GaussianMessage)) # Neutral 'one' message
    c_gam = ConstantNode(uninformative(GammaMessage)) # Neutral 'one' message
    Edge(m_prior.out, m_eq_nodes[1].interfaces[1])
    Edge(s_prior.out, gam_eq_nodes[1].interfaces[1])
    Edge(m_eq_nodes[end].interfaces[2], c_m.out)
    Edge(gam_eq_nodes[end].interfaces[2], c_gam.out)

    return(g_nodes, obs_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges)
end

function initializeLinearCompositeNode(msgs::Array{Message})
    # Set up a linear composite node and prepare the marginals
    # A MockNode is connected for each argument message
    #
    #        [a_in][b_in] [s_in]
    #            |  |       |
    #       |-----------------|
    #       |    |  |       | |
    #       |    |  -->[N]<-- |
    #       |    |      |     |
    #       |    v      v     |
    #[in1]--|-->[a]--->[+]----|-->[out]
    #       |                 |
    #       |-----------------|

    lin_node = LinearCompositeNode(true)
    for intf = 1:5
        edge = Edge(MockNode(GaussianMessage()).out, lin_node.interfaces[intf])
        if msgs[intf] != nothing
            edge.marginal = msgs[intf]
        end
    end    
    return lin_node
end

function initializeLinearCompositeNodeChain()
    # Set up a chain of Gaussian nodes for regression parameter estimation
    #
    #  [a_0]-->[=]------ ... ->[=]------->[C_a]
    #          |               |        
    #  [b_0]---|>[=]---- ... --|>[=]----->[C_b]
    #          |  |            |  |     
    #  [s_0]---|--|>[=]- ... --|--|>[=]-->[C_s]
    #          |  |  |         |  |  |  
    #          v  v  v         v  v  v  
    #         |-------|       |-------| 
    #         |   L   |       |   L   | 
    #         |-------|       |-------| 
    #           ^   |           ^   |   
    #           |   v           |   v   
    #          x_1 y_1         x_n y_n

    # prepare samples
    true_s = 2.0
    true_a = 3.0
    true_b = 5.0
    n_samples = 20
    x = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0]
    y = [6.1811923622357625,7.917496269084679,11.286102016681964,14.94255088702814,16.82264686442818,19.889355802073506,23.718253510300688,28.18105443765643,27.72075206943362,32.15921446069328,34.97262678800721,38.86444301740928,40.79138365100076,45.84963364094473,47.818481172238165,51.51027620022872,52.623019301773,53.91583839744505,58.14426361122961,59.895517438500164]
    # x = [0.0:(n_samples-1.0)]
    # y = true_a*x + true_b + sqrt(true_s)*randn(n_samples)

    # show some stuff
    # plot(x,y)
    # println("True slope a: $(true_a)")
    # println("True offset b: $(true_b)")
    # println("True noise variance s: $(true_s)")
    # println("-------")
    # pars = linreg(x, y)
    # println("Sample slope a: $(pars[2])")
    # println("Sample offset b: $(pars[1])")
    # println("Sample variance s: $(var((x - pars[1])/pars[2]))")
    # println("-------")

    # Pre-assign arrays for later reference
    lin_nodes = Array(LinearCompositeNode, n_samples)
    a_eq_nodes = Array(EqualityNode, n_samples)
    b_eq_nodes = Array(EqualityNode, n_samples)
    s_eq_nodes = Array(EqualityNode, n_samples)
    x_nodes = Array(ConstantNode, n_samples)
    y_nodes = Array(ConstantNode, n_samples)
    a_eq_edges = Array(Edge, n_samples)
    b_eq_edges = Array(Edge, n_samples)
    s_eq_edges = Array(Edge, n_samples)
    x_edges = Array(Edge, n_samples)
    y_edges = Array(Edge, n_samples)

    # Build graph
    for section=1:n_samples
        lin_node = LinearCompositeNode()
        a_eq_node = EqualityNode()
        b_eq_node = EqualityNode()
        s_eq_node = EqualityNode()
        x_node = ConstantNode(GaussianMessage(m = [x[section]], V = [0.0001]))
        y_node = ConstantNode(GaussianMessage(m = [y[section]], V = [0.0001]))
        # Save to array
        lin_nodes[section] = lin_node
        a_eq_nodes[section] = a_eq_node
        b_eq_nodes[section] = b_eq_node
        s_eq_nodes[section] = s_eq_node
        x_nodes[section] = x_node
        y_nodes[section] = y_node
        # Connect section within
        a_eq_edges[section] = Edge(a_eq_node.interfaces[3], lin_node.a_in)
        b_eq_edges[section] = Edge(b_eq_node.interfaces[3], lin_node.b_in)
        s_eq_edges[section] = Edge(s_eq_node.interfaces[3], lin_node.s_in)
        x_edges[section] = Edge(x_node.out, lin_node.in1)
        y_edges[section] = Edge(lin_node.out, y_node.out)
        # Preset marginals
        setMarginal!(a_eq_edges[section], uninformative(GaussianMessage)) # uninformative
        setMarginal!(b_eq_edges[section], uninformative(GaussianMessage))
        setMarginal!(s_eq_edges[section], uninformative(InverseGammaMessage))
        setMarginal!(x_edges[section], GaussianMessage(m = [x[section]], V = [0.0001])) # samples
        setMarginal!(y_edges[section], GaussianMessage(m = [y[section]], V = [0.0001]))

        if section > 1 # Connect sections
            Edge(a_eq_nodes[section-1].interfaces[2], a_eq_nodes[section].interfaces[1])
            Edge(b_eq_nodes[section-1].interfaces[2], b_eq_nodes[section].interfaces[1])
            Edge(s_eq_nodes[section-1].interfaces[2], s_eq_nodes[section].interfaces[1])
        end
    end
    # Attach beginning and end nodes
    a_0 = ConstantNode(GaussianMessage(m=[0.0], V=[100.0])) # priors
    b_0 = ConstantNode(GaussianMessage(m=[0.0], V=[100.0]))
    s_0 = ConstantNode(InverseGammaMessage(a=0.001, b=0.001))
    C_a = ConstantNode(uninformative(GaussianMessage)) # uninformative
    C_b = ConstantNode(uninformative(GaussianMessage))
    C_s = ConstantNode(uninformative(InverseGammaMessage))
    # connect
    Edge(a_0, a_eq_nodes[1].interfaces[1])
    Edge(b_0, b_eq_nodes[1].interfaces[1])
    Edge(s_0, s_eq_nodes[1].interfaces[1])
    Edge(a_eq_nodes[end].interfaces[2], C_a.out)
    Edge(b_eq_nodes[end].interfaces[2], C_b.out)
    Edge(s_eq_nodes[end].interfaces[2], C_s.out)

    return(lin_nodes, a_eq_nodes, b_eq_nodes, s_eq_nodes, a_eq_edges, b_eq_edges, s_eq_edges, x_edges, y_edges)
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
