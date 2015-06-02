# This file contains intgration helper functions for constructing and validating context graphs

#############
# Mocks
#############

type MockNode <: Node
    # MockNode is a node with an arbitrary function, that when created
    # initiates a message on all its interfaces.
    # Interface 1 is named :out
    id::Symbol
    interfaces::Array{Interface, 1}
    i::Dict{Symbol, Interface}

    function MockNode(num_interfaces::Int=1; id=ForneyLab.generateNodeId(MockNode))
        self = new(id, Array(Interface, num_interfaces), Dict{Symbol, Interface}())
        !haskey(current_graph.n, id) ? current_graph.n[id] = self : error("Node id $(id) already present")

        for interface_id = 1:num_interfaces
            self.interfaces[interface_id] = Interface(self)
        end
        
        self.i[:out] = self.interfaces[1]

        return(self)
    end
end
function MockNode(message::Message, num_interfaces::Int=1; kwargs...)
    self = MockNode(num_interfaces; kwargs...)
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
    node1 = MockNode()
    node2 = MockNode()
    return node1, node2
end

function initializePairOfTerminalNodes(d1::ProbabilityDistribution=GaussianDistribution(), d2::ProbabilityDistribution=GaussianDistribution())
    # Two connected TerminalNodes
    #
    # [T]-->[T]

    FactorGraph()
    t1 = TerminalNode(d1)
    t2 = TerminalNode(d2)
    Edge(t1.i[:out], t2.i[:out])
    t1.i[:out].message = Message(d1)
    t2.i[:out].message = Message(d2)
    return t1, t2
end

function initializePairOfNodes(; A=[1.0], msg_gain_1=Message(DeltaDistribution(2.0)), msg_gain_2=Message(DeltaDistribution(3.0)), msg_terminal=Message(DeltaDistribution(1.0)))
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
    node1 = TerminalNode(DeltaDistribution(reshape([3.0], 1, 1)))
    node2 = FixedGainNode([2.0])
    node3 = FixedGainNode([2.0])
    Edge(node1.i[:out], node2.i[:in])
    Edge(node2.i[:out], node3.i[:in])
    Edge(node3.i[:out], MockNode().i[:out])
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
    driver      = FixedGainNode(A, id=:driver)
    inhibitor   = FixedGainNode(B, id=:inhibitor)
    noise       = TerminalNode(GaussianDistribution(m=noise_m, V=noise_V), id=:noise)
    add         = AdditionNode(id=:adder)
    Edge(add.i[:out], inhibitor.i[:in])
    Edge(inhibitor.i[:out], driver.i[:in])
    Edge(driver.i[:out], add.i[:in1])
    Edge(noise.i[:out], add.i[:in2])
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
    Edge(c1.i[:out], add.i[:in1])
    Edge(c2.i[:out], add.i[:in2])
    Edge(add.i[:out], equ.interfaces[1])
    Edge(c3.i[:out], equ.interfaces[2])
    Edge(c3.i[:out], equ.interfaces[2])
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
    t1 = TerminalNode(id=:t1)
    a1 = FixedGainNode(id=:a1)
    g1 = GaussianNode(form=:moment, id=:g1)
    t2 = TerminalNode(id=:t2)
    add1 = AdditionNode(id=:add1)
    g2 = GaussianNode(id=:g2, V=0.1)
    Edge(t1.i[:out], a1.i[:in])
    Edge(a1.i[:out], g1.i[:mean])
    Edge(t2.i[:out], g1.i[:variance], InverseGammaDistribution)
    Edge(g1.i[:out], add1.i[:in1])
    Edge(add1.i[:out], g2.i[:mean])
    Edge(g2.i[:out], add1.i[:in2])
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
    t1 = TerminalNode(id=:t1)
    a1 = FixedGainNode(id=:a1)
    g1 = GaussianNode(form=:moment, id=:g1)
    t2 = TerminalNode(id=:t2)
    t3 = TerminalNode(id=:t3)
    Edge(t1.i[:out], a1.i[:in])
    Edge(a1.i[:out], g1.i[:mean])
    Edge(t2.i[:out], g1.i[:variance], InverseGammaDistribution)
    Edge(g1.i[:out], t3.i[:out])
    return (t1, a1, g1, t2, t3)
end

function initializeGaussianFactoringGraph()
    # [T]<--[A]<--[N]

    FactorGraph()
    t = TerminalNode(id=:t)
    gain = FixedGainNode(id=:gain)
    gauss = GaussianNode(m=1.0, V=0.5, id=:gauss)
    Edge(gauss.i[:out], gain.i[:in])
    Edge(gain.i[:out], t.i[:out])
    return (t, gain, gauss)
end

function initializeSimpleFactoringGraph()
    # [T]<--[A]<--[T]

    FactorGraph()
    t1 = TerminalNode(id=:t1)
    gain = FixedGainNode(id=:gain)
    t2 = TerminalNode(id=:t2)
    Edge(t2.i[:out], gain.i[:in])
    Edge(gain.i[:out], t1.i[:out])
    return (t1, gain, t2)
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
            Edge(MockNode(msg).i[:out], add_node.interfaces[interface_count])
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
            Edge(MockNode(msg).i[:out], eq_node.interfaces[interface_count])
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
    #     c_node
    #    -------
    #    |     |
    # |--|-[+]-|--|
    #    |     |
    #    | ... |

    FactorGraph()
    c_node = GainAdditionNode([1.0], false)
    node = TerminalNode()
    return(c_node, node)
end

function initializeGainAdditionNode(A::Array, msgs::Array{Any})
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
    gac_node = GainAdditionNode(A)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().i[:out], gac_node.interfaces[interface_count])
        else
            Edge(MockNode(msg).i[:out], gac_node.interfaces[interface_count])
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
    #     c_node
    #    -------
    #    |     |
    # |--|-[=]-|--|
    #    |     |
    #    | ... |

    FactorGraph()
    c_node = GainEqualityNode([1.0], false)
    node = TerminalNode()
    return(c_node, node)
end

function initializeGainEqualityNode(A::Array, msgs::Array{Any})
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
    gec_node = GainEqualityNode(A)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().i[:out], gec_node.interfaces[interface_count])
        else
            Edge(MockNode(msg).i[:out], gec_node.interfaces[interface_count])
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
    node = GaussianNode(form=:precision)
    edges = Array(Edge, 3)
    edges[1] = Edge(MockNode().i[:out], node.i[:mean], GaussianDistribution)
    edges[1].tail.message = Message(GaussianDistribution())
    edges[1].head.message = Message(GaussianDistribution())
    edges[2] = Edge(MockNode().i[:out], node.i[:precision], GammaDistribution)
    edges[2].tail.message = Message(GammaDistribution())
    edges[2].head.message = Message(GammaDistribution())
    edges[3] = Edge(node.i[:out], MockNode().i[:out], GaussianDistribution)
    edges[3].tail.message = Message(GaussianDistribution())
    if y_type == Float64
        edges[3].head.message = Message(DeltaDistribution(1.0))
    elseif y_type == GaussianDistribution
        edges[3].head.message = Message(GaussianDistribution())
    else
        error("Can't handle y_type $(y_type)")
    end
    return (node, edges)
end

function initializeBufferGraph()
    # Background
    g = FactorGraph()
    node_t1 = TerminalNode()
    node_t2 = TerminalNode()
    e = Edge(node_t1, node_t2)
    return (node_t1, node_t2, e, g)
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
        g_node = GaussianNode(; id=symbol("g_node_$(section)"), form=:precision)
        m_eq_node = EqualityNode(; id=symbol("m_eq_$(section)")) # Equality node chain for mean
        gam_eq_node = EqualityNode(; id=symbol("s_eq_$(section)")) # Equality node chain for variance
        y_node = TerminalNode(GaussianDistribution(m=y[section], V=tiny()), id=symbol("c_obs_$(section)")) # Observed y values are stored in terminal node
        g_nodes[section] = g_node
        m_eq_nodes[section] = m_eq_node
        gam_eq_nodes[section] = gam_eq_node
        y_nodes[section] = y_node
        q_y_edges[section] = Edge(g_node.i[:out], y_node.i[:out], GaussianDistribution)
        q_m_edges[section] = Edge(m_eq_node.interfaces[3], g_node.i[:mean], GaussianDistribution)
        q_gam_edges[section] = Edge(gam_eq_node.interfaces[3], g_node.i[:precision], GammaDistribution)
        if section > 1 # Connect sections
            Edge(m_eq_nodes[section-1].interfaces[2], m_eq_nodes[section].interfaces[1], GaussianDistribution)
            Edge(gam_eq_nodes[section-1].interfaces[2], gam_eq_nodes[section].interfaces[1], GammaDistribution)
        end
    end
    # Attach beginning and end nodes
    m_0 = TerminalNode(GaussianDistribution(m=0.0, V=100.0)) # Prior
    gam_0 = TerminalNode(GammaDistribution(a=0.01, b=0.01)) # Prior
    m_N = TerminalNode(vague(GaussianDistribution))
    gam_N = TerminalNode(vague(GammaDistribution))
    Edge(m_0.i[:out], m_eq_nodes[1].interfaces[1])
    Edge(gam_0.i[:out], gam_eq_nodes[1].interfaces[1])
    Edge(m_eq_nodes[end].interfaces[2], m_N.i[:out])
    Edge(gam_eq_nodes[end].interfaces[2], gam_N.i[:out])

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

    g_node = GaussianNode(; id=:g_node, form=:precision)
    m_eq_node = EqualityNode(; id=:m_eq_node) # Equality node chain for mean
    gam_eq_node = EqualityNode(; id=:gam_eq_node) # Equality node chain for variance
    y_node = TerminalNode(GaussianDistribution(), id=:y_node) # Observed y values are stored in terminal node
    y_edge = Edge(g_node.i[:out], y_node.i[:out], GaussianDistribution)
    m_edge = Edge(m_eq_node.interfaces[3], g_node.i[:mean], GaussianDistribution)
    gam_edge = Edge(gam_eq_node.interfaces[3], g_node.i[:precision], GammaDistribution)

    # Attach beginning and end nodes
    m_0_node = TerminalNode(GaussianDistribution(m=0.0, V=100.0); id=:m_0_node) # Prior
    gam_0_node = TerminalNode(GammaDistribution(a=0.01, b=0.01); id=:gam_0_node) # Prior
    m_N_node = TerminalNode(vague(GaussianDistribution); id=:m_N_node) # Neutral 'one' message
    gam_N_node = TerminalNode(vague(GammaDistribution); id=:gam_N_node) # Neutral 'one' message
    m_0_eq_edge = Edge(m_0_node.i[:out], m_eq_node.interfaces[1], GaussianDistribution)
    gam_0_eq_edge = Edge(gam_0_node.i[:out], gam_eq_node.interfaces[1], GammaDistribution)
    m_N_eq_edge = Edge(m_eq_node.interfaces[2], m_N_node.i[:out], GaussianDistribution)
    gam_N_eq_edge = Edge(gam_eq_node.interfaces[2], gam_N_node.i[:out], GammaDistribution)

    return (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge)
end


#############
# Validations
#############

function ==(x::ScheduleEntry, y::ScheduleEntry)
    if is(x, y) return true end
    ((x.interface == y.interface) && (x.message_calculation_rule == y.message_calculation_rule)) || (return false)
    (isdefined(x, :post_processing) == isdefined(y, :post_processing)) || (return false)
    if isdefined(x, :post_processing)
        (x.post_processing == y.post_processing) || (return false)
    end
    return true
end

function testInterfaceConnections(node1::FixedGainNode, node2::TerminalNode)
    # Helper function for node comparison

    # Check that nodes are properly connected
    @fact typeof(node1.interfaces[1].message.payload) <: DeltaDistribution => true
    @fact typeof(node2.interfaces[1].message.payload) <: DeltaDistribution => true
    @fact mean(node1.interfaces[1].message.payload) => 2.0
    @fact mean(node2.interfaces[1].message.payload) => 1.0
    @fact mean(node1.interfaces[1].partner.message.payload) => 1.0
    @fact mean(node2.interfaces[1].partner.message.payload) => 2.0
    # Check that pointers are initiatized correctly
    @fact mean(node1.i[:out].message.payload) => 3.0
    @fact mean(node2.i[:out].message.payload) => 1.0
    @fact mean(node1.i[:in].partner.message.payload) => 1.0
    @fact mean(node2.i[:out].partner.message.payload) => 2.0
end

function validateOutboundMessage(node::Node, outbound_interface_id::Int, inbound_messages::Array, correct_outbound_value::ProbabilityDistribution, update_function::Function=ForneyLab.sumProduct!)
    (rule, msg) = update_function(node, outbound_interface_id, inbound_messages...)
    @fact node.interfaces[outbound_interface_id].message => msg
    @fact node.interfaces[outbound_interface_id].message.payload => correct_outbound_value

    return node.interfaces[outbound_interface_id].message
end