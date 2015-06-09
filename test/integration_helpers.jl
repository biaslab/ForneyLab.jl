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
        !haskey(current_graph.n, id) || error("Node id $(id) already present")
        current_graph.n[id] = self
 
        for interface_index = 1:num_interfaces
            self.interfaces[interface_index] = Interface(self)
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

    g = FactorGraph()
    MockNode(id=:node1)
    MockNode(id=:node2)

    return g
end

function initializePairOfTerminalNodes(d1::ProbabilityDistribution=GaussianDistribution(), d2::ProbabilityDistribution=GaussianDistribution())
    # Two connected TerminalNodes
    #
    # [T]-->[T]

    g = FactorGraph()
    TerminalNode(d1, id=:t1)
    TerminalNode(d2, id=:t2)
    Edge(n(:t1).i[:out], n(:t2).i[:out])
    n(:t1).i[:out].message = Message(d1)
    n(:t2).i[:out].message = Message(d2)

    return g
end

function initializePairOfNodes(; A=[1.0], msg_gain_1=Message(DeltaDistribution(2.0)), msg_gain_2=Message(DeltaDistribution(3.0)), msg_terminal=Message(DeltaDistribution(1.0)))
    # Helper function for initializing an unconnected pair of nodes
    #
    # |--[A]--|
    #
    # |--[C]

    g = FactorGraph()
    FixedGainNode(A, id=:node1)
    n(:node1).interfaces[1].message = msg_gain_1
    n(:node1).interfaces[2].message = msg_gain_2
    TerminalNode(msg_terminal.payload, id=:node2)
    n(:node2).interfaces[1].message = msg_terminal

    return g
end

function initializeChainOfNodes()
    # Chain of three nodes
    #
    #  1     2     3
    # [T]-->[A]-->[B]-->[M]

    g = FactorGraph()
    TerminalNode(DeltaDistribution(reshape([3.0], 1, 1)), id=:node1)
    FixedGainNode([2.0], id=:node2)
    FixedGainNode([2.0], id=:node3)
    Edge(n(:node1).i[:out], n(:node2).i[:in])
    Edge(n(:node2).i[:out], n(:node3).i[:in])
    Edge(n(:node3).i[:out], MockNode().i[:out])

    return g
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

    g = FactorGraph()
    FixedGainNode(A, id=:driver)
    FixedGainNode(B, id=:inhibitor)
    TerminalNode(GaussianDistribution(m=noise_m, V=noise_V), id=:noise)
    AdditionNode(id=:add)
    Edge(n(:add).i[:out], n(:inhibitor).i[:in])
    Edge(n(:inhibitor).i[:out], n(:driver).i[:in])
    Edge(n(:driver).i[:out], n(:add).i[:in1])
    Edge(n(:noise).i[:out], n(:add).i[:in2])

    return g
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

    g = FactorGraph()
    TerminalNode(GaussianDistribution(), id=:c1)
    TerminalNode(GaussianDistribution(), id=:c2)
    TerminalNode(GaussianDistribution(m=-2.0, V=3.0), id=:c3)
    AdditionNode(id=:add)
    EqualityNode(id=:equ)
    # Edges from left to right
    Edge(n(:c1).i[:out], n(:add).i[:in1])
    Edge(n(:c2).i[:out], n(:add).i[:in2])
    Edge(n(:add).i[:out], n(:equ).interfaces[1])
    Edge(n(:c3).i[:out], n(:equ).interfaces[2])
    Edge(n(:c3).i[:out], n(:equ).interfaces[2])

    return g
end

function initializeFactoringGraph()
    # Set up a graph to test factorize function
    #             [T]
    #              |     -----
    #              v     v   |
    # [T]-->[A]-->[N]-->[+] [N]
    #                    |   ^
    #                    -----

    g = FactorGraph()
    TerminalNode(id=:t1)
    FixedGainNode(id=:a1)
    GaussianNode(form=:moment, id=:g1)
    TerminalNode(id=:t2)
    AdditionNode(id=:add1)
    GaussianNode(id=:g2, V=0.1)
    Edge(n(:t1).i[:out], n(:a1).i[:in])
    Edge(n(:a1).i[:out], n(:g1).i[:mean])
    Edge(n(:t2).i[:out], n(:g1).i[:variance], InverseGammaDistribution)
    Edge(n(:g1).i[:out], n(:add1).i[:in1])
    Edge(n(:add1).i[:out], n(:g2).i[:mean])
    Edge(n(:g2).i[:out], n(:add1).i[:in2])

    return g
end

function initializeFactoringGraphWithoutLoop()
    # Set up a graph to test factorize function
    #             [T]
    #              |
    #              v
    # [T]-->[A]-->[N]-->[T]
    #
    #

    g = FactorGraph()
    TerminalNode(id=:t1)
    FixedGainNode(id=:a1)
    GaussianNode(form=:moment, id=:g1)
    TerminalNode(id=:t2)
    TerminalNode(id=:t3)
    Edge(n(:t1).i[:out], n(:a1).i[:in])
    Edge(n(:a1).i[:out], n(:g1).i[:mean])
    Edge(n(:t2).i[:out], n(:g1).i[:variance], InverseGammaDistribution)
    Edge(n(:g1).i[:out], n(:t3).i[:out])

    return g
end

function initializeGaussianFactoringGraph()
    # [T]<--[A]<--[N]

    g = FactorGraph()
    TerminalNode(id=:t)
    FixedGainNode(id=:gain)
    GaussianNode(m=1.0, V=0.5, id=:gauss)
    Edge(n(:gauss).i[:out], n(:gain).i[:in])
    Edge(n(:gain).i[:out], n(:t).i[:out])

    return g
end

function initializeSimpleFactoringGraph()
    # [T]<--[A]<--[T]

    g = FactorGraph()
    TerminalNode(id=:t1)
    FixedGainNode(id=:gain)
    TerminalNode(id=:t2)
    Edge(n(:t2).i[:out], n(:gain).i[:in])
    Edge(n(:gain).i[:out], n(:t1).i[:out])

    return g
end

function initializeAdditionNode(msgs::Array{Any})
    # Set up an addition node and prepare the messages
    # A MockNode is connected for each argument message.
    #
    # [M]-->[+]<--[M]
    #        |

    g = FactorGraph()
    AdditionNode(id=:add_node)
    interface_count = 1
    for msg in msgs
        if msg != nothing
            Edge(MockNode(msg).i[:out], n(:add_node).interfaces[interface_count])
        end
        interface_count += 1
    end

    return g
end

function initializeEqualityNode(msgs::Array{Any})
    # Set up an equality node and prepare the messages
    # A MockNode is connected for each argument message
    #
    # [M]-->[=]<--[M] (as many incoming edges as length(msgs))
    #        |

    g = FactorGraph()
    EqualityNode(length(msgs), id=:eq_node)
    interface_count = 1
    for msg in msgs
        if msg!=nothing
            Edge(MockNode(msg).i[:out], n(:eq_node).interfaces[interface_count])
        end
        interface_count += 1
    end
    return g
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

    g = FactorGraph()
    GainAdditionNode([1.0], id=:c_node)
    TerminalNode(id=:node)

    return g
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

    g = FactorGraph()
    GainAdditionNode(A, id=:gac_node)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().i[:out], n(:gac_node).interfaces[interface_count])
        else
            Edge(MockNode(msg).i[:out], n(:gac_node).interfaces[interface_count])
        end
        interface_count += 1
    end

    return g
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

    g = FactorGraph()
    GainEqualityNode([1.0], id=:c_node)
    TerminalNode(id=:node)

    return g
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

    g = FactorGraph()
    GainEqualityNode(A, id=:gec_node)
    interface_count = 1
    for msg in msgs
        if msg == nothing
            Edge(MockNode().i[:out], n(:gec_node).interfaces[interface_count])
        else
            Edge(MockNode(msg).i[:out], n(:gec_node).interfaces[interface_count])
        end
        interface_count += 1
    end

    return g
end

function initializeGaussianNode(; y_type::DataType=Float64)
    # Initialize a Gaussian node
    #
    #    mean   precision
    #  [M]-->[N]<--[M]
    #         |
    #         v y
    #        [M]

    g = FactorGraph()
    GaussianNode(form=:precision, id=:node)
    Edge(MockNode().i[:out], n(:node).i[:mean], GaussianDistribution, id=:edge1)
    e(:edge1).tail.message = Message(GaussianDistribution())
    e(:edge1).head.message = Message(GaussianDistribution())
    Edge(MockNode().i[:out], n(:node).i[:precision], GammaDistribution, id=:edge2)
    e(:edge2).tail.message = Message(GammaDistribution())
    e(:edge2).head.message = Message(GammaDistribution())
    Edge(n(:node).i[:out], MockNode().i[:out], GaussianDistribution, id=:edge3)
    e(:edge3).tail.message = Message(GaussianDistribution())
    if y_type == Float64
        e(:edge3).head.message = Message(DeltaDistribution(1.0))
    elseif y_type == GaussianDistribution
        e(:edge3).head.message = Message(GaussianDistribution())
    else
        error("Can't handle y_type $(y_type)")
    end

    return g
end

function initializeBufferGraph()
    g = FactorGraph()
    TerminalNode(id=:node_t1)
    TerminalNode(id=:node_t2)
    Edge(n(:node_t1), n(:node_t2), id=:e)

    return g
end

function initializeWrapGraph()
    g = FactorGraph()
    TerminalNode(id=:t1)
    TerminalNode(id=:t2)
    TerminalNode(id=:t3)
    TerminalNode(id=:t4)
    AdditionNode(id=:add1)
    AdditionNode(id=:add2)
    Edge(n(:t1), n(:add1).i[:in1])
    Edge(n(:t2), n(:add1).i[:in2])
    Edge(n(:t3), n(:add2).i[:in1])
    Edge(n(:t4), n(:add2).i[:in2])
    Edge(n(:add1).i[:out], n(:add2).i[:out])

    return g
end

function initializeCompositeGraph()
    # Build the internal graph
    g = FactorGraph()
    t_constant = TerminalNode(3.0)
    t_in = TerminalNode(id=:in)
    t_out = TerminalNode(id=:out)
    a = AdditionNode(id=:adder)
    Edge(t_in, a.i[:in1])
    Edge(t_constant, a.i[:in2])
    Edge(a.i[:out], t_out)
    
    return (g, t_in, t_out)
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

    g = FactorGraph()
    # Initial settings
    n_samples = length(y) # Number of observed samples

    # Build graph
    for sec=1:n_samples
        GaussianNode(form=:precision, id=s(:g,sec))
        EqualityNode(id=s(:m_eq,sec)) # Equality node chain for mean
        EqualityNode(id=s(:gam_eq,sec)) # Equality node chain for precision
        TerminalNode(GaussianDistribution(m=y[sec], V=tiny()), id=s(:y,sec)) # Observed y values are stored in terminal node
        Edge(n(s(:g,sec)).i[:out], n(s(:y,sec)).i[:out], GaussianDistribution, id=s(:q_y,sec))
        Edge(n(s(:m_eq,sec)).i[3], n(s(:g,sec)).i[:mean], GaussianDistribution, id=s(:q_m,sec))
        Edge(n(s(:gam_eq,sec)).i[3], n(s(:g,sec)).i[:precision], GammaDistribution, id=s(:q_gam,sec))
        if sec > 1 # Connect sections
            Edge(n(s(:m_eq,sec-1)).i[2], n(s(:m_eq,sec)).i[1], GaussianDistribution)
            Edge(n(s(:gam_eq,sec-1)).i[2], n(s(:gam_eq,sec)).i[1], GammaDistribution)
        end
    end
    # Attach beginning and end nodes
    TerminalNode(GaussianDistribution(m=0.0, V=100.0), id=:m0) # Prior
    TerminalNode(GammaDistribution(a=0.01, b=0.01), id=:gam0) # Prior
    TerminalNode(vague(GaussianDistribution), id=:mN)
    TerminalNode(vague(GammaDistribution), id=:gamN)
    Edge(n(:m0).i[:out], n(:m_eq1).i[1])
    Edge(n(:gam0).i[:out], n(:gam_eq1).i[1])
    Edge(n(s(:m_eq,n_samples)).i[2], n(:mN).i[:out])
    Edge(n(s(:gam_eq,n_samples)).i[2], n(:gamN).i[:out])

    return g
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

    g = FactorGraph()

    GaussianNode(id=:g, form=:precision)
    EqualityNode(id=:m_eq) # Equality node chain for mean
    EqualityNode(id=:gam_eq) # Equality node chain for variance
    TerminalNode(GaussianDistribution(), id=:y) # Observed y values are stored in terminal node
    Edge(n(:g).i[:out], n(:y).i[:out], GaussianDistribution, id=:y)
    Edge(n(:m_eq).i[3], n(:g).i[:mean], GaussianDistribution, id=:m)
    Edge(n(:gam_eq).i[3], n(:g).i[:precision], GammaDistribution, id=:gam)

    # Attach beginning and end nodes
    TerminalNode(GaussianDistribution(m=0.0, V=100.0), id=:m0) # Prior
    TerminalNode(GammaDistribution(a=0.01, b=0.01), id=:gam0) # Prior
    TerminalNode(vague(GaussianDistribution), id=:mN) # Neutral 'one' message
    TerminalNode(vague(GammaDistribution), id=:gamN) # Neutral 'one' message
    Edge(n(:m0).i[:out], n(:m_eq).i[1], GaussianDistribution, id=:m0)
    Edge(n(:gam0).i[:out], n(:gam_eq).i[1], GammaDistribution, id=:gam0)
    Edge(n(:m_eq).i[2], n(:mN).i[:out], GaussianDistribution, id=:mN)
    Edge(n(:gam_eq).i[2], n(:gamN).i[:out], GammaDistribution, id=:gamN)

    return g
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

function validateOutboundMessage(node::Node, outbound_interface_index::Int, inbound_messages::Array, correct_outbound_value::ProbabilityDistribution, update_function::Function=ForneyLab.sumProduct!)
    (rule, msg) = update_function(node, outbound_interface_index, inbound_messages...)
    @fact node.interfaces[outbound_interface_index].message => msg
    @fact node.interfaces[outbound_interface_index].message.payload => correct_outbound_value

    return node.interfaces[outbound_interface_index].message
end