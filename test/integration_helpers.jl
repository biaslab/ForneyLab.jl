# This file contains intgration helper functions for constructing and validating context graphs

import Base.==

#############
# Mocks
#############

type MockNode <: Node
    # MockNode is an arbitrary node without update functions
    # The last interface is called :out

    id::Symbol
    interfaces::Array{Interface, 1}
    i::Dict{Symbol, Interface}

    function MockNode(num_interfaces::Int=1; id=ForneyLab.generateNodeId(MockNode))
        self = new(id, Array(Interface, num_interfaces), Dict{Symbol, Interface}())
        addNode!(currentGraph(), self)

        for interface_index = 1:num_interfaces
            self.interfaces[interface_index] = Interface(self)
        end

        self.i[:out] = self.interfaces[end]

        return(self)
    end
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
    GainNode(gain=A, id=:node1)
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
    TerminalNode(DeltaDistribution(3.0), id=:node1)
    GainNode(gain=[2.0], id=:node2)
    GainNode(gain=[2.0], id=:node3)
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
    GainNode(gain=A, id=:driver)
    GainNode(gain=B, id=:inhibitor)
    TerminalNode(GaussianDistribution(m=noise_m, V=noise_V), id=:noise)
    AdditionNode(id=:add)
    Edge(n(:add).i[:out], n(:inhibitor).i[:in])
    Edge(n(:inhibitor).i[:out], n(:driver).i[:in])
    Edge(n(:driver).i[:out], n(:add).i[:in1])
    Edge(n(:noise).i[:out], n(:add).i[:in2])

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
    GainNode(gain=[1.0], id=:a1)
    GaussianNode(form=:variance, id=:g1)
    TerminalNode(id=:t2)
    AdditionNode(id=:add1)
    GaussianNode(id=:g2, V=0.1)
    Edge(n(:t1).i[:out], n(:a1).i[:in])
    Edge(n(:a1).i[:out], n(:g1).i[:mean])
    Edge(n(:t2).i[:out], n(:g1).i[:variance])
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
    GainNode(gain=[1.0], id=:a1)
    GaussianNode(form=:variance, id=:g1)
    TerminalNode(id=:t2)
    TerminalNode(id=:t3)
    Edge(n(:t1).i[:out], n(:a1).i[:in])
    Edge(n(:a1).i[:out], n(:g1).i[:mean], id=:q_mean)
    Edge(n(:t2).i[:out], n(:g1).i[:variance], id=:q_var)
    Edge(n(:g1).i[:out], n(:t3).i[:out], id=:q_out)

    return g
end

function initializeGaussianFactoringGraph()
    # [T]<--[A]<--[N]

    g = FactorGraph()
    TerminalNode(id=:t)
    GainNode(gain=[1.0], id=:gain)
    GaussianNode(m=1.0, V=0.5, id=:gauss)
    Edge(n(:gauss).i[:out], n(:gain).i[:in])
    Edge(n(:gain).i[:out], n(:t).i[:out])

    return g
end

function initializeSimpleFactoringGraph()
    # [T]<--[A]<--[T]

    g = FactorGraph()
    TerminalNode(id=:t1)
    GainNode(gain=[1.0], id=:gain)
    TerminalNode(id=:t2)
    Edge(n(:t2).i[:out], n(:gain).i[:in])
    Edge(n(:gain).i[:out], n(:t1).i[:out])

    return g
end

function initializeAdditionNode(values=[GaussianDistribution(), GaussianDistribution(), GaussianDistribution()])
    # Set up an addition node    #
    # [T]-->[+]<--[T]
    #        |
    #       [T]

    g = FactorGraph()
    AdditionNode(id=:add_node)
    for (id, value) in enumerate(values)
        Edge(TerminalNode(value).i[:out], n(:add_node).interfaces[id])
    end

    return g
end

function initializeGaussianNode(; y::ProbabilityDistribution=GaussianDistribution())
    # Initialize a Gaussian node
    #
    #    mean   precision
    #  [T]-->[N]<--[T]
    #         |
    #         v y
    #        [T]

    g = FactorGraph()
    GaussianNode(form=:precision, id=:node)
    Edge(TerminalNode(GaussianDistribution()).i[:out], n(:node).i[:mean], id=:edge1)
    Edge(TerminalNode(GammaDistribution()).i[:out], n(:node).i[:precision], id=:edge2)
    Edge(n(:node).i[:out], TerminalNode(y).i[:out], id=:edge3)

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
    t_constant = TerminalNode(DeltaDistribution(3.0))
    t_in = TerminalNode(DeltaDistribution(), id=:in)
    t_out = TerminalNode(DeltaDistribution(), id=:out)
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
        GaussianNode(form=:precision, id=:g*sec)
        EqualityNode(id=:m_eq*sec) # Equality node chain for mean
        EqualityNode(id=:gam_eq*sec) # Equality node chain for precision
        TerminalNode(y[sec], id=:y*sec) # Observed y values are stored in terminal node
        Edge(n(:g*sec).i[:out], n(:y*sec).i[:out], id=:q_y*sec)
        Edge(n(:m_eq*sec).i[3], n(:g*sec).i[:mean], id=:q_m*sec)
        Edge(n(:gam_eq*sec).i[3], n(:g*sec).i[:precision], id=:q_gam*sec)
        if sec > 1 # Connect sections
            Edge(n(:m_eq*(sec-1)).i[2], n(:m_eq*sec).i[1])
            Edge(n(:gam_eq*(sec-1)).i[2], n(:gam_eq*sec).i[1])
        end
    end
    # Attach beginning and end nodes
    TerminalNode(vague(GaussianDistribution), id=:m0) # Prior
    TerminalNode(GammaDistribution(a=1.0-tiny, b=tiny), id=:gam0) # Unifirm prior
    TerminalNode(vague(GaussianDistribution), id=:mN)
    TerminalNode(GammaDistribution(a=1.0-tiny, b=tiny), id=:gamN) # Uniform
    Edge(n(:m0).i[:out], n(:m_eq1).i[1])
    Edge(n(:gam0).i[:out], n(:gam_eq1).i[1])
    Edge(n(:m_eq*n_samples).i[2], n(:mN).i[:out])
    Edge(n(:gam_eq*n_samples).i[2], n(:gamN).i[:out])

    return g
end

function initializeMvGaussianNodeChain(y::Array{Float64, 2})
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
    n_samples = size(y, 1) # Number of observed samples

    # Build graph
    for sec=1:n_samples
        GaussianNode(form=:precision, id=:g*sec)
        EqualityNode(id=:m_eq*sec) # Equality node chain for mean
        EqualityNode(id=:gam_eq*sec) # Equality node chain for precision
        TerminalNode(MvGaussianDistribution(m=vec(y[sec, :]), V=tiny*eye(2)), id=:y*sec) # Observed y values are stored in terminal node
        Edge(n(:g*sec).i[:out], n(:y*sec).i[:out], id=:q_y*sec)
        Edge(n(:m_eq*sec).i[3], n(:g*sec).i[:mean], id=:q_m*sec)
        Edge(n(:gam_eq*sec).i[3], n(:g*sec).i[:precision], id=:q_gam*sec)
        if sec > 1 # Connect sections
            Edge(n(:m_eq*(sec-1)).i[2], n(:m_eq*sec).i[1])
            Edge(n(:gam_eq*(sec-1)).i[2], n(:gam_eq*sec).i[1])
        end
    end
    # Attach beginning and end nodes
    TerminalNode(vague(MvGaussianDistribution{2}), id=:m0) # Prior
    TerminalNode(vague(WishartDistribution{2}), id=:gam0) # Unifirm prior
    TerminalNode(vague(MvGaussianDistribution{2}), id=:mN)
    TerminalNode(vague(WishartDistribution{2}), id=:gamN) # Uniform
    Edge(n(:m0).i[:out], n(:m_eq1).i[1])
    Edge(n(:gam0).i[:out], n(:gam_eq1).i[1])
    Edge(n(:m_eq*n_samples).i[2], n(:mN).i[:out])
    Edge(n(:gam_eq*n_samples).i[2], n(:gamN).i[:out])

    return g
end

#############
# Validations
#############

function ==(x::ScheduleEntry, y::ScheduleEntry)
    if is(x, y) return true end
    ((x.outbound_interface_id == y.outbound_interface_id) && (x.node == y.node) && (x.rule == y.rule)) || (return false)
    (isdefined(x, :post_processing) == isdefined(y, :post_processing)) || (return false)
    if isdefined(x, :post_processing)
        (x.post_processing == y.post_processing) || (return false)
    end
    return true
end

function testInterfaceConnections(node1::GainNode, node2::TerminalNode)
    # Helper function for node comparison

    # Check that nodes are properly connected
    @fact typeof(node1.interfaces[1].message.payload) <: DeltaDistribution --> true
    @fact typeof(node2.interfaces[1].message.payload) <: DeltaDistribution --> true
    @fact mean(node1.interfaces[1].message.payload) --> 2.0
    @fact mean(node2.interfaces[1].message.payload) --> 1.0
    @fact mean(node1.interfaces[1].partner.message.payload) --> 1.0
    @fact mean(node2.interfaces[1].partner.message.payload) --> 2.0
    # Check that pointers are initiatized correctly
    @fact mean(node1.i[:out].message.payload) --> 3.0
    @fact mean(node2.i[:out].message.payload) --> 1.0
    @fact mean(node1.i[:in].partner.message.payload) --> 1.0
    @fact mean(node2.i[:out].partner.message.payload) --> 2.0
end

function validateOutboundMessage(node::Node, outbound_interface_index::Int, inbound_messages::Array, correct_outbound_dist::ProbabilityDistribution, update_function::Function=ForneyLab.sumProductRule!)
    # Preset an outbound distribution on which the update may operate
    if typeof(correct_outbound_dist) <: DeltaDistribution
        outbound_dist = DeltaDistribution()
    elseif typeof(correct_outbound_dist) <: MvDeltaDistribution
        outbound_dist = MvDeltaDistribution(zeros(dimensions(correct_outbound_dist)))
    else
        outbound_dist = vague(typeof(correct_outbound_dist))
    end

    # Perform the update and verify the result
    dist = update_function(node, Val{outbound_interface_index}, outbound_dist, inbound_messages...)
    @fact dist --> correct_outbound_dist

    if dist != correct_outbound_dist
        # Print full parameters if distribution is incorrect
        println("Occured:")
        for name in fieldnames(outbound_dist)
            println("$(name): $(getfield(outbound_dist, name))")
        end
        println("")
    end
end
