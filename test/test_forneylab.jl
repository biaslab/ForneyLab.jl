# This file contains the general ForneyLab tests.
# Tests for specific node and message types are
# found in test_nodes.jl and test_messages.jl 

module TestForneyLab

using FactCheck
using ForneyLab

include("test_style.jl") # Test style conventions on source files
include("integration_helpers.jl") # Helper file for integration tests, contains backgrounds and validations
include("test_helpers.jl") # Tests for ForneyLab helper methods

#####################
# Unit tests
#####################
facts("General ProbabilityDistribution unit tests") do
    for probdist_type in subtypes(ProbabilityDistribution)
        context("$(probdist_type) should have a default constructor and a == operator") do
            @fact probdist_type()==probdist_type() => true
        end
    end
end

facts("General node properties unit tests") do
    for node_type in [subtypes(Node), subtypes(CompositeNode)]
        if node_type!=CompositeNode && node_type!=MockNode
            context("$(node_type) properties should include interfaces and name") do
                @fact typeof(node_type().interfaces) => Array{Interface, 1} # Check for interface array
                @fact length(node_type().interfaces) >= 1 => true # Check length of interface array
                @fact typeof(node_type().name) => ASCIIString
            end

            context("$(node_type) constructor should assign a name") do
                my_node = node_type(;name="my_name")
                @fact my_node.name => "my_name"
            end

            context("$(node_type) constructor should assign interfaces to itself") do
                my_node = node_type()
                for interface_id in 1:length(my_node.interfaces)
                    # Check if the node interfaces couple back to the same node
                    @fact my_node.interfaces[interface_id].node => my_node
                end
            end

            context("$(node_type) should have at least 1 updateNodeMessage!() method") do
                @fact contains(string(methods(ForneyLab.updateNodeMessage!)), string("::", node_type)) => true
            end
        end
    end

    for node_type in [subtypes(CompositeNode)]
        context("$(node_type) should have property use_composite_update_rules") do
            @fact node_type().use_composite_update_rules => true || false
        end
    end
end

facts("CalculateMessage!() unit tests") do
    context("calculateMessage!() should throw an error if the specified interface does not belong to the specified node") do
        (node1, node2) = initializePairOfNodes()
        @fact_throws calculateMessage!(node1.out, node2)
    end

    context("calculateMessage!() should throw an error if one or more interfaces have no partner") do
        node = FixedGainNode()
        @fact_throws calculateMessage!(node.out)
    end
end

facts("Graph level unit tests") do
    context("FactorGraph() should initialize a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) => FactorGraph
        @fact length(fg.factorization) => 1
        @fact typeof(fg.factorization[1]) => Subgraph
        @fact typeof(fg.edge_to_subgraph) => Dict{Edge, Subgraph}
        @fact current_graph => fg # Global should be set
    end

    context("Subgraph() should initialize a subgraph and add it to the current graph") do
        sg = Subgraph()
        @fact typeof(sg) => Subgraph
        @fact typeof(sg.internal_schedule) => Schedule
        @fact typeof(sg.external_schedule) => ExternalSchedule
        graph = getCurrentGraph()
        @fact graph.factorization[2] => sg
    end

    context("getCurrentGraph() should return a current graph object") do
        FactorGraph() # Reset
        @fact getCurrentGraph() => current_graph
        my_graph = getCurrentGraph()  # Make local pointer to global variable
        @fact typeof(my_graph) => FactorGraph
    end

    context("setCurrentGraph() should set a new current graph object") do
        my_first_graph = FactorGraph() # Reset
        my_second_graph = FactorGraph()
        @fact my_first_graph == current_graph => false
        setCurrentGraph(my_first_graph)
        @fact my_first_graph == current_graph => true
    end

    context("getSubgraph(edge) should return the subgraph where edge is internal") do
        @fact true => false
    end
end

# Node and message specific tests are in separate files
include("distributions/test_gaussian.jl")
include("distributions/test_gamma.jl")
include("distributions/test_inverse_gamma.jl")
include("distributions/test_normal_gamma.jl")
include("distributions/test_students_t.jl")
include("nodes/test_addition.jl")
include("nodes/test_terminal.jl")
include("nodes/test_equality.jl")
include("nodes/test_fixed_gain.jl")
include("nodes/test_gaussian.jl")
include("nodes/composite/test_gain_addition.jl")
include("nodes/composite/test_gain_equality.jl")
include("nodes/composite/test_linear.jl")

#####################
# Integration tests
#####################

facts("Connections between nodes integration tests") do
    context("Nodes can directly be coupled through interfaces by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.interfaces[1].partner = node2.interfaces[1]
        node2.interfaces[1].partner = node1.interfaces[1]
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can directly be coupled through interfaces by using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        node1.in1.partner = node2.out
        node2.out.partner = node1.in1
        testInterfaceConnections(node1, node2)
    end

    context("Nodes can be coupled by edges by using the interfaces array") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Edge constructor should add edge and nodes to current subgraph") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        graph = getCurrentGraph()
        @fact edge in graph.factorization[1].internal_edges => true
        @fact node1 in graph.factorization[1].nodes => true
        @fact node2 in graph.factorization[1].nodes => true
    end

    context("Edge constructor should add edge to edge_to_subgraph mapping") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.interfaces[1], node1.interfaces[1]) # Edge from node 2 to node 1
        graph = getCurrentGraph()
        @fact graph.edge_to_subgraph[edge] => graph.factorization[1]
    end

    context("Nodes can be coupled by edges using the explicit interface names") do
        (node1, node2) = initializePairOfNodes()
        # Couple the interfaces that carry GeneralMessage
        edge = Edge(node2.out, node1.in1) # Edge from node 2 to node 1
        testInterfaceConnections(node1, node2)
    end

    context("Edge should throw an error when two interfaces on the same node are connected") do
        node = FixedGainNode()
        # Connect output directly to input
        @fact_throws Edge(node.interfaces[2], node.interfaces[1])
    end

    context("Edge constructor should write the expected message value types to the interfaces") do
        (node1, node2) = initializePairOfMockNodes()
        edge = Edge(node1.out, node2.out, GaussianDistribution, Float64)
        @fact edge.tail.message_payload_type => GaussianDistribution
        @fact edge.head.message_payload_type => Float64
    end

    context("Edge construction should couple interfaces to edge") do
        (node1, node2) = initializePairOfMockNodes()
        @fact node1.out.edge => nothing
        @fact node2.out.edge => nothing
        edge = Edge(node1.out, node2.out)
        @fact node1.out.edge => edge
        @fact node2.out.edge => edge
    end

    context("Edge should couple standard to composite nodes") do
        comp_node = GainEqualityCompositeNode()
        node = TerminalNode()
        edge = Edge(node.out, comp_node.in1)
        @fact comp_node.equality_node.interfaces[1].partner => node.out
        @fact comp_node.equality_node.interfaces[1].edge => edge
    end

end

facts("Graph level integration tests") do
    context("getNodes() should return an array of all nodes in the graph") do
        nodes = initializeLoopyGraph()
        found_nodes = getNodes(getCurrentGraph())
        @fact length(found_nodes) => length(nodes) # FactorGraph test
        for node in nodes
            @fact node in found_nodes => true
        end

        found_nodes = getNodes(getCurrentGraph().factorization[1]) # Subgraph test
        @fact length(found_nodes) => length(nodes)
        for node in nodes
            @fact node in found_nodes => true
        end
    end

    context("getEdges() should get all edges internal (optionally external as well) to the argument node set") do
        nodes = initializeLoopyGraph()
        @fact getEdges(Set{Node}({nodes[1], nodes[2]}), include_external=false) => Set{Edge}({nodes[1].in1.edge})
        @fact getEdges(Set{Node}({nodes[1], nodes[2]})) => Set{Edge}({nodes[1].in1.edge, nodes[4].in1.edge, nodes[4].out.edge})
    end

    context("conformSubGraph!() should complete a subgraph with nodes and external edges based in its internal edges") do
        my_graph = FactorGraph()
        # On empty subgraph
        my_subgraph = my_graph.factorization[1]
        @fact length(my_subgraph.internal_edges) => 0
        ForneyLab.conformSubgraph!(my_subgraph)
        @fact length(my_subgraph.nodes) => 0
        @fact length(my_subgraph.external_edges) => 0
        # Initialize a subgraph
        node1 = MockNode()
        node2 = MockNode(2)
        node3 = MockNode() 
        edge1 = Edge(node1.out, node2.interfaces[1])
        edge2 = Edge(node2.interfaces[2], node3.out)
        @fact length(my_subgraph.internal_edges) => 2
        ForneyLab.conformSubgraph!(my_subgraph)
        @fact length(my_subgraph.nodes) => 3
        @fact length(my_subgraph.external_edges) => 0
        # Subgraph with external edges
        new_subgraph = Subgraph(Set{Node}(), Set{Edge}({edge2}), Set{Edge}(), Array(Interface, 0), Array(Node, 0))
        @fact length(new_subgraph.internal_edges) => 1
        ForneyLab.conformSubgraph!(new_subgraph)
        @fact length(new_subgraph.nodes) => 2
        @fact length(new_subgraph.external_edges) => 1
    end

    context("addChildNodes!() should add composite node's child nodes to the node array") do
        node = initializeGainEqualityCompositeNode(eye(1), false, [Message(GaussianDistribution()), Message(GaussianDistribution()), nothing])
        @fact ForneyLab.addChildNodes!(Set{Node}({node})) => Set{Node}({node, node.equality_node, node.fixed_gain_node})
    end

    context("factorize!() should include argument edges in a new subgraph") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph()
        factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
        graph = getCurrentGraph()
        @fact graph.factorization[1].nodes => Set{Node}({driver, inhibitor, noise, add})
        @fact graph.factorization[1].internal_edges => Set{Edge}({add.out.edge, add.in1.edge, add.in2.edge})
        @fact graph.factorization[1].external_edges => Set{Edge}({inhibitor.out.edge})
        @fact graph.factorization[2].nodes => Set{Node}({driver, inhibitor})
        @fact graph.factorization[2].internal_edges => Set{Edge}({inhibitor.out.edge})
        @fact graph.factorization[2].external_edges => Set{Edge}({add.in1.edge, add.out.edge})
    end

    context("factorize!() should update the edge_to_subgraph mapping for the graph") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph()
        factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
        graph = getCurrentGraph()
        @fact graph.edge_to_subgraph[add.out.edge] => graph.factorization[1]
        @fact graph.edge_to_subgraph[add.in1.edge] => graph.factorization[1]
        @fact graph.edge_to_subgraph[add.in2.edge] => graph.factorization[1]
        @fact graph.edge_to_subgraph[inhibitor.out.edge] => graph.factorization[2]
    end

    context("factorizeMeanField!() should output a mean field factorized graph") do
        data = [1.0, 1.0, 1.0]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        graph = getCurrentGraph()
        factorizeMeanField!(graph)
        gam_set = Set{Edge}()
        for gam_eq_node in gam_eq_nodes
            for interface in gam_eq_node.interfaces
                push!(gam_set, interface.edge)
            end
        end
        m_set = Set{Edge}()
        for m_eq_node in m_eq_nodes
            for interface in m_eq_node.interfaces
                push!(m_set, interface.edge)
            end
        end
        factorized_internal_edges = [sg.internal_edges for sg in graph.factorization]
        
        @fact m_set in factorized_internal_edges => true
        @fact gam_set in factorized_internal_edges => true
        @fact Set{Edge}({q_y_edges[1]}) in factorized_internal_edges => true
        @fact Set{Edge}({q_y_edges[2]}) in factorized_internal_edges => true
        @fact Set{Edge}({q_y_edges[3]}) in factorized_internal_edges => true
    end

    context("pushRequiredInbound!() should add the proper message/marginal") do
        @fact true => false
    end
end

facts("calculateMessage!() integration tests") do
    context("calculateMessage!() should return and write back an output message") do
        (gain, terminal) = initializePairOfNodes(A=[2.0], msg_gain_1=nothing, msg_gain_2=nothing, msg_terminal=Message(3.0))
        Edge(terminal.out, gain.in1, Float64, Array{Float64, 2})
        gain.out.message_payload_type = Array{Float64, 2} # Expect a matrix
        @fact gain.out.message => nothing
        # Request message on node for which the input is unknown
        msg = calculateMessage!(gain.out)
        @fact msg => gain.out.message # Returned message should be identical to message stored on interface.
        @fact typeof(gain.out.message.payload) => Array{Float64, 2}
        @fact gain.out.message.payload => reshape([6.0], 1, 1)
    end

    context("calculateMessage!() should recursively calculate required inbound message") do
        # Define three nodes in series
        (node1, node2, node3) = initializeChainOfNodes()
        @fact node3.out.message => nothing
        # Request message on node for which the input is unknown
        node3.out.message_payload_type = Array{Float64, 2} # Expect a matrix
        calculateMessage!(node3.out)
        @fact typeof(node3.out.message.payload) => Array{Float64, 2}
        @fact node3.out.message.payload => reshape([12.0], 1, 1)
    end

    context("calculateMessage!() should throw an error when there is an unbroken loop") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
        @fact_throws calculateMessage!(driver.out)
    end
end

facts("generateSchedule() and executeSchedule() integration tests") do

    # Begin of graph context
    context("Graph context") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)

        context("generateSchedule() should throw an error when there is an unbroken loop") do
            @fact_throws generateSchedule(driver.out)
        end

        # Initial message
        setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
        setMessage!(add.out, Message(GaussianDistribution()))

        context("generateSchedule() should auto-generate a feasible schedule") do
            # Generate schedule automatically
            schedule = generateSchedule(add.in2) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact inhibitor.out in schedule => true
            @fact driver.out    in schedule => true
            @fact inhibitor.in1 in schedule => true
            @fact driver.in1    in schedule => true
            @fact add.in2       in schedule => true
            @fact add.in1       in schedule => false
            @fact add.out       in schedule => false
            @fact noise.out     in schedule => false
            # Validate correct relative order in schedule
            @fact findfirst(schedule, inhibitor.out)    < findfirst(schedule, driver.out)   => true
            @fact findfirst(schedule, driver.out)       < findfirst(schedule, add.in2)      => true
            @fact findfirst(schedule, driver.in1)       < findfirst(schedule, inhibitor.in1)=> true
            @fact findfirst(schedule, inhibitor.in1)    < findfirst(schedule, add.in2)      => true
        end

        context("generateSchedule() should correctly complete a partial schedule") do
            # Generate a schedule that first passes clockwise through the cycle and then counterclockwise
            schedule = generateSchedule([driver.out, add.in2]) # Message towards noise factor
            # All (but just) required calculations should be in the schedule
            @fact schedule[1] => inhibitor.out
            @fact schedule[2] => driver.out
            @fact schedule[3] => driver.in1
            @fact schedule[4] => inhibitor.in1
            @fact schedule[5] => add.in2
        end

        context("executeSchedule() should correctly execute a schedule and return the result of the last step") do
            schedule = generateSchedule(add.in2)
            dist = ensureMVParametrization!(executeSchedule(schedule).payload)
            @fact dist => add.in2.message.payload
            @fact isApproxEqual(dist.m, [2.0]) => true
            @fact isApproxEqual(dist.V, reshape([1.5], 1, 1)) => true
        end

        context("generateSchedule!() should generate an internal and external schedule when called on a subgraph") do
            (driver, inhibitor, noise, add) = initializeLoopyGraph()
            factorize!(Set{Edge}({inhibitor.out.edge})) # Put this edge in a different subgraph
            graph = getCurrentGraph()
            for subgraph in graph.factorization
                #graphViz(subgraph)
                generateSchedule!(subgraph)
                #println(subgraph.internal_schedule)
                #println(subgraph.external_schedule)
            end
            # Order differs between test runs, validity of order is not checked
            @fact Set{Interface}(graph.factorization[1].internal_schedule) => Set{Interface}([noise.out, inhibitor.in1, add.in1, driver.out, add.out])
            @fact Set{Interface}(graph.factorization[2].internal_schedule) => Set{Interface}([driver.in1, inhibitor.out])
            @fact Set{Node}(graph.factorization[1].external_schedule) => Set{Node}([driver, inhibitor])
            @fact Set{Node}(graph.factorization[2].external_schedule) => Set{Node}([inhibitor, driver])
        end
    end # End of graph context

    context("executeSchedule() should work as expeced in loopy graphs") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
        setMessage!(driver.out, Message(GaussianDistribution()))
        schedule = generateSchedule(driver.out)
        for count = 1:100
            executeSchedule(schedule)
        end
        @fact typeof(driver.out.message) => Message{GaussianDistribution}
        @fact ensureMVParametrization!(driver.out.message.payload).m => [100.0] # For stop conditions at 100 cycles deep
    end

    context("executeSchedule() should be called repeatedly until convergence") do
        (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[1.1], B=[0.1], noise_m=0.0, noise_V=0.1)
        # Now set a breaker message and check that it works
        breaker_message = Message(GaussianDistribution(m=10.0, V=100.0))
        setMessage!(driver.out, breaker_message)
        prev_dist = deepcopy(breaker_message.payload)
        converged = false
        schedule = generateSchedule(driver.out)
        while !converged
            dist = ensureMVParametrization!(executeSchedule(schedule).payload)
            converged = isApproxEqual(prev_dist.m, dist.m)
            prev_dist = deepcopy(dist)
        end
        @fact isApproxEqual(driver.out.message.payload.m, [0.0]) => true
    end
end

facts("clearMessage!(), clearMessages!() integration tests") do
    (driver, inhibitor, noise, add) = initializeLoopyGraph(A=[2.0], B=[0.5], noise_m=1.0, noise_V=0.1)
    setMessage!(add.in1, Message(GaussianDistribution(m=2.0, V=0.5)))
    setMessage!(add.out, Message(GaussianDistribution()))
    schedule = generateSchedule(add.in2)
    executeSchedule(schedule)
    clearMessage!(add.in2)
    @fact add.in2.message => nothing
    clearMessages!(add)
    @fact add.in1.message => nothing
    @fact add.out.message => nothing
end

# Vmp test
include("test_vmp.jl")

try
    # Try to load user-defined extensions tests
    include("$(Main.FORNEYLAB_EXTENSION_DIR)/test/test_forneylab_extensions.jl")
end

end # module TestForneyLab