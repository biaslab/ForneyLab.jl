#####################
# Unit tests
#####################

facts("LinearCompositeNode unit tests") do
    context("LinearCompositeNode() should initialize a LinearCompositeNode with 5 interfaces") do
        node = LinearCompositeNode()
        @fact typeof(node) => LinearCompositeNode
        @fact length(node.interfaces) => 5
        @fact node.in1 => node.interfaces[1]
        @fact node.a_in => node.interfaces[2]
        @fact node.b_in => node.interfaces[3]
        @fact node.s_in => node.interfaces[4]
        @fact node.out => node.interfaces[5]
        @fact node.variational => true # default variational to true
        @fact node.use_composite_update_rules => true # default use_composite_update_rules to true
    end
end

#####################
# Integration tests
#####################

facts("LinearCompositeNode integration tests") do
    context("LinearCompositeNode should propagate a backward variational message to in1") do
        #lin_node = initializeLinearCompositeNode()
    end

    context("LinearCompositeNode should propagate a backward variational message to a_in") do
    end

    context("LinearCompositeNode should propagate a backward variational message to b_in") do
    end

    context("LinearCompositeNode should propagate a backward variational message to s_in") do
    end

    context("LinearCompositeNode should propagate a forward variational message to out") do
    end

    context("LinearCompositeNode should perform regression") do
        (lin_nodes, a_eq_nodes, b_eq_nodes, s_eq_nodes, a_eq_edges, b_eq_edges, s_eq_edges, x_edges, y_edges) = initializeLinearCompositeNodeChain()

        # Equality chain schedules
        left_update_run_a = generateSchedule(a_eq_nodes[1].interfaces[1]) # a
        right_update_run_a = generateSchedule(a_eq_nodes[end].interfaces[2])
        downward_a = map(x -> x.interfaces[3], a_eq_nodes)
        left_update_run_b = generateSchedule(b_eq_nodes[1].interfaces[1]) # b
        right_update_run_b = generateSchedule(b_eq_nodes[end].interfaces[2])
        downward_b = map(x -> x.interfaces[3], b_eq_nodes)
        left_update_run_s = generateSchedule(s_eq_nodes[1].interfaces[1]) # s
        right_update_run_s = generateSchedule(s_eq_nodes[end].interfaces[2])
        downward_s = map(x -> x.interfaces[3], s_eq_nodes)
        # Update for samples
        node_update = Array(Interface, 0)
        for node = lin_nodes
            push!(node_update, node.out.partner)
            push!(node_update, node.in1.partner)
            push!(node_update, node.a_in)
            push!(node_update, node.b_in)
            push!(node_update, node.s_in)
        end
        # Put it all together
        sumproduct_schedule = [left_update_run_a, right_update_run_a, downward_a, left_update_run_b, right_update_run_b, downward_b, left_update_run_s, right_update_run_s, downward_s, node_update]

        # Marginal updates
        marginal_schedule = [a_eq_edges, b_eq_edges, s_eq_edges, x_edges, y_edges]

        # Perform vmp updates
        n_its = 50
        for iter = 1:n_its
            executeSchedule(sumproduct_schedule)
            executeSchedule(marginal_schedule)
        end
        executeSchedule(sumproduct_schedule) # One last time to ensure all calculations have propagated through the equality chains

        # Save outcome
        ensureMVParametrization!(a_eq_nodes[end].interfaces[2].message)
        ensureMVParametrization!(b_eq_nodes[end].interfaces[2].message)
        a_out = a_eq_nodes[end].interfaces[2].message
        b_out = b_eq_nodes[end].interfaces[2].message
        s_out = s_eq_nodes[end].interfaces[2].message

        # Print
        println("a estimate mean $(a_out.m[1]) and variance $(a_out.V[1, 1])")
        println("b estimate mean $(b_out.m[1]) and variance $(b_out.V[1, 1])")
        println("s estimate mean $(mean(s_out)) and variance $(var(s_out))")
    end
end