#####################
# Integration tests
#####################

facts("VMP implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Integration test for the VMP implementation by trying to estimate the mean and variance of a Gaussian
        # and comparing the outcome against the outcome of the Infer.NET framework

        # Initialize chain
        (g_nodes, obs_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges) = initializeGaussianNodeChain()

        # Initialize schedules
        # Sumproduct updates
        # Update for mean equality chain
        left_update_run_m = generateSchedule(m_eq_nodes[1].interfaces[1])
        right_update_run_m = generateSchedule(m_eq_nodes[end].interfaces[2])
        downward_m = map(x -> x.interfaces[3], m_eq_nodes)
        # Update for variance equality chain
        left_update_run_gam = generateSchedule(gam_eq_nodes[1].interfaces[1])
        right_update_run_gam = generateSchedule(gam_eq_nodes[end].interfaces[2])
        downward_gam = map(x -> x.interfaces[3], gam_eq_nodes)
        # Update for samples
        upward_y = Array(Interface, 0)
        for node = g_nodes
            push!(upward_y, node.out.partner)
            push!(upward_y, node.in1)
            push!(upward_y, node.in2)
        end
        # Put it all together
        sumproduct_schedule = [left_update_run_m, right_update_run_m, downward_m, left_update_run_gam, right_update_run_gam, downward_gam, upward_y]

        # Marginal updates
        marginal_schedule = [q_m_edges, q_gam_edges]

        # Perform vmp updates
        n_its = 50
        for iter = 1:n_its
            executeSchedule(sumproduct_schedule)
            executeSchedule(marginal_schedule)
        end
        executeSchedule(sumproduct_schedule) # One last time to ensure all calculations have propagated through the equality chains

        # Save outcome
        ensureMVParametrization!(m_eq_nodes[end].interfaces[2].message.value)

        # Check the results against the outcome of Infer.NET
        accuracy = 4 # number of decimals accuracy
        m_out = m_eq_nodes[end].interfaces[2].message.value
        gam_out = gam_eq_nodes[end].interfaces[2].message.value
        @fact round(m_out.m[1], accuracy) => round(4.37638750753, accuracy)
        @fact round(m_out.V[1, 1], accuracy) => round(0.101492691239, accuracy)
        @fact round(mean(gam_out), accuracy) => round(0.984292623332, accuracy)
        @fact round(var(gam_out), accuracy) => round(0.1933796344, accuracy)
    end

    context("LinearCompositeNode linear regression parameter estimation") do
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
        n_its = 200
        for iter = 1:n_its
            executeSchedule(sumproduct_schedule)
            executeSchedule(marginal_schedule)
        end
        executeSchedule(sumproduct_schedule) # One last time to ensure all calculations have propagated through the equality chains

        # Save outcome
        ensureMVParametrization!(a_eq_nodes[end].interfaces[2].message.value)
        ensureMVParametrization!(b_eq_nodes[end].interfaces[2].message.value)
        a_out = a_eq_nodes[end].interfaces[2].message.value
        b_out = b_eq_nodes[end].interfaces[2].message.value
        s_out = s_eq_nodes[end].interfaces[2].message.value

        # Print
        println("a estimate mean $(a_out.m[1]) and variance $(a_out.V[1, 1])")
        println("b estimate mean $(b_out.m[1]) and variance $(b_out.V[1, 1])")
        println("s estimate mean $(mean(s_out)) and variance $(var(s_out))")
    end
end