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
            push!(upward_y, node.mean)
            push!(upward_y, node.precision)
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
        ensureMVParametrization!(m_eq_nodes[end].interfaces[2].message.payload)

        # Check the results against the outcome of Infer.NET
        accuracy = 3 # number of decimals accuracy
        m_out = m_eq_nodes[end].interfaces[2].message.payload
        gam_out = gam_eq_nodes[end].interfaces[2].message.payload
        @fact round(m_out.m[1], accuracy) => round(4.37638750753, accuracy)
        @fact round(m_out.V[1, 1], accuracy+1) => round(0.101492691239, accuracy+1)
        @fact round(mean(gam_out), accuracy+1) => round(0.984292623332, accuracy+1)
        @fact round(var(gam_out), accuracy+1) => round(0.1933796344, accuracy+1)
    end

    context("LinearCompositeNode linear regression parameter estimation") do
        (lin_nodes, a_eq_nodes, b_eq_nodes, gam_eq_nodes, a_eq_edges, b_eq_edges, gam_eq_edges, x_edges, y_edges) = initializeLinearCompositeNodeChain()

        # Equality chain schedules
        left_update_run_a = generateSchedule(a_eq_nodes[1].interfaces[1]) # a
        right_update_run_a = generateSchedule(a_eq_nodes[end].interfaces[2])
        downward_a = map(x -> x.interfaces[3], a_eq_nodes)
        left_update_run_b = generateSchedule(b_eq_nodes[1].interfaces[1]) # b
        right_update_run_b = generateSchedule(b_eq_nodes[end].interfaces[2])
        downward_b = map(x -> x.interfaces[3], b_eq_nodes)
        left_update_run_gam = generateSchedule(gam_eq_nodes[1].interfaces[1]) # s
        right_update_run_gam = generateSchedule(gam_eq_nodes[end].interfaces[2])
        downward_gam = map(x -> x.interfaces[3], gam_eq_nodes)
        # Update for samples
        node_update = Array(Interface, 0)
        for node = lin_nodes
            push!(node_update, node.out.partner)
            push!(node_update, node.in1.partner)
            push!(node_update, node.slope)
            push!(node_update, node.offset)
            push!(node_update, node.noise)
        end
        # Put it all together
        sumproduct_schedule = [left_update_run_a, right_update_run_a, downward_a, left_update_run_b, right_update_run_b, downward_b, left_update_run_gam, right_update_run_gam, downward_gam, node_update]

        # Marginal updates
        marginal_schedule = [a_eq_edges, b_eq_edges, gam_eq_edges, x_edges, y_edges]

        # Perform vmp updates
        n_its = 100
        for iter = 1:n_its
            executeSchedule(sumproduct_schedule)
            executeSchedule(marginal_schedule)
        end
        executeSchedule(sumproduct_schedule) # One last time to ensure all calculations have propagated through the equality chains

        # Check the results against the outcome of Infer.NET
        ensureMVParametrization!(a_eq_nodes[end].interfaces[2].message.payload)
        ensureMVParametrization!(b_eq_nodes[end].interfaces[2].message.payload)
        a_out = a_eq_nodes[end].interfaces[2].message.payload
        b_out = b_eq_nodes[end].interfaces[2].message.payload
        gam_out = gam_eq_nodes[end].interfaces[2].message.payload

        accuracy = 1

        @fact round(mean(a_out), accuracy)[1] => round(2.92642601384, accuracy)
        @fact round(var(a_out), accuracy+4)[1, 1] => round(0.000493670181134, accuracy+4)
        @fact round(mean(b_out), accuracy)[1] => round(5.85558752435, accuracy)
        @fact round(var(b_out), accuracy+2)[1, 1] => round(0.0609314195382, accuracy+2)
        @fact round(mean(gam_out), accuracy+1) => round(0.820094703716, accuracy+1)
        @fact round(var(gam_out), accuracy+2) => round(0.0671883439624, accuracy+2)
    end
end