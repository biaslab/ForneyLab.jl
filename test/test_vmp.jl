#####################
# Integration tests
#####################

facts("VMP implementation integration tests") do
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
    ensureMWParametrization!(m_eq_nodes[end].interfaces[2].message)

    # Check the results against the outcome of Infer.NET
    accuracy = 4 # number of decimals accuracy
    m_out = m_eq_nodes[end].interfaces[2].message
    gam_out = gam_eq_nodes[end].interfaces[2].message
    @fact round(m_out.m[1], accuracy) => round(4.37638750753, accuracy)
    @fact round(m_out.V[1, 1], accuracy) => round(0.101492691239, accuracy)
    @fact round(mean(gam_out), accuracy) => round(0.984292623332, accuracy)
    @fact round(var(gam_out), accuracy) => round(0.1933796344, accuracy)
end