#####################
# Integration tests
#####################

facts("Naive VMP implementation integration tests") do
    context("Gaussian node mean precision estimation") do
        # Integration test for the VMP implementation by trying to estimate the mean and precision of a Gaussian
        # and comparing the outcome against the outcome of the Infer.NET framework

        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)

        # Initialize schedules
        # Sumproduct and marginal schedule for q(y)
        sumproduct_y_schedule = Array(Interface, 0)
        for section = 1:n_sections 
            push!(sumproduct_y_schedule, y_nodes[section].out) # Backward y message
            push!(sumproduct_y_schedule, g_nodes[section].out) # Forward y message
        end
        marginal_y_schedule = q_y_edges

        # Sumproduct and marginal schedule for q(m)
        sumproduct_m_schedule = Array(Interface, 0)
        for section = 1:n_sections; push!(sumproduct_m_schedule, g_nodes[section].mean); end # Backward
        push!(sumproduct_m_schedule, m_eq_nodes[1].interfaces[1].partner) # Prior message        
        for section = 1:n_sections; push!(sumproduct_m_schedule, m_eq_nodes[section].interfaces[2]); end # Forward run
        push!(sumproduct_m_schedule, m_eq_nodes[end].interfaces[2].partner) # Terminal message        
        for section = n_sections:-1:1; push!(sumproduct_m_schedule, m_eq_nodes[section].interfaces[1]); end # Backward run
        for section = 1:n_sections; push!(sumproduct_m_schedule, m_eq_nodes[section].interfaces[3]); end # Forward (downward run)
        marginal_m_schedule = q_m_edges

        # Sumproduct and marginal schedule for q(gam)
        sumproduct_gam_schedule = Array(Interface, 0)
        for section = 1:n_sections; push!(sumproduct_gam_schedule, g_nodes[section].precision); end # Backward
        push!(sumproduct_gam_schedule, gam_eq_nodes[1].interfaces[1].partner) # Prior message        
        for section = 1:n_sections; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[2]); end # Forward run
        push!(sumproduct_gam_schedule, gam_eq_nodes[end].interfaces[2].partner) # Terminal message        
        for section = n_sections:-1:1; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[1]); end # Backward run
        for section = 1:n_sections; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[3]); end # Forward (downward run)
        marginal_gam_schedule = q_gam_edges

        # Perform vmp updates
        n_its = 50
        # q(y) update
        # We need to execute the q(y) updates only once, because sample values do not change.
        executeSchedule(sumproduct_y_schedule)
        executeSchedule(marginal_y_schedule)
        for iter = 1:n_its
            # q(m) update
            executeSchedule(sumproduct_m_schedule)
            executeSchedule(marginal_m_schedule)
            # q(gam) update
            executeSchedule(sumproduct_gam_schedule)
            executeSchedule(marginal_gam_schedule)
        end
        # One last time to ensure all calculations have propagated through the equality chains
        executeSchedule(sumproduct_m_schedule)
        executeSchedule(sumproduct_gam_schedule)

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
        #true_gam = 0.5
        #true_a = 3.0
        #true_b = 5.0
        x = [0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,17.0,18.0,19.0]
        y = [6.1811923622357625,7.917496269084679,11.286102016681964,14.94255088702814,16.82264686442818,19.889355802073506,23.718253510300688,28.18105443765643,27.72075206943362,32.15921446069328,34.97262678800721,38.86444301740928,40.79138365100076,45.84963364094473,47.818481172238165,51.51027620022872,52.623019301773,53.91583839744505,58.14426361122961,59.895517438500164]
        (lin_nodes, a_eq_nodes, b_eq_nodes, gam_eq_nodes, a_eq_edges, b_eq_edges, gam_eq_edges, x_edges, y_edges) = initializeLinearCompositeNodeChain(x, y)
        n_sections = length(y)

        # Sumproduct and marginal schedule for q(a)
        sumproduct_a_schedule = Array(Interface, 0)
        for section = 1:n_sections; push!(sumproduct_a_schedule, lin_nodes[section].slope); end # Backward
        push!(sumproduct_a_schedule, a_eq_nodes[1].interfaces[1].partner) # Prior message        
        for section = 1:n_sections; push!(sumproduct_a_schedule, a_eq_nodes[section].interfaces[2]); end # Forward run
        push!(sumproduct_a_schedule, a_eq_nodes[end].interfaces[2].partner) # Terminal message        
        for section = n_sections:-1:1; push!(sumproduct_a_schedule, a_eq_nodes[section].interfaces[1]); end # Backward run
        for section = 1:n_sections; push!(sumproduct_a_schedule, a_eq_nodes[section].interfaces[3]); end # Forward (downward run)
        marginal_a_schedule = a_eq_edges

        # Sumproduct and marginal schedule for q(b)
        sumproduct_b_schedule = Array(Interface, 0)
        for section = 1:n_sections; push!(sumproduct_b_schedule, lin_nodes[section].offset); end # Backward
        push!(sumproduct_b_schedule, b_eq_nodes[1].interfaces[1].partner) # Prior message        
        for section = 1:n_sections; push!(sumproduct_b_schedule, b_eq_nodes[section].interfaces[2]); end # Forward run
        push!(sumproduct_b_schedule, b_eq_nodes[end].interfaces[2].partner) # Terminal message        
        for section = n_sections:-1:1; push!(sumproduct_b_schedule, b_eq_nodes[section].interfaces[1]); end # Backward run
        for section = 1:n_sections; push!(sumproduct_b_schedule, b_eq_nodes[section].interfaces[3]); end # Forward (downward run)
        marginal_b_schedule = b_eq_edges

        # Sumproduct and marginal schedule for q(gam)
        sumproduct_gam_schedule = Array(Interface, 0)
        for section = 1:n_sections; push!(sumproduct_gam_schedule, lin_nodes[section].noise); end # Backward
        push!(sumproduct_gam_schedule, gam_eq_nodes[1].interfaces[1].partner) # Prior message        
        for section = 1:n_sections; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[2]); end # Forward run
        push!(sumproduct_gam_schedule, gam_eq_nodes[end].interfaces[2].partner) # Terminal message        
        for section = n_sections:-1:1; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[1]); end # Backward run
        for section = 1:n_sections; push!(sumproduct_gam_schedule, gam_eq_nodes[section].interfaces[3]); end # Forward (downward run)
        marginal_gam_schedule = gam_eq_edges

        # x and y marginals are already set upon initialization. These need to be set only once because samples remain the same.

        # Perform vmp updates
        n_its = 100
        for iter = 1:n_its
            executeSchedule(sumproduct_a_schedule)
            executeSchedule(marginal_a_schedule)
            executeSchedule(sumproduct_b_schedule)
            executeSchedule(marginal_b_schedule)
            executeSchedule(sumproduct_gam_schedule)
            executeSchedule(marginal_gam_schedule)
        end
        executeSchedule(sumproduct_a_schedule) # One last time to ensure all calculations have propagated through the equality chains
        executeSchedule(sumproduct_b_schedule)
        executeSchedule(sumproduct_gam_schedule)

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

facts("Structured VMP implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        n_samples = length(data)
        (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge) = initializeGaussianNodeChainForSvmp(data)
        n_its = 50 # Number of VMP iterations

        # Update schedules for q(y) subgraph
        y_sumproduct_schedule = [g_node.out, y_node.out]
        q_y_marginal_schedule = [y_edge]

        # Update schedules for q(m, gam) subgraph
        m_gam_sumproduct_schedule = [m_0_node.out, m_N_node.out, m_eq_node.interfaces[3], gam_0_node.out, gam_N_node.out, gam_eq_node.interfaces[3], g_node.mean, m_eq_node.interfaces[2], g_node.precision, gam_eq_node.interfaces[2]]
        q_m_gam_marginal_schedule = [g_node] # Not used in practice, but here for completeness. We do not really need to calculate the downward y message, because its mean will collapse to the sample value.

        for sample = 1:n_samples
            # Reset
            # Preset uninformative ('one') messages
            setMarginal!(g_node, uninformative(NormalGammaDistribution))
            setMarginal!(y_edge, uninformative(Float64)) # 1
            #y_node.value = data[sample] # Set sample

            # Do the VMP 'iterations'
            # TODO: does this reduce to the standard sumproduct rule for a fixed sample?
            #executeSchedule(y_sumproduct_schedule) # y updates only once, because sample value does not change
            #executeSchedule(q_y_marginal_schedule)
            y_edge.marginal=data[sample] # The quick and dirty way
            executeSchedule(m_gam_sumproduct_schedule) # m_gam updates only once, because q(y) does not change
            #executeSchedule(q_m_gam_marginal_schedule) # for show

            # Switch posterior to prior for next sample
            m_0_node.value = deepcopy(m_eq_node.interfaces[2].message.payload)
            gam_0_node.value = deepcopy(gam_eq_node.interfaces[2].message.payload)
        end

        println("sample mean = $(mean(data))")
        println("sample variance = $(var(data))")
        println("m estimate = $(mean(m_eq_node.interfaces[2].message.payload)) var: $(var(m_eq_node.interfaces[2].message.payload))")
        println("var estimate = $(mean(gam_eq_node.interfaces[2].message.payload)) var: $(var(gam_eq_node.interfaces[2].message.payload))")
        @fact true => false
    end
end