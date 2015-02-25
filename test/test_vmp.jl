# #####################
# # Integration tests
# #####################

facts("Naive VMP implementation integration tests") do
    context("Gaussian node mean precision estimation") do
        # Integration test for the VMP implementation by trying to estimate the mean and precision of a Gaussian
        # and comparing the outcome against the outcome of the Infer.NET framework

        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        (g_nodes, y_nodes, m_eq_nodes, gam_eq_nodes, q_m_edges, q_gam_edges, q_y_edges) = initializeGaussianNodeChain(data)
        n_sections = length(data)

        # Apply mean field factorization
        factorize!()

        graph = currentGraph()
        for subgraph in graph.factorization
            generateSchedule!(subgraph) # Generate internal and external schedule automatically
        end

        setVagueMarginals!()

        # Perform vmp updates
        n_its = 50

        subgraph_m = subgraph(g_nodes[1].mean.edge)
        subgraph_gam = subgraph(g_nodes[1].precision.edge)
        subgraphs_y = [subgraph(e) for e in q_y_edges]
        for iter = 1:n_its
            for sg_y in subgraphs_y
                execute(sg_y)
            end
            # q(m) update
            execute(subgraph_m)
            # q(gam) update
            execute(subgraph_gam)
        end
        # One last time to ensure all calculations have propagated through the equality chains
        execute(ForneyLab.convert(Schedule, [m_eq_nodes[end].interfaces[2]]))
        execute(ForneyLab.convert(Schedule, [gam_eq_nodes[end].interfaces[2]]))

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
end

facts("Structured VMP implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Initialize chain
        # 100 samples drawn from N(mean 5.0, prec 0.5): data = randn(100)*(1/sqrt(0.5))+5.0
        true_mean = 5.0
        true_prec = 0.5
        data = [2.9133739230396,6.776946270056516,5.699964237997796,3.0119520131244513,4.993687477925646,2.881070146797452,7.629830860404964,5.954041063509354,7.5532738533379575,3.6611337404678705,4.511142300904889,5.550456117201564,8.034701344789596,6.6650044853725,7.101840837092003,4.6008371450675325,6.346078557872094,4.332171347481462,3.3241857462173927,3.5173617282549094,7.041611210311107,5.650982184976964,5.409847595551146,5.983741217058288,6.955368609201553,3.7551413767655166,6.777625369831803,3.3221445669751453,3.958075930250893,3.7782063708759006,6.248367587394229,7.706396857185407,5.925106466541943,7.275126285133447,2.7894263038295897,6.301796475025749,4.944969659867805,2.406699646675323,3.297436847112211,6.128679686897025,5.607333293256828,3.5895918562291813,6.811148920896203,4.859402517744455,2.5918075356111885,5.76730643469031,5.78370631320422,5.834672226384856,4.883023830265342,6.6521709869249745,2.155223456720972,5.9361238868926,4.732878170437103,5.888299163098336,4.90977472267389,4.306006320861194,6.179382449782395,2.412314907046015,5.164360962519157,3.817047666470755,5.163951789662665,6.449495630973551,4.304708322895846,3.790402855120055,5.42744802571948,3.7155725574003826,5.718747174625624,6.246789859516861,5.100705318199726,4.46915729683993,5.316181918934593,2.4373233015460936,5.3718738155266355,6.894455289139208,6.17653704887158,5.730586963905106,3.911294862495409,5.638772864526039,3.5131576213804454,4.994037568420812,3.835497990047835,4.19953408648465,3.5664907542111877,5.067659198961311,4.131824295081574,6.583043832829379,6.440797611033075,6.615011690694574,4.789467095878956,3.0398417343094426,6.140798845980758,3.8314793542388528,6.003274190743673,4.959480584969705,4.908735288499479,6.892347993289387,6.144780127407775,3.136896288776731,4.185693744866867,3.164021612264378]
        n_samples = length(data)
        n_its = 10
        (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge) = initializeGaussianNodeChainForSvmp(data)

        # Structured factorization
        factorize!(Set{Edge}([y_edge]))

        graph = currentGraph()
        for subgraph in graph.factorization
            generateSchedule!(subgraph) # Generate internal and external schedule automatically
        end

        setVagueMarginals!()

        y_subgraph = subgraph(g_node.out.edge)
        m_gam_subgraph = subgraph(g_node.mean.edge)
        for sample = 1:n_samples
            # Reset
            y_node.value = GaussianDistribution(m = data[sample], W=10.0) # Small variance on sample
            # Reset vague ('one') messages
            setVagueMarginals!()

            # Do the VMP iterations
            for it = 1:n_its
                execute(m_gam_subgraph)
                execute(y_subgraph)
            end
            # Propagate through chain
            execute(ForneyLab.convert(Schedule, [g_node.mean, m_eq_node.interfaces[2]]))
            execute(ForneyLab.convert(Schedule, [g_node.precision, gam_eq_node.interfaces[2]]))

            # Switch posterior to prior for next sample
            m_0_node.value = deepcopy(m_eq_node.interfaces[2].message.payload)
            gam_0_node.value = deepcopy(gam_eq_node.interfaces[2].message.payload)
        end

        m_m = mean(m_eq_node.interfaces[2].message.payload)[1]
        m_sigma = sqrt(var(m_eq_node.interfaces[2].message.payload)[1])
        gam_m = mean(gam_eq_node.interfaces[2].message.payload)
        gam_sigma = sqrt(var(gam_eq_node.interfaces[2].message.payload))
    
        # Check for small enough sigma
        @fact m_sigma < 0.2 => true
        @fact gam_sigma < 0.2 => true
        # Check for correctness of estimation withing 1 sigma range
        @fact m_m-m_sigma < true_mean < m_m+m_sigma => true
        @fact gam_m-gam_sigma < true_prec < gam_m+gam_sigma => true
    end
end