facts("step() function for VMP") do
    context("step should accept and execute a number of iterations for VMP") do
        data = [GaussianDistribution(m=2.0, V=tiny())]
        g = FactorGraph()
        g_node = GaussianNode(form="precision")
        t_out = TerminalNode(name="t_out")
        t_mean = TerminalNode(GaussianDistribution(), name="t_mean")
        t_var = TerminalNode(GammaDistribution(), name="t_var")
        Edge(g_node.out, t_out, GaussianDistribution)
        Edge(t_mean, g_node.mean, GaussianDistribution)
        Edge(t_var, g_node.precision, GammaDistribution)

        setReadBuffer(t_out, data)
        mean_out = setWriteBuffer(g_node.mean.edge)
        prec_out = setWriteBuffer(g_node.precision.edge)

        algo = VMP.Algorithm(n_iterations=10)

        step(algo)

        @fact round(mean_out[end].W[1,1], 2) => 1.79
        @fact round(mean_out[end].xi[1], 2) => 1.57
        @fact prec_out[end].a => 1.5
        @fact round(prec_out[end].b, 2) => 1.91
    end
end

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

        m_buffer = setWriteBuffer(m_eq_nodes[end].interfaces[2])
        gam_buffer = setWriteBuffer(gam_eq_nodes[end].interfaces[2])

        # Apply mean field factorization
        algo = VMP.Algorithm(n_iterations=20)

        # Perform vmp updates
        step(algo)

        # Check the results against the outcome of Infer.NET
        accuracy = 3 # number of decimals accuracy
        m_out = m_buffer[end]
        gam_out = gam_buffer[end]
        @fact round(mean(m_out)[1], accuracy) => round(4.37638750753, accuracy)
        @fact round(var(m_out)[1, 1], accuracy+1) => round(0.101492691239, accuracy+1)
        @fact round(mean(gam_out), accuracy+1) => round(0.984292623332, accuracy+1)
        @fact round(var(gam_out), accuracy+1) => round(0.1933796344, accuracy+1)
    end

    context("Sigmoid node variational estimation") do
        # Data
        y = [1 0 0 1 0 1 0 1 0 0 1 0 0 1 1]
        data = BetaDistribution[]
        for y_k in y
            if y_k == 1
                push!(data, BetaDistribution(a=1.0, b=0.0))
            else
                push!(data, BetaDistribution(a=0.0, b=1.0))
            end
        end

        graph = initializeSigmoidSlice()

        setTimeWrap(node("theta_k"), node("theta_k_min"))
        setReadBuffer(node("y"), data)
        state_buff = setWriteBuffer(node("eq").interfaces[2])

        # Initialize algo for a mean field factorization
        algo = VMP.Algorithm(n_iterations = 5)

        # Set prior over theta
        node("theta_k_min").value = GaussianDistribution(m=0.0, V=100.0)

        # Perform inference
        #
        # We need to make sure that the downward variational message is computed first,
        # otherwise the sigmoid node can't compute the upward variational message.
        run(algo)

        @fact state_buff[end].m[1] => -0.03730307502629462     
        @fact state_buff[end].W[1,1] => 80.3323823509323      
    end
end

facts("Structured VMP implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Initialize chain
        # 100 samples drawn from N(mean 5.0, prec 0.5): data = randn(100)*(1/sqrt(0.5))+5.0
        true_mean = 5.0
        true_prec = 0.5
        data = [2.9133739230396,6.776946270056516,5.699964237997796,3.0119520131244513,4.993687477925646,2.881070146797452,7.629830860404964,5.954041063509354,7.5532738533379575,3.6611337404678705,4.511142300904889,5.550456117201564,8.034701344789596,6.6650044853725,7.101840837092003,4.6008371450675325,6.346078557872094,4.332171347481462,3.3241857462173927,3.5173617282549094,7.041611210311107,5.650982184976964,5.409847595551146,5.983741217058288,6.955368609201553,3.7551413767655166,6.777625369831803,3.3221445669751453,3.958075930250893,3.7782063708759006,6.248367587394229,7.706396857185407,5.925106466541943,7.275126285133447,2.7894263038295897,6.301796475025749,4.944969659867805,2.406699646675323,3.297436847112211,6.128679686897025,5.607333293256828,3.5895918562291813,6.811148920896203,4.859402517744455,2.5918075356111885,5.76730643469031,5.78370631320422,5.834672226384856,4.883023830265342,6.6521709869249745,2.155223456720972,5.9361238868926,4.732878170437103,5.888299163098336,4.90977472267389,4.306006320861194,6.179382449782395,2.412314907046015,5.164360962519157,3.817047666470755,5.163951789662665,6.449495630973551,4.304708322895846,3.790402855120055,5.42744802571948,3.7155725574003826,5.718747174625624,6.246789859516861,5.100705318199726,4.46915729683993,5.316181918934593,2.4373233015460936,5.3718738155266355,6.894455289139208,6.17653704887158,5.730586963905106,3.911294862495409,5.638772864526039,3.5131576213804454,4.994037568420812,3.835497990047835,4.19953408648465,3.5664907542111877,5.067659198961311,4.131824295081574,6.583043832829379,6.440797611033075,6.615011690694574,4.789467095878956,3.0398417343094426,6.140798845980758,3.8314793542388528,6.003274190743673,4.959480584969705,4.908735288499479,6.892347993289387,6.144780127407775,3.136896288776731,4.185693744866867,3.164021612264378]
        d_data = [GaussianDistribution(m=d_k, W=10.0) for d_k in data]
        n_samples = length(data)
        (g_node, y_node, m_0_node, gam_0_node, m_N_node, gam_N_node, m_eq_node, gam_eq_node, m_edge, gam_edge, y_edge) = initializeGaussianNodeChainForSvmp(data)

        m_buffer = setWriteBuffer(m_eq_node.interfaces[2])
        gam_buffer = setWriteBuffer(gam_eq_node.interfaces[2])
        y_buffer = setReadBuffer(y_node, d_data)
        setTimeWrap(m_N_node, m_0_node)
        setTimeWrap(gam_N_node, gam_0_node)

        # Structured factorization
        algo = VMP.Algorithm([Set{Edge}([y_edge])], n_iterations=10)

        run(algo)

        m_m = mean(m_buffer[end])[1]
        m_sigma = sqrt(var(m_buffer[end])[1,1])
        gam_m = mean(gam_buffer[end])
        gam_sigma = sqrt(var(gam_buffer[end]))

        # Check for small enough sigma
        @fact m_sigma < 0.2 => true
        @fact gam_sigma < 0.2 => true
        # Check for correctness of estimation withing 1 sigma range
        @fact m_m-m_sigma < true_mean < m_m+m_sigma => true
        @fact gam_m-gam_sigma < true_prec < gam_m+gam_sigma => true
    end
end