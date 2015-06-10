facts("Call step() for VMP algorithm") do
    data = [GaussianDistribution(m=2.0, V=tiny())]
    g = FactorGraph()
    GaussianNode(form=:precision,id=:g_node)
    TerminalNode(id=:t_out)
    TerminalNode(GaussianDistribution(), id=:t_mean)
    TerminalNode(GammaDistribution(), id=:t_var)
    Edge(n(:g_node).i[:out], n(:t_out), GaussianDistribution)
    Edge(n(:t_mean), n(:g_node).i[:mean], GaussianDistribution)
    Edge(n(:t_var), n(:g_node).i[:precision], GammaDistribution)

    attachReadBuffer(n(:t_out), data)
    mean_out = attachWriteBuffer(n(:g_node).i[:mean].edge)
    prec_out = attachWriteBuffer(n(:g_node).i[:precision].edge)

    algo = VMP.Algorithm(n_iterations=10)

    step(algo)

    ForneyLab.ensureXiWParametrization!(mean_out[end])
    @fact round(mean_out[end].W[1,1], 2) => 1.79
    @fact round(mean_out[end].xi[1], 2) => 1.57
    @fact prec_out[end].a => 1.5
    @fact round(prec_out[end].b, 2) => 1.91
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
        initializeGaussianNodeChain(data)
        n_sections = length(data)

        m_buffer = attachWriteBuffer(n(s(:m_eq,n_sections)).i[2])
        gam_buffer = attachWriteBuffer(n(s(:gam_eq,n_sections)).i[2])

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
        initializeGaussianNodeChainForSvmp(data)

        m_buffer = attachWriteBuffer(n(:m_eq).interfaces[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq).interfaces[2])
        y_buffer = attachReadBuffer(n(:y), d_data)
        Wrap(n(:mN), n(:m0))
        Wrap(n(:gamN), n(:gam0))

        # Structured factorization
        algo = VMP.Algorithm(Set{Edge}([e(:y)]), n_iterations=10)

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