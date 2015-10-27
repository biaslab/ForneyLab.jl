facts("Call step() for VMP algorithm") do
    data = [GaussianDistribution(m=2.0, V=tiny)]
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
    @fact round(mean_out[end].W[1,1], 2) --> 1.79
    @fact round(mean_out[end].xi[1], 2) --> 1.57
    @fact prec_out[end].a --> 1.5
    @fact round(prec_out[end].b, 2) --> 1.91
end

#####################
# Integration tests
#####################

facts("Naive VMP implementation integration tests") do
    context("Gaussian node mean precision batch estimation") do
        # Integration test for the VMP implementation by trying to estimate the mean and precision of a Gaussian
        # and comparing the outcome against the outcome of the Infer.NET framework

        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        initializeGaussianNodeChain(data)
        n_sections = length(data)

        m_buffer = attachWriteBuffer(n(:m_eq*n_sections).i[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq*n_sections).i[2])

        # Apply mean field factorization
        algo = VMP.Algorithm(n_iterations=50)

        # Perform vmp updates
        step(algo)

        # Check the results against the outcome of similar Infer.NET script
        m_out = m_buffer[end]
        gam_out = gam_buffer[end]

        @fact round(mean(m_out)[1], 3) --> 4.381
        @fact round(var(m_out)[1, 1], 5) --> 0.08313
        @fact round(gam_out.a, 3)  --> 6.000
        @fact round(1/gam_out.b, 4) --> 0.2005 # Scale
    end

    context("Gaussian node mean precision online estimation") do
        # Integration test for the VMP implementation by trying to estimate the mean and precision of a Gaussian
        # and comparing the outcome against the outcome of the Infer.NET framework

        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        initializeGaussianNodeChain([1.0]) # Initialize a chain with length 1
        Wrap(n(:mN), n(:m0))
        Wrap(n(:gamN), n(:gam0))

        attachReadBuffer(n(:y1), data)
        m_buffer = attachWriteBuffer(n(:m_eq1).i[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq1).i[2])

        # Apply mean field factorization
        algo = VMP.Algorithm(n_iterations=50)

        # Perform vmp updates
        run(algo)

        # Check the results against the outcome of Infer.NET
        m_out = m_buffer[end]
        gam_out = gam_buffer[end]
        @fact round(mean(m_out)[1], 3) --> 4.941
        @fact round(var(m_out)[1, 1], 3) --> 0.000 # Infer.net works with uniform gamma priors and initialization, which makes the variance collapse, (Jeffrey's prior solves this)
        @fact round(gam_out.a, 3) --> 6.000
        @fact round(1/gam_out.b, 4) --> 0.1628 # Scale
    end
end

facts("Structured VMP implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Initialize chain
        # Samples drawn from N(mean 5.0, prec 0.5): data = randn(100)*(1/sqrt(0.5))+5.0
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        d_data = [DeltaDistribution(d_k) for d_k in data]
        # d_data = [GaussianDistribution(m=d_k, W=10.0) for d_k in data]
        initializeGaussianNodeChain([0.0]) # Initialize a length 1 chain
        Wrap(n(:mN), n(:m0))
        Wrap(n(:gamN), n(:gam0))

        # n(:mN).value = vague(GaussianDistribution)
        # n(:m0).value = vague(GaussianDistribution)
        # n(:gamN).value = vague(GammaDistribution)
        # n(:gam0).value = vague(GammaDistribution)

        attachReadBuffer(n(:y1), d_data)
        m_buffer = attachWriteBuffer(n(:m_eq1).i[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq1).i[2])

        # Structured factorization
        algo = VMP.Algorithm(Set{Edge}([ForneyLab.e(:q_y1)]), n_iterations=10)

        run(algo)

        m_out = m_buffer[end]
        gam_out = gam_buffer[end]
        # Reference values from first run
        # TODO: obtain proper reference values
        @fact round(mean(m_out)[1], 3) --> 4.521
        @fact round(var(m_out)[1, 1], 3) --> 0.873 # Uniform gamma priors make the variance collapse
        @fact round(gam_out.a, 3) --> 6.000
        @fact round(1/gam_out.b, 5) --> 0.04549 # Scale
    end
end
