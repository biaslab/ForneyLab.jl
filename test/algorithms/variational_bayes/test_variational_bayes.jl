#####################
# Integration tests
#####################

facts("VariationalBayes should initialize a pre, iterative and post schedule") do
    # Initialize chain
    data = [0.0, 0.0]
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    m_buffer = attachWriteBuffer(n(:m_eq*n_sections).i[2]) # Propagate to the end
    gam_buffer = attachWriteBuffer(n(:gam_eq*n_sections).i[2])

    # Apply mean field factorization
    rf = RecognitionFactorization()
    factorizeMeanField()
    for sec in 1:n_sections
        initialize(eg(:q_y*sec), vague(Gaussian))
        initialize(eg(:q_m*sec), vague(Gaussian))
        initialize(eg(:q_gam*sec), vague(Gamma))
    end

    # Construct algorithm
    algo = VariationalBayes(n_iterations=50)

    # Verify correct schedule lengths
    @fact length(algo.recognition_factorization.subgraphs[1].internal_pre_schedule) --> 2
    @fact length(algo.recognition_factorization.subgraphs[1].internal_iterative_schedule) --> 6
    @fact length(algo.recognition_factorization.subgraphs[1].internal_post_schedule) --> 1
    @fact length(algo.recognition_factorization.subgraphs[2].internal_pre_schedule) --> 2
    @fact length(algo.recognition_factorization.subgraphs[2].internal_iterative_schedule) --> 6
    @fact length(algo.recognition_factorization.subgraphs[2].internal_post_schedule) --> 1
    @fact length(algo.recognition_factorization.subgraphs[3].internal_pre_schedule) --> 1
    @fact length(algo.recognition_factorization.subgraphs[3].internal_iterative_schedule) --> 0 # Message is skipped
    @fact length(algo.recognition_factorization.subgraphs[3].internal_post_schedule) --> 0
    @fact length(algo.recognition_factorization.subgraphs[4].internal_pre_schedule) --> 1
    @fact length(algo.recognition_factorization.subgraphs[4].internal_iterative_schedule) --> 0 # Message is skipped
    @fact length(algo.recognition_factorization.subgraphs[4].internal_post_schedule) --> 0
end

facts("VariationalBayes should collect the proper inbound types as dependent on the recognition factorization") do
    # Mean field factorized Gaussian node
    initializeGaussianNode()

    rf = RecognitionFactorization()
    factorizeMeanField()
    initialize(eg(:edge1), vague(Gaussian))
    initialize(eg(:edge2), vague(Gamma))
    initialize(eg(:edge3), vague(Gaussian))

    algo = VariationalBayes()

    @fact algo.recognition_factorization.subgraphs[3].internal_pre_schedule[1].outbound_type --> Gaussian
    @fact algo.recognition_factorization.subgraphs[2].internal_pre_schedule[1].outbound_type --> Gamma
    @fact algo.recognition_factorization.subgraphs[1].internal_pre_schedule[1].outbound_type --> Gaussian

    @fact algo.recognition_factorization.subgraphs[3].internal_iterative_schedule[1].inbound_types --> [Gaussian, Gamma, Void]
    @fact algo.recognition_factorization.subgraphs[3].internal_iterative_schedule[1].outbound_type --> Gaussian
    @fact algo.recognition_factorization.subgraphs[2].internal_iterative_schedule[1].inbound_types --> [Gaussian, Void, Gaussian]
    @fact algo.recognition_factorization.subgraphs[2].internal_iterative_schedule[1].outbound_type --> Gamma
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[1].inbound_types --> [Void, Gamma, Gaussian]
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[1].outbound_type --> Gaussian

    # Structurally factorized
    initializeGaussianNode()

    rf = RecognitionFactorization()
    factor([eg(:edge1), eg(:edge2)])
    factor(eg(:edge3))
    initialize([eg(:edge1), eg(:edge2)], vague(NormalGamma))
    initialize(eg(:edge3), vague(Gaussian))

    algo = VariationalBayes()
    
    @fact algo.recognition_factorization.subgraphs[2].internal_pre_schedule[1].outbound_type --> Gaussian
    @fact algo.recognition_factorization.subgraphs[1].internal_pre_schedule[1].outbound_type --> Gaussian
    @fact algo.recognition_factorization.subgraphs[1].internal_pre_schedule[2].outbound_type --> Gamma

    @fact algo.recognition_factorization.subgraphs[2].internal_iterative_schedule[1].inbound_types --> [NormalGamma, NormalGamma, Void]
    @fact algo.recognition_factorization.subgraphs[2].internal_iterative_schedule[1].outbound_type --> Gaussian
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[1].inbound_types --> [Void,Message{Gamma},Gaussian]
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[1].outbound_type --> StudentsT
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[2].inbound_types --> [Message{Gaussian},Void,Gaussian]
    @fact algo.recognition_factorization.subgraphs[1].internal_iterative_schedule[2].outbound_type --> Gamma
end

facts("Naive vmp implementation integration tests") do
    context("Call step() for VMP algorithm") do
        data = [Gaussian(m=2.0, V=tiny)]
        g = FactorGraph()
        GaussianNode(form=:precision,id=:g_node)
        TerminalNode(id=:t_out)
        TerminalNode(Gaussian(), id=:t_mean)
        TerminalNode(Gamma(), id=:t_prec)
        Edge(n(:g_node).i[:out], n(:t_out), id=:y)
        Edge(n(:t_mean), n(:g_node).i[:mean], id=:m)
        Edge(n(:t_prec), n(:g_node).i[:precision], id=:gam)

        attachReadBuffer(n(:t_out), data)
        mean_out = attachWriteBuffer(n(:g_node).i[:mean].edge)
        prec_out = attachWriteBuffer(n(:g_node).i[:precision].edge)

        rf = RecognitionFactorization()
        factorizeMeanField()
        initialize(eg(:y), vague(Gaussian))
        initialize(eg(:m), vague(Gaussian))
        initialize(eg(:gam), vague(Gamma))

        algo = VariationalBayes(n_iterations=10)
        prepare!(algo)
        step(algo)

        ForneyLab.ensureParameters!(mean_out[end], (:xi, :W))
        @fact round(mean_out[end].W[1,1], 2) --> 1.79
        @fact round(mean_out[end].xi[1], 2) --> 1.57
        @fact prec_out[end].a --> 1.5
        @fact round(prec_out[end].b, 2) --> 1.91
    end

    context("Gaussian node mean precision batch estimation") do
        # Integration test for the vmp implementation by trying to estimate the mean and precision of a Gaussian

        # Initialize chain
        # Fixed observations drawn from N(5.0, 2.0)
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        initializeGaussianNodeChain(data)
        n_sections = length(data)

        m_buffer = attachWriteBuffer(n(:m_eq*n_sections).i[2]) # Propagate to the end
        gam_buffer = attachWriteBuffer(n(:gam_eq*n_sections).i[2])

        # Apply mean field factorization
        rf = RecognitionFactorization()
        factorizeMeanField()
        for sec in 1:n_sections
            initialize(eg(:q_y*sec), vague(Gaussian))
            initialize(eg(:q_m*sec), vague(Gaussian))
            initialize(eg(:q_gam*sec), vague(Gamma))
        end

        # Construct algorithm
        algo = VariationalBayes(n_iterations=50)

        # Perform vmp updates
        run(algo)

        m_out = m_buffer[end]
        gam_out = gam_buffer[end]

        @fact round(mean(m_out)[1], 3) --> 4.381
        @fact round(var(m_out)[1, 1], 5) --> 0.08313
        @fact round(gam_out.a, 3)  --> 6.000
        @fact round(1/gam_out.b, 4) --> 0.2005 # Scale
    end

    context("Gaussian node mean precision batch estimation for multivariate distributions") do
        # Fixed observations drawn from N(5.0, 2.0)
        data_1 = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        # Fixed observations drawn from N(2.0, 0.5)
        data_2 = [2.317991577739536,2.8244768760034407,1.2501129068467165,2.5729664094889424,3.05374531248249,1.5149856277603246,2.3119227037528614,2.0264643318813644,1.6248999854457278,0.7425070466631876]
        data = hcat(data_1, data_2)
        
        initializeMvGaussianNodeChain(data)
        n_sections = size(data, 1)

        m_buffer = attachWriteBuffer(n(:m_eq*n_sections).i[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq*n_sections).i[2])

        # Apply mean field factorization
        rf = RecognitionFactorization()
        factorizeMeanField()
        for sec in 1:n_sections
            initialize(eg(:q_y*sec), vague(MvGaussian{2}))
            initialize(eg(:q_m*sec), vague(MvGaussian{2}))
            initialize(eg(:q_gam*sec), vague(Wishart{2}))
        end

        # Construct algorithm
        algo = VariationalBayes(n_iterations=50)

        # Perform vmp updates
        run(algo)

        m_out = m_buffer[end]
        gam_out = gam_buffer[end]

        @fact round(mean(m_out), 5) --> [4.38083, 2.02401]
        @fact round(cov(m_out), 5) --> [0.15241 -0.03349; -0.03349 0.08052]
        @fact round(mean(gam_out), 5)  --> [1.03159 0.42907; 0.42907 1.95259]
        @fact round(var(gam_out), 5) --> [0.21283  0.21984; 0.21984 0.76252]
    end

    context("Gaussian node mean precision online estimation") do
        # Integration test for the vmp implementation by trying to estimate the mean and precision of a Gaussian

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
        rf = RecognitionFactorization()
        factorizeMeanField()
        initialize(eg(:q_y1), vague(Gaussian))
        initialize(eg(:q_m1), vague(Gaussian))
        initialize(eg(:q_gam1), vague(Gamma))

        # Construct algorithm
        algo = VariationalBayes(n_iterations=50)

        # Perform vmp updates
        run(algo)

        m_out = m_buffer[end]
        gam_out = gam_buffer[end]
        @fact round(mean(m_out)[1], 3) --> 4.941
        @fact round(var(m_out)[1, 1], 3) --> 0.000
        @fact round(gam_out.a, 3) --> 6.000
        @fact round(1/gam_out.b, 4) --> 0.1628 # Scale
    end
end

facts("Structured vmp implementation integration tests") do
    context("Gaussian node joint mean variance estimation") do
        # Initialize chain
        # Samples drawn from N(mean 5.0, prec 0.5): data = randn(100)*(1/sqrt(0.5))+5.0
        data = [4.9411489951651735,4.4083330961647595,3.535639074214823,2.1690761263145855,4.740705436131505,5.407175878845115,3.6458623443189957,5.132115496214244,4.485471215629411,5.342809672818667]
        d_data = [Delta(d_k) for d_k in data]
        # d_data = [Gaussian(m=d_k, W=10.0) for d_k in data]
        initializeGaussianNodeChain([0.0]) # Initialize a length 1 chain
        Wrap(n(:mN), n(:m0))
        Wrap(n(:gamN), n(:gam0))

        attachReadBuffer(n(:y1), d_data)
        m_buffer = attachWriteBuffer(n(:m_eq1).i[2])
        gam_buffer = attachWriteBuffer(n(:gam_eq1).i[2])

        # Structured factorization
        rf = RecognitionFactorization()
        factor([eg(:q_m1), eg(:q_gam1)])
        factor(eg(:q_y1))
        initialize([eg(:q_m1), eg(:q_gam1)], vague(NormalGamma))
        initialize(eg(:q_y1), vague(Gaussian))

        # Construct algorithm
        algo = VariationalBayes(n_iterations=10)

        run(algo)

        m_out = m_buffer[end]
        gam_out = gam_buffer[end]
        # Reference values from first run
        @fact round(mean(m_out)[1], 3) --> 4.521
        @fact round(var(m_out)[1, 1], 3) --> 0.873 # Uniform gamma priors make the variance collapse
        @fact round(gam_out.a, 3) --> 6.000
        @fact round(1/gam_out.b, 5) --> 0.04549 # Scale
    end
end