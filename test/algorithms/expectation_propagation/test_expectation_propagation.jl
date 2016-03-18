facts("ExpectationPropagation algorithm integration tests") do
    # Perform expectation propagation in a graph consisting
    # of SigmoidNode sections:
    #
    # ... ---[=]--- ... X
    #         |
    #         |
    #        [σ]
    #         |
    #         |  Y
    #
    # We assume X ~ 𝓝(μ,σ2) and we observe Y ∈ [false, true]
    # The goal is to infer X from the Y observations.

    ###############
    # Generate data
    ###############

    const NUM_SECTIONS = 50
    g_gen = FactorGraph()
    X = TerminalNode(GaussianDistribution(m=-0.5, V=0.1), id=:t_X)
    sig = SigmoidNode(id=:sig)
    Y = TerminalNode(DeltaDistribution(false), id=:t_Y)
    Edge(X, sig.i[:real])
    Edge(sig.i[:bin], Y)
    generating_distributions = [GaussianDistribution(m=-0.5, V=0.1) for i = 1:NUM_SECTIONS]
    attachReadBuffer(X, generating_distributions)
    y_predictions = attachWriteBuffer(Y.i[:out].partner)
    algo = SumProduct()
    run(algo)
    samples = map(sample, y_predictions)

    ###############
    # Build graph
    ###############

    g = FactorGraph()
    X_prior = TerminalNode(GaussianDistribution(m=0.0, V=20.0), id=:t_X_prior)
    prev_section_connector = X_prior.i[:out]
    sites = Vector{Tuple{Interface, DataType}}()
    for section = 1:NUM_SECTIONS
        equ = EqualityNode(id=symbol("equ$(section)"))
        sig = SigmoidNode(id=symbol("sig$(section)"))
        y = TerminalNode(samples[section], id=symbol("t_Y$(section)"))
        Edge(prev_section_connector, equ.interfaces[1])
        Edge(equ.interfaces[2], sig.i[:real], id=symbol("X$(section)"))
        Edge(sig.i[:bin], y)
        prev_section_connector = equ.interfaces[3]
        push!(sites, (sig.i[:real], GaussianDistribution))
    end
    Edge(prev_section_connector, TerminalNode(vague(GaussianDistribution), id=:t_X_term), id=:X_marg) # Terminate equality chain

    # ForneyLab.draw(g)

    ###############
    # EP algorithm
    ###############
    #setVerbosity(true)

    # Define callback function to collect the posterior dist. of X after each iteration for later reference
    X_marg = Vector{GaussianDistribution}()
    function log_results()
        push!(X_marg, deepcopy(ensureParameters!(edge(:X1).marginal, (:m, :V))))
        return (length(X_marg) >= 5) # terminate the EP algorithm after 5 iterations
    end

    ep_algo = ExpectationPropagation(prev_section_connector, sites, n_iterations=10, callback=log_results)

    sitelist = [site for (site, distribution) in sites]
    context("ExpectationPropagation construction") do
        @fact ep_algo.n_iterations --> 10
        @fact is(ep_algo.callback, log_results) --> true
        @fact typeof(ep_algo.iterative_schedule) --> Schedule
        for entry in ep_algo.iterative_schedule
            if entry.node.interfaces[entry.outbound_interface_id] in sitelist
                @fact is(entry.rule, expectationRule!) --> true
            else
                @fact is(entry.rule, sumProductRule!) --> true
            end
        end
        @fact typeof(ep_algo.post_convergence_schedule) --> Schedule
        @fact length(ep_algo.post_convergence_schedule) --> 1
        for entry in ep_algo.post_convergence_schedule
            @fact is(entry.rule, sumProductRule!) --> true
        end
    end

    run(ep_algo)

    ###############
    # Compare sample mean with predictive mean after inference
    ###############

    context("Inference results") do
        @fact length(X_marg) --> 5
        predictive = BernoulliDistribution()
        sumProductRule!(n(:sig1), Val{2}, predictive, n(:equ1).interfaces[2].message, nothing)
        @fact predictive.p --> roughly(mean(samples), atol=0.1)
    end


    ###############
    # EP in a graph with wraps
    ###############
    FactorGraph()
    EqualityNode(id=:eq)
    SigmoidNode(id=:sig)
    TerminalNode(DeltaDistribution(false), id=:Y)
    PriorNode(GaussianDistribution(m=-0.5, V=0.1), id=:X0)
    TerminalNode(vague(GaussianDistribution), id=:XN)
    Y_dists = ProbabilityDistribution[samples[section] for section in 1:NUM_SECTIONS]
    Y_buffer = attachReadBuffer(n(:Y), deepcopy(Y_dists))
    # Connect
    Edge(n(:X0),n(:eq).i[1])
    Edge(n(:eq).i[2], n(:sig).i[:real], id=:X)
    Edge(n(:eq).i[3], n(:XN))
    Edge(n(:sig).i[:bin], n(:Y))
    Wrap(n(:XN),n(:X0), block_size=NUM_SECTIONS)

    sites = Vector{Tuple{Interface, DataType}}()
    push!(sites, (n(:sig).i[:real], GaussianDistribution))

    context("Throws an error when the wrap is not implemented correct") do
            algo = ExpectationPropagation(currentGraph(), sites; n_iterations=1)
    end

    sitelist = [site for (site, distribution) in sites]
    context("ExpectationPropagation construction with wraps") do
        @fact algo.n_iterations --> 1
        @fact typeof(algo.iterative_schedule) --> Schedule
        @fact length(algo.iterative_schedule) --> 5
        for entry in algo.iterative_schedule
            if entry.node.interfaces[entry.outbound_interface_id] in sitelist
                @fact is(entry.rule, expectationRule!) --> true
            else
                @fact is(entry.rule, sumProductRule!) --> true
            end
        end
        @fact typeof(algo.post_convergence_schedule) --> Schedule
        @fact length(algo.post_convergence_schedule) --> 2
        for entry in algo.post_convergence_schedule
            @fact is(entry.rule, sumProductRule!) --> true
        end
    end



end
