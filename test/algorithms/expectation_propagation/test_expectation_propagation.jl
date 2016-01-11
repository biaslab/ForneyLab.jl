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
    samples = attachWriteBuffer(Y.i[:out].partner)
    algo = SumProduct(post_processing_functions=Dict(sig.i[:bin] => sample))
    run(algo)

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
        Edge(prev_section_connector, equ.interfaces[1], GaussianDistribution)
        Edge(equ.interfaces[2], sig.i[:real], GaussianDistribution, id=symbol("X$(section)"))
        Edge(sig.i[:bin], y, BernoulliDistribution)
        prev_section_connector = equ.interfaces[3]
        push!(sites, (sig.i[:real], GaussianDistribution))
    end
    Edge(prev_section_connector, TerminalNode(vague(GaussianDistribution), id=:t_X_term), GaussianDistribution, id=:X_marg) # Terminate equality chain

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

    ep_algo = ExpectationPropagation(sites, num_iterations=10, callback=log_results)

    sitelist = [site for (site, distribution) in sites]
    context("ExpectationPropagation construction") do
        @fact ep_algo.num_iterations --> 10
        @fact is(ep_algo.callback, log_results) --> true
        @fact typeof(ep_algo.schedule) --> Schedule
        for entry in ep_algo.schedule
            if entry.node.interfaces[entry.outbound_interface_id] in sitelist
                @fact is(entry.rule, ep!) --> true
            else
                @fact is(entry.rule, sumProduct!) --> true
            end
        end
    end

    #show(ep_algo.schedule)
    run(ep_algo)

    ###############
    # Compare sample mean with predictive mean after inference
    ###############

    context("Inference results") do
        @fact length(X_marg) --> 5
        predictive = BernoulliDistribution()
        sumProduct!(n(:sig1), Val{2}, n(:equ1).interfaces[2].message, nothing, predictive)
        @fact predictive.p --> roughly(mean([sample.m for sample in samples]), atol=0.1)
    end
end
