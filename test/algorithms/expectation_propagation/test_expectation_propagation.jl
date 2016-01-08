facts("ExpectationPropagation algorithm integration tests") do
    # TODO: Still some mistake around here

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
    # Build graph
    ###############

    const NUM_SECTIONS = 50
    g = FactorGraph()
    X_prior = TerminalNode(GaussianDistribution(m=0.0, V=20.0), id=:t_X_prior)
    prev_section_connector = X_prior.i[:out]
    sites = Vector{Interface}()
    for section = 1:NUM_SECTIONS
        equ = EqualityNode(id=symbol("equ$(section)"))
        sig = SigmoidNode(id=symbol("sig$(section)"))
        y = TerminalNode(DeltaDistribution(false), id=symbol("t_Y$(section)"))
        Edge(prev_section_connector, equ.interfaces[1], GaussianDistribution)
        Edge(equ.interfaces[2], sig.i[:real], GaussianDistribution, id=symbol("X$(section)"))
        Edge(sig.i[:bin], y, BernoulliDistribution)
        prev_section_connector = equ.interfaces[3]
        push!(sites, sig.i[:real])
    end
    Edge(prev_section_connector, TerminalNode(vague(GaussianDistribution), id=:t_X_term), GaussianDistribution, id=:X_marg) # Terminate equality chain

    # ForneyLab.draw(g)

    ###############
    # Generate data
    ###############

    generating_distribution = GaussianDistribution(m=-0.5, V=0.1)
    edge(:X1).tail.message = Message(generating_distribution)
    run(SumProduct(n(:sig1).i[:bin]))
    bin_dist = n(:sig1).i[:bin].message.payload
    samples = Vector{Float64}()
    for section = 1:NUM_SECTIONS
        s = sample(bin_dist)
        push!(samples, float(s))
        n(symbol("t_Y$(section)")).value = DeltaDistribution(s)
    end
    clearMessages!(g)


    ###############
    # EP algorithm
    ###############
    #setVerbosity(true)

    # Define callback function to collect the posterior dist. of X after each iteration
    X_marg = Vector{GaussianDistribution}()
    function log_results()
        run(SumProduct(edge(:X_marg).tail))
        push!(X_marg, deepcopy(ensureParameters!(edge(:X_marg).tail.message.payload, (:m, :V))))
        return (length(X_marg) >= 5) # terminate the EP algorithm after 5 iterations
    end

    ep_algo = ExpectationPropagation(sites, num_iterations=10, callback=log_results)

    context("ExpectationPropagation construction") do
        @fact ep_algo.num_iterations --> 10
        @fact is(ep_algo.callback, log_results) --> true
        @fact typeof(ep_algo.schedule) --> Schedule
        for schedule_entry in ep_algo.schedule
            if schedule_entry.interface in sites
                @fact is(schedule_entry.message_calculation_rule, ep!) --> true
            else
                @fact is(schedule_entry.message_calculation_rule, sumProduct!) --> true
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
        run(SumProduct(n(:sig1).i[:bin]))
        predictive = n(:sig1).i[:bin].message.payload
        @fact predictive.p --> roughly(mean(samples), atol=0.1)
    end
end
