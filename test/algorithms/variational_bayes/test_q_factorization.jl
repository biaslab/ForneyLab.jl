#####################
# Integration tests
#####################

facts("QFactorization integration tests") do
    context("extend() should extend a set of edges to envelope deterministic nodes") do
        initializeFactoringGraph()
        cluster = ForneyLab.extend(Set{Edge}([n(:g1).i[:out].edge, n(:g1).i[:mean].edge]))
        @fact cluster --> Set{Edge}(Edge[n(:t1).i[:out].edge, n(:a1).i[:out].edge, n(:g1).i[:out].edge, n(:add1).i[:out].edge, n(:g2).i[:out].edge])
    end

    context("factorize!()") do
        context("Should include argument edges in a new subgraph") do
            initializeFactoringGraph()
            f = ForneyLab.factorize!(Set{Edge}([n(:t2).i[:out].edge]))
            @fact f.factors[2].internal_edges --> Set{Edge}(Edge[n(:t2).i[:out].edge])
            @fact f.factors[1].internal_edges --> Set{Edge}(Edge[n(:t1).i[:out].edge, n(:a1).i[:out].edge, n(:g1).i[:out].edge, n(:add1).i[:out].edge, n(:g2).i[:out].edge])
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            initializeFactoringGraph()
            f = ForneyLab.factorize!(Set{Edge}([n(:t2).i[:out].edge]))
            @fact f.edge_to_subgraph[n(:t1).i[:out].edge] --> f.factors[1]
            @fact f.edge_to_subgraph[n(:a1).i[:out].edge] --> f.factors[1]
            @fact f.edge_to_subgraph[n(:g1).i[:out].edge] --> f.factors[1]
            @fact f.edge_to_subgraph[n(:t2).i[:out].edge] --> f.factors[2]
            @fact f.edge_to_subgraph[n(:add1).i[:out].edge] --> f.factors[1]
            @fact f.edge_to_subgraph[n(:g2).i[:out].edge] --> f.factors[1]
        end

        context("Should not factorize a GaussianNode with fixed mean and variance") do
            initializeGaussianFactoringGraph()
            f = ForneyLab.factorize!(Set{Edge}([n(:t).i[:out].edge]))
            @fact length(f.factors) --> 1
            @fact f.edge_to_subgraph[n(:t).i[:out].edge] --> f.factors[1]
            @fact f.edge_to_subgraph[n(:gauss).i[:out].edge] --> f.factors[1]
        end

        context("Should retain node ids") do
            initializeSimpleFactoringGraph()
            f = ForneyLab.factorize!(Set{Edge}([n(:t1).i[:out].edge]))
            @fact node(:t1).id --> :t1
        end
    end

    context("factorize() should output a mean field factorized graph") do
        data = Float64[1.0, 1.0, 1.0]
        initializeGaussianNodeChain(data)
        f = ForneyLab.factorize(Dict(
                eg(:q_m*(1:3)) => GaussianDistribution,
                eg(:q_gam*(1:3)) => GammaDistribution,
                eg(:q_y*(1:3)) => GaussianDistribution))

        gam_set = Set{Edge}()
        for gam_eq_node in n(:gam_eq*(1:3))
            for interface in gam_eq_node.interfaces
                push!(gam_set, interface.edge)
            end
        end
        m_set = Set{Edge}()
        for m_eq_node in n(:m_eq*(1:3))
            for interface in m_eq_node.interfaces
                push!(m_set, interface.edge)
            end
        end
        @fact length(f.factors) --> 5

        @fact f.edge_to_subgraph[n(:g1).i[:mean].edge].internal_edges --> m_set
        @fact f.edge_to_subgraph[n(:g1).i[:precision].edge].internal_edges --> gam_set
        @fact f.edge_to_subgraph[n(:g1).i[:out].edge].internal_edges --> Set{Edge}([eg(:q_y1)])
        @fact f.edge_to_subgraph[n(:g2).i[:out].edge].internal_edges --> Set{Edge}([eg(:q_y2)])
        @fact f.edge_to_subgraph[n(:g3).i[:out].edge].internal_edges --> Set{Edge}([eg(:q_y3)])
    end
end

facts("initializeVagueQDistributions() should set vague marginals at the appropriate places") do
    data = [1.0, 1.0, 1.0]

    # MF case
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    recognition_distribution_types = Dict(
        eg(:q_m*(1:3)) => GaussianDistribution,
        eg(:q_gam*(1:3)) => GammaDistribution,
        eg(:q_y*(1:3)) => GaussianDistribution)

    f = ForneyLab.factorize(recognition_distribution_types)

    ForneyLab.generateVariationalBayesSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = ForneyLab.initializeVagueQDistributions(f, recognition_distribution_types) # Initialize vague q distributions

    m_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
    gam_subgraph = f.edge_to_subgraph[n(:g1).i[:precision].edge]
    y1_subgraph = f.edge_to_subgraph[n(:g1).i[:out].edge]
    y2_subgraph = f.edge_to_subgraph[n(:g2).i[:out].edge]
    y3_subgraph = f.edge_to_subgraph[n(:g3).i[:out].edge]

    @fact length(qs) --> 9
    @fact qs[(n(:g1), m_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g2), m_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g3), m_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g1), gam_subgraph)].distribution --> vague(GammaDistribution)
    @fact qs[(n(:g2), gam_subgraph)].distribution --> vague(GammaDistribution)
    @fact qs[(n(:g3), gam_subgraph)].distribution --> vague(GammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g2), y2_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g3), y3_subgraph)].distribution --> vague(GaussianDistribution)

    # Structured case
    data = [1.0]
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    recognition_distribution_types = Dict(
        [eg(:q_m1), eg(:q_gam1)].' => NormalGammaDistribution,
        eg(:q_y1) => GaussianDistribution)

    f = ForneyLab.factorize(recognition_distribution_types)

    ForneyLab.generateVariationalBayesSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = ForneyLab.initializeVagueQDistributions(f, recognition_distribution_types)

    m_gam_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
    y1_subgraph = f.edge_to_subgraph[n(:g1).i[:out].edge]

    @fact length(qs) --> 2
    @fact qs[(n(:g1), m_gam_subgraph)].distribution --> vague(NormalGammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution --> vague(GaussianDistribution)
end

facts("resetQDistributions!() should reset already present q distributions to vague") do
    data = [1.0]

    # MF case
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    algo = VariationalBayes(Dict(
                eg(:q_m1) => GaussianDistribution,
                eg(:q_gam1) => GammaDistribution,
                eg(:q_y1) => GaussianDistribution))

    f = algo.factorization
    qs = algo.q_distributions

    m_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
    gam_subgraph = f.edge_to_subgraph[n(:g1).i[:precision].edge]
    y1_subgraph = f.edge_to_subgraph[n(:g1).i[:out].edge]

    # Distributions before step() should be vague
    @fact qs[(n(:g1), m_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g1), gam_subgraph)].distribution --> vague(GammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution --> vague(GaussianDistribution)

    prepare!(algo)
    step(algo)

    # Distributions after step() should not be vague
    @fact qs[(n(:g1), m_subgraph)].distribution == vague(GaussianDistribution) --> false
    @fact qs[(n(:g1), gam_subgraph)].distribution == vague(GammaDistribution) --> false
    @fact qs[(n(:g1), y1_subgraph)].distribution == vague(GaussianDistribution) --> false

    ForneyLab.resetQDistributions!(qs)

    # Distributions after resetting should be vague again
    @fact qs[(n(:g1), m_subgraph)].distribution --> vague(GaussianDistribution)
    @fact qs[(n(:g1), gam_subgraph)].distribution --> vague(GammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution --> vague(GaussianDistribution)
end
