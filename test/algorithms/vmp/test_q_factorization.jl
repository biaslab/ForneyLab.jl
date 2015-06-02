#####################
# Integration tests
#####################

facts("QFactorization integration tests") do
    context("extend() should extend a set of edges to envelope deterministic nodes") do
        initializeFactoringGraph()
        cluster = VMP.extend(Set{Edge}([n(:g1).i[:out].edge, n(:g1).i[:mean].edge]))
        @fact cluster => Set{Edge}({n(:t1).i[:out].edge, n(:a1).i[:out].edge, n(:g1).i[:out].edge, n(:add1).i[:out].edge, n(:g2).i[:out].edge})
    end

    context("factorize!()") do
        context("Should include argument edges in a new subgraph") do
            initializeFactoringGraph()
            f = VMP.factorize!(Set{Edge}([n(:t2).i[:out].edge]))
            @fact f.factors[2].internal_edges => Set{Edge}({n(:t2).i[:out].edge})
            @fact f.factors[1].internal_edges => Set{Edge}({n(:t1).i[:out].edge, n(:a1).i[:out].edge, n(:g1).i[:out].edge, n(:add1).i[:out].edge, n(:g2).i[:out].edge})
        end

        context("Should update the edge_to_subgraph mapping for the graph") do
            initializeFactoringGraph()
            f = VMP.factorize!(Set{Edge}([n(:t2).i[:out].edge]))
            @fact f.edge_to_subgraph[n(:t1).i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[n(:a1).i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[n(:g1).i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[n(:t2).i[:out].edge] => f.factors[2]
            @fact f.edge_to_subgraph[n(:add1).i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[n(:g2).i[:out].edge] => f.factors[1]
        end

        context("Should not factorize a GaussianNode with fixed mean and variance") do
            initializeGaussianFactoringGraph()
            f = VMP.factorize!(Set{Edge}([n(:t).i[:out].edge]))
            @fact length(f.factors) => 1
            @fact f.edge_to_subgraph[n(:t).i[:out].edge] => f.factors[1]
            @fact f.edge_to_subgraph[n(:gauss).i[:out].edge] => f.factors[1]
        end

        context("Should retain node ids") do
            initializeSimpleFactoringGraph()
            f = VMP.factorize!(Set{Edge}([n(:t1).i[:out].edge]))
            @fact node(:t1).id => :t1
        end
    end

    context("factorize() should output a mean field factorized graph") do
        data = [1.0, 1.0, 1.0]
        initializeGaussianNodeChain(data)
        f = VMP.factorize()
        gam_set = Set{Edge}()
        for gam_eq_node in [n(:gam_eq1), n(:gam_eq2), n(:gam_eq3)]
            for interface in gam_eq_node.interfaces
                push!(gam_set, interface.edge)
            end
        end
        m_set = Set{Edge}()
        for m_eq_node in [n(:m_eq1), n(:m_eq2), n(:m_eq3)]
            for interface in m_eq_node.interfaces
                push!(m_set, interface.edge)
            end
        end
        @fact length(f.factors) => 5

        @fact f.edge_to_subgraph[n(:g1).i[:mean].edge].internal_edges => m_set 
        @fact f.edge_to_subgraph[n(:g1).i[:precision].edge].internal_edges => gam_set 
        @fact f.edge_to_subgraph[n(:g1).i[:out].edge].internal_edges => Set{Edge}([e(:q_y1)]) 
        @fact f.edge_to_subgraph[n(:g2).i[:out].edge].internal_edges => Set{Edge}([e(:q_y2)])
        @fact f.edge_to_subgraph[n(:g3).i[:out].edge].internal_edges => Set{Edge}([e(:q_y3)])
    end
end

facts("vagueQDistributions() should set vague marginals at the appropriate places") do
    data = [1.0, 1.0, 1.0]

    # MF case
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    f = VMP.factorize()
    VMP.generateSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = VMP.vagueQDistributions(f) # Initialize vague q distributions

    m_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
    gam_subgraph = f.edge_to_subgraph[n(:g1).i[:precision].edge]
    y1_subgraph = f.edge_to_subgraph[n(:g1).i[:out].edge]
    y2_subgraph = f.edge_to_subgraph[n(:g2).i[:out].edge]
    y3_subgraph = f.edge_to_subgraph[n(:g3).i[:out].edge]

    @fact length(qs) => 9
    @fact qs[(n(:g1), m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(n(:g2), m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(n(:g3), m_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(n(:g1), gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(n(:g2), gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(n(:g3), gam_subgraph)].distribution => vague(GammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(n(:g2), y2_subgraph)].distribution => vague(GaussianDistribution)
    @fact qs[(n(:g3), y3_subgraph)].distribution => vague(GaussianDistribution)

    # Structured case
    data = [1.0]
    initializeGaussianNodeChain(data)
    n_sections = length(data)

    f = VMP.QFactorization()
    for edge in [e(:q_y1)]
        f = VMP.factorize!(Set{Edge}({edge}), f)
    end
    VMP.generateSchedule!(f) # Generate and store internal and external schedules on factorization subgraphs
    qs = VMP.vagueQDistributions(f)

    m_gam_subgraph = f.edge_to_subgraph[n(:g1).i[:mean].edge]
    y1_subgraph = f.edge_to_subgraph[n(:g1).i[:out].edge]

    @fact length(qs) => 2
    @fact qs[(n(:g1), m_gam_subgraph)].distribution => vague(NormalGammaDistribution)
    @fact qs[(n(:g1), y1_subgraph)].distribution => vague(GaussianDistribution)
end