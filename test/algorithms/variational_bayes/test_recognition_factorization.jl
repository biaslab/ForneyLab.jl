#####################
# Unit tests
#####################

facts("RecognitionFactorization(::FactorGraph) should initialize a new recognition factorization") do
    g = FactorGraph()
    rf = RecognitionFactorization(g)
    @fact rf.graph --> g
    @fact rf.subgraphs --> Dict{Symbol, Subgraph}()
    @fact rf.subgraph_id_counter --> 0
    @fact rf.edge_to_subgraph --> Dict{Edge, Subgraph}()
    @fact rf.node_subgraph_to_internal_edges --> Dict{Tuple{Node, Subgraph}, Set{Edge}}()
    @fact rf.initial_recognition_distributions --> Partitioned{ProbabilityDistribution,0}(ProbabilityDistribution[])
    @fact rf.recognition_distributions --> Partitioned{ProbabilityDistribution,0}(ProbabilityDistribution[])
    @fact rf.node_subgraph_to_recognition_distribution --> Dict{Tuple{Node, Subgraph}, ProbabilityDistribution}()
    @fact rf.default_update_order --> Subgraph[]
end

facts("RecognitionFactorization() should initialize a new recognition factorization and create a current factorization") do
    g = FactorGraph()
    rf = RecognitionFactorization()
    @fact is(currentRecognitionFactorization(), rf) --> true
end


#####################
# Integration tests
#####################

facts("extend() should extend a set of edges to envelope deterministic nodes") do
    #               [T2]   [T3]
    #                |      | 
    #                v      v 
    # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

    g = initializeFactoringGraph()

    @fact ForneyLab.extend(eg(:t1_a1)) --> Set([eg(:t1_a1), eg(:a1_g1)])
    @fact ForneyLab.extend(eg(:g1_g2)) --> Set([eg(:g1_g2)])
end

facts("factor() should specify a new custom subgraph for the recognition factorization and initialize lookup tables") do
    context("Mean field case") do
        #               [T2]   [T3]
        #     here       |      |
        #      v         v      v
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

        g = initializeFactoringGraph()
        rf = RecognitionFactorization()
        factor(eg(:t1_a1)) # Should automatically bind to current graph and create a RecognitionFactorization instance

        # Verify factorization properties
        @fact typeof(rf) --> RecognitionFactorization
        @fact rf.graph --> g
        @fact length(rf.subgraphs) --> 1
        @fact rf.subgraph_id_counter --> 1

        sg1 = sg(:subgraph1)
        @fact length(rf.edge_to_subgraph) --> 2
        @fact rf.edge_to_subgraph[eg(:t1_a1)] --> sg1
        @fact rf.edge_to_subgraph[eg(:a1_g1)] --> sg1

        @fact length(rf.node_subgraph_to_internal_edges) --> 1
        @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg1)] --> Set([eg(:a1_g1)])

        # Verify subgraph properties
        @fact sg1.internal_edges --> Set([eg(:t1_a1), eg(:a1_g1)])
        @fact sg1.external_edges --> Set([eg(:t2_g1), eg(:g1_g2)])
        @fact sg1.nodes_connected_to_external_edges --> [n(:g1)]

        #               [T2]   [T3]
        #                |      |
        #                v      v
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]
        #                    ^
        #                   here

        factor(eg(:g1_g2), id=:second_sg)

        # Verify factorization properties
        @fact length(rf.subgraphs) --> 2 # A subgraph should be added
        @fact rf.subgraph_id_counter --> 1 # Custom name, so identifier counter should still be one

        sg2 = sg(:second_sg)
        @fact length(rf.edge_to_subgraph) --> 3
        @fact rf.edge_to_subgraph[eg(:g1_g2)] --> sg2

        @fact length(rf.node_subgraph_to_internal_edges) --> 3
        @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg1)] --> Set([eg(:a1_g1)])
        @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg2)] --> Set([eg(:g1_g2)])
        @fact rf.node_subgraph_to_internal_edges[(n(:g2), sg2)] --> Set([eg(:g1_g2)])

        # Verify subgraph properties
        @fact sg2.internal_edges --> Set([eg(:g1_g2)])
        @fact sg2.external_edges --> Set([eg(:a1_g1), eg(:t2_g1), eg(:t3_g2), eg(:g2_t4)])
        @fact sg2.nodes_connected_to_external_edges --> [n(:g1), n(:g2)]
    end

    context("Structured case") do
        #               [T2]   [T3]
        #                |      | < here
        #                v      v   v
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

        g = initializeFactoringGraph()
        rf = RecognitionFactorization()
        factor(Set([eg(:t3_g2), eg(:g2_t4)])) # Should automatically bind to current graph and create a RecognitionFactorization instance

        # Verify factorization properties
        @fact length(rf.subgraphs) --> 1

        sg1 = sg(:subgraph1)
        @fact length(rf.edge_to_subgraph) --> 2
        @fact rf.edge_to_subgraph[eg(:t3_g2)] --> sg1
        @fact rf.edge_to_subgraph[eg(:g2_t4)] --> sg1

        @fact length(rf.node_subgraph_to_internal_edges) --> 1
        @fact rf.node_subgraph_to_internal_edges[(n(:g2), sg1)] --> Set([eg(:t3_g2), eg(:g2_t4)])

        # Verify subgraph properties
        @fact sg1.internal_edges --> Set([eg(:t3_g2), eg(:g2_t4)])
        @fact sg1.external_edges --> Set([eg(:g1_g2)])
        @fact sg1.nodes_connected_to_external_edges --> [n(:g2)]
    end
end

facts("factorizeMeanField() should specify a mean-field factorization") do
    #               [T2]   [T3]
    #                |      | 
    #                v      v 
    # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

    g = initializeFactoringGraph()
    rf = RecognitionFactorization()
    factorizeMeanField()

    @fact rf.graph --> currentGraph()

    @fact length(rf.subgraphs) --> 5

    # Test edge_to_subgraph
    @fact length(rf.edge_to_subgraph) --> 6
    sg_t1_a1 = rf.edge_to_subgraph[eg(:t1_a1)]
    sg_a1_g1 = rf.edge_to_subgraph[eg(:a1_g1)]
    sg_t2_g1 = rf.edge_to_subgraph[eg(:t2_g1)]
    sg_g1_g2 = rf.edge_to_subgraph[eg(:g1_g2)]
    sg_t3_g2 = rf.edge_to_subgraph[eg(:t3_g2)]
    sg_g2_t4 = rf.edge_to_subgraph[eg(:g2_t4)]
    @fact is(sg_t1_a1, sg_a1_g1) --> true
    @fact is(sg_a1_g1, sg_t2_g1) --> false
    @fact is(sg_t2_g1, sg_g1_g2) --> false
    @fact is(sg_g1_g2, sg_t3_g2) --> false
    @fact is(sg_t3_g2, sg_g2_t4) --> false

    # Test node_subgraph_to_internal_edges
    @fact length(rf.node_subgraph_to_internal_edges) --> 6
    @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg_a1_g1)] --> Set([eg(:a1_g1)])
    @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg_t2_g1)] --> Set([eg(:t2_g1)])
    @fact rf.node_subgraph_to_internal_edges[(n(:g1), sg_g1_g2)] --> Set([eg(:g1_g2)])
    @fact rf.node_subgraph_to_internal_edges[(n(:g2), sg_g1_g2)] --> Set([eg(:g1_g2)])
    @fact rf.node_subgraph_to_internal_edges[(n(:g2), sg_t3_g2)] --> Set([eg(:t3_g2)])
    @fact rf.node_subgraph_to_internal_edges[(n(:g2), sg_g2_t4)] --> Set([eg(:g2_t4)])
end

facts("initialize() should specify a recognition distribution for a node-subgraph combination") do
    context("Mean field case") do
        #               [T2]   [T3]
        #     here       |      |
        #      v         v      v
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

        g = initializeFactoringGraph()
        rf = RecognitionFactorization()

        factor(eg(:t1_a1)) # Should automatically bind to current graph and create a RecognitionFactorization instance
        initialize(eg(:a1_g1), vague(Gaussian))

        @fact rf.initial_recognition_distributions --> Partitioned([vague(Gaussian)])
        @fact rf.recognition_distributions --> Partitioned([vague(Gaussian)])
        @fact is(rf.initial_recognition_distributions, rf.recognition_distributions) --> false
        
        sg1 = sg(:subgraph1)
        @fact rf.node_subgraph_to_recognition_distribution[(n(:g1), sg1)] --> rf.recognition_distributions.factors[1]
    end

    context("Structured case") do
        #               [T2]   [T3]
        #                |      | < here
        #                v      v   v
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

        g = initializeFactoringGraph()
        rf = RecognitionFactorization()

        factor(Set([eg(:t3_g2), eg(:g2_t4)])) # Should automatically bind to current graph and create a RecognitionFactorization instance
        initialize(Set([eg(:t3_g2), eg(:g2_t4)]), vague(MvGaussian{2}))

        @fact rf.initial_recognition_distributions --> Partitioned([vague(MvGaussian{2})])
        @fact rf.recognition_distributions --> Partitioned([vague(MvGaussian{2})])
        @fact is(rf.initial_recognition_distributions, rf.recognition_distributions) --> false
        
        sg1 = sg(:subgraph1)
        @fact rf.node_subgraph_to_recognition_distribution[(n(:g2), sg1)] --> rf.recognition_distributions.factors[1]
        @fact rf.node_subgraph_to_recognition_distribution[(n(:g2), sg1)] --> rf.initial_recognition_distributions.factors[1]
    end

    context("Invalid cases") do
        #               [T2]   [T3]
        #                |      |
        #                v      v   
        # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

        g = initializeFactoringGraph()
        rf = RecognitionFactorization()

        # Try to set distribution before partitioning
        @fact_throws initialize(eg(:a1_g1), vague(Gaussian))

        # Try to set distribution over edges on separate subgraphs 
        factor(eg(:a1_g1))
        factor(eg(:t2_g1))
        @fact_throws initialize(Set([eg(:a1_g1), eg(:t2_g1)]), vague(MvGaussian{2}))
    end
end

facts("verifyProperFactorization should return an error if a factorization is improper and return true otherwise") do
    #               [T2]   [T3]
    #                |      |
    #                v      v   
    # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

    g = initializeFactoringGraph()
    rf = RecognitionFactorization()

    factor(Set([eg(:a1_g1), eg(:t2_g1)]))

    # Partial factorization
    @fact_throws ForneyLab.verifyProper(rf)

    factor(Set([eg(:t3_g2), eg(:g2_t4)]))
    factor(eg(:g1_g2))

    initialize(Set([eg(:a1_g1), eg(:t2_g1)]), vague(MvGaussian{2}))

    # Partially set recognition distributions
    @fact_throws ForneyLab.verifyProper(rf)

    initialize(Set([eg(:t3_g2), eg(:g2_t4)]), vague(MvGaussian{2}))
    initialize(eg(:g1_g2), vague(Gaussian))

    # Finished recognition factorization specification
    @fact ForneyLab.verifyProper(rf) --> true
end


facts("resetRecognitionDistributions!() should reset recognition distributions to vague") do
    #               [T2]   [T3]
    #                |      |
    #                v      v   
    # [T1]-->[A1]-->[N1]-->[N2]-->[T4]

    g = initializeFactoringGraph()
    rf = RecognitionFactorization()

    factorizeMeanField()
    initialize(eg(:a1_g1), vague(Gaussian))
    initialize(eg(:t2_g1), vague(Gaussian))
    initialize(eg(:g1_g2), vague(Gaussian))
    initialize(eg(:t3_g2), vague(Gaussian))
    initialize(eg(:g2_t4), vague(Gaussian))
    @fact rf.initial_recognition_distributions --> Partitioned([vague(Gaussian), vague(Gaussian), vague(Gaussian), vague(Gaussian), vague(Gaussian)])

    # Change something so reset will have an effect
    rf.recognition_distributions.factors[1].V = 1.0
    rf.recognition_distributions.factors[2].V = 1.0
    rf.recognition_distributions.factors[3].V = 1.0
    rf.recognition_distributions.factors[4].V = 1.0
    rf.recognition_distributions.factors[5].V = 1.0
    @fact rf.recognition_distributions --> Partitioned([Gaussian(m=0.0,V=1.0), Gaussian(m=0.0,V=1.0), Gaussian(m=0.0,V=1.0), Gaussian(m=0.0,V=1.0), Gaussian(m=0.0,V=1.0)])

    ForneyLab.resetRecognitionDistributions!(rf)
    @fact rf.recognition_distributions --> Partitioned([vague(Gaussian), vague(Gaussian), vague(Gaussian), vague(Gaussian), vague(Gaussian)])
end
