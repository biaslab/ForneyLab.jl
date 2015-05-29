facts("VMP.collectInbounds() tests") do
    context("collectInbounds() should add the proper message/marginal") do
        # Mean field factorized Gaussian node
        (node, edges) = initializeGaussianNode()
        algo = VMP.Algorithm()
        @fact VMP.collectInbounds(node.i[:mean]) => (1, [nothing, vague(GammaDistribution), vague(GaussianDistribution)])
        @fact VMP.collectInbounds(node.i[:precision]) => (2, [vague(GaussianDistribution), nothing, vague(GaussianDistribution)])
        @fact VMP.collectInbounds(node.i[:out]) => (3, [vague(GaussianDistribution), vague(GammaDistribution), nothing])

        # Structurally factorized
        (node, edges) = initializeGaussianNode()
        algo = VMP.Algorithm(Set{Edge}({node.i[:out].edge})) # Split off extensions of these groups into separate subgraphs
        @fact VMP.collectInbounds(node.i[:mean]) => (1, [nothing, Message(GammaDistribution()), vague(GaussianDistribution)])
        @fact VMP.collectInbounds(node.i[:precision]) => (2, [Message(GaussianDistribution()), nothing, vague(GaussianDistribution)])
        @fact VMP.collectInbounds(node.i[:out]) => (3, [vague(NormalGammaDistribution), vague(NormalGammaDistribution), nothing])
    end
end

# Test VMP specific functionality
include("test_q_distribution.jl")
include("test_q_factorization.jl")
include("test_subgraph.jl")
include("test_generate_schedule.jl")
include("test_message_passing.jl")
