module InferenceAlgorithmTest

using Test
using ForneyLab

@testset "InferenceAlgorithm" begin
    algo = InferenceAlgorithm()
    @test algo.posterior_factorization == currentPosteriorFactorization()
    @test currentInferenceAlgorithm() === algo
end

end # module