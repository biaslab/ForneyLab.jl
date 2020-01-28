module AlgorithmTest

using Test
using ForneyLab

@testset "Algorithm" begin
    algo = Algorithm()
    @test algo.recognition_factorization == currentRecognitionFactorization()
    @test currentAlgorithm() === algo
end

end # module