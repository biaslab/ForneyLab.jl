module AlgorithmTest

using Test
using ForneyLab

@testset "Algorithm" begin
    algo = Algorithm()
    @test algo.recognition_factors == Dict{Symbol, RecognitionFactor}()
    @test algo.edge_to_recognition_factor == Dict{Edge, RecognitionFactor}()
    @test currentAlgorithm() === algo
end

end # module