module RecognitionFactorizationTest

using Test
using ForneyLab

@testset "RecognitionFactorization" begin
    rfz = RecognitionFactorization()
    @test rfz.recognition_factors == Dict{Symbol, RecognitionFactor}()
    @test rfz.edge_to_recognition_factor == Dict{Edge, RecognitionFactor}()
    @test rfz.node_edge_to_cluster == Dict{Tuple{FactorNode, Edge}, Cluster}()
    @test currentRecognitionFactorization() === rf
end

end # module