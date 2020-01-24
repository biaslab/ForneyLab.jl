module RecognitionFactorTest

using Test
using ForneyLab

using ForneyLab: nodesConnectedToExternalEdges, Cluster, condense, flatten

@testset "RecognitionFactor" begin
    g = FactorGraph()
    @RV m ~ GaussianMeanVariance(constant(0.0), constant(1.0))
    @RV w ~ Gamma(constant(1.0), constant(1.0))
    y = Variable[]
    for i = 1:3
        @RV y_i ~ GaussianMeanPrecision(m, w)
        placeholder(y_i, :y, index=i)
        push!(y, y_i)
    end

    algo = Algorithm()
    q_m = RecognitionFactor(m)
    @test q_m.id == :recognitionfactor_1
    @test q_m.internal_edges == edges(m)
    @test algo.recognition_factors[:recognitionfactor_1] === q_m

    q_w = RecognitionFactor(w)
    @test q_w.id == :recognitionfactor_2
    @test q_w.internal_edges == edges(w)
    @test algo.recognition_factors[:recognitionfactor_2] === q_w

    # Joint factorizations
    q_m_w = RecognitionFactor([m, w])
    @test q_m_w.id == :recognitionfactor_3
    @test q_m_w.internal_edges == edges(Set([m, w]))
    @test algo.recognition_factors[:recognitionfactor_3] === q_m_w

    q_y = RecognitionFactor(y)
    @test q_y.id == :recognitionfactor_4
    @test q_y.internal_edges == edges(Set(y))
    @test algo.recognition_factors[:recognitionfactor_4] === q_y
end

@testset "setTargets!" begin
    @test true == false
end

end # module