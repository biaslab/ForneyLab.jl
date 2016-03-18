#####################
# Unit tests
#####################

facts("MatrixDeltaDistribution unit tests") do
    context("MatrixDeltaDistribution() should initialize a matrix delta distribution") do
        @fact MatrixDeltaDistribution().M --> [1.0].'
        @fact typeof(MatrixDeltaDistribution([2.0].')) --> MatrixDeltaDistribution{Float64, 1, 1}
        @fact MatrixDeltaDistribution([2.0].').M --> [2.0].'
        @fact typeof(MatrixDeltaDistribution(eye(2))) --> MatrixDeltaDistribution{Float64, 2, 2}
        @fact MatrixDeltaDistribution(eye(2)).M --> eye(2)
        @fact MatrixDeltaDistribution(diageye(2)).M --> diageye(2)
        @fact typeof(MatrixDeltaDistribution()) --> MatrixDeltaDistribution{Float64, 1, 1}
    end

    context("MatrixDeltaDistribution can be sampled") do
        @fact sample(MatrixDeltaDistribution([2.0].')) --> [2.0].'
    end

    context("There should be no such thing as vague(MatrixDeltaDistribution)") do
        @fact_throws vague(MatrixDeltaDistribution{Float64, 2, 2})
    end

    context("Matrices should convert to MatrixDeltaDistribution") do
        @fact convert(ProbabilityDistribution, eye(2)) --> MatrixDeltaDistribution(eye(2))
    end
end
