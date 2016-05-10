#####################
# Unit tests
#####################

facts("MatrixDelta unit tests") do
    context("MatrixDelta() should initialize a matrix delta distribution") do
        @fact MatrixDelta().M --> [1.0].'
        @fact typeof(MatrixDelta([2.0].')) --> MatrixDelta{Float64, 1, 1}
        @fact MatrixDelta([2.0].').M --> [2.0].'
        @fact typeof(MatrixDelta(eye(2))) --> MatrixDelta{Float64, 2, 2}
        @fact MatrixDelta(eye(2)).M --> eye(2)
        @fact MatrixDelta(diageye(2)).M --> diageye(2)
        @fact typeof(MatrixDelta()) --> MatrixDelta{Float64, 1, 1}
        @fact dimensions(MatrixDelta{Float64,2,2}) --> (2, 2)
        @fact dimensions(MatrixDelta(eye(2))) --> (2, 2)
    end

    context("MatrixDelta can be sampled") do
        @fact sample(MatrixDelta([2.0].')) --> [2.0].'
    end

    context("There should be no such thing as vague(MatrixDelta)") do
        @fact_throws vague(MatrixDelta{Float64, 2, 2})
    end

    context("Matrices should convert to MatrixDelta") do
        @fact convert(ProbabilityDistribution, eye(2)) --> MatrixDelta(eye(2))
    end
end
