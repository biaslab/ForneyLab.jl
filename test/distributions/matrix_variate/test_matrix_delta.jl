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
        @fact pdf(MatrixDelta(eye(3)), eye(3)) --> 1.0
        @fact pdf(MatrixDelta(eye(3)), eye(4)) --> 0.0
        @fact pdf(MatrixDelta(eye(3)), 2.0*eye(3)) --> 0.0
    end

    context("MatrixDelta can be sampled") do
        @fact sample(MatrixDelta([2.0].')) --> [2.0].'
    end

    context("Product of two MatrixDeltas") do
        @fact MatrixDelta([2.0 1.0;2.0 1.0]) * MatrixDelta([2.0 1.0;2.0 1.0]) --> MatrixDelta([2.0 1.0;2.0 1.0])
        @fact_throws MatrixDelta([1.0 2.0; 1.0 2.0]) * MatrixDelta([2.0 1.0; 1.0 1.0])
        @fact_throws MethodError MatrixDelta([2.0 5.0 6.0;2.0 4.0 1.0]) * MatrixDelta([2.0 1.0;2.0 1.0])
        @fact ForneyLab.prod!(MatrixDelta([2.0 1.0;2.0 1.0]), nothing) --> MatrixDelta([2.0 1.0;2.0 1.0])
    end

    context("unsafeDetLogMean() should return correct result") do
        @fact ForneyLab.unsafeDetLogMean(MatrixDelta(eye(2))) --> 0.0
        @fact ForneyLab.unsafeDetLogMean(MatrixDelta([1.0 0.0;0.0 2.0])) --> log(2)
    end

    context("Matrices should convert to MatrixDelta") do
        @fact convert(ProbabilityDistribution, eye(2)) --> MatrixDelta(eye(2))
    end
end
