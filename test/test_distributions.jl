# Probability distribution test
facts("ProbabilityDistribution unit tests") do
    for probdist_type in [subtypes(UnivariateProbabilityDistribution); subtypes(MultivariateProbabilityDistribution); subtypes(MatrixvariateProbabilityDistribution)]
        context("$(probdist_type) should have a default constructor and a == operator") do
            @fact probdist_type()==probdist_type() --> true
        end
        context("$(probdist_type) should provide a isProper method") do
            @fact isProper(probdist_type()) --> true
        end
    end

    # Univarite PDF can be multiplied with DeltaDistribution
    for dtype in subtypes(UnivariateProbabilityDistribution)
        (dtype==DeltaDistribution) && continue
        context("$(dtype) can be multiplied with a DeltaDistribution") do
            @fact length(methods(ForneyLab.prod!, [dtype, DeltaDistribution, Any])) --> greater_than_or_equal(1)
        end
    end

    # Multivariate PDF can be multiplied with MvDeltaDistribution
    for dtype in subtypes(MultivariateProbabilityDistribution)
        (dtype==MvDeltaDistribution) && continue
        context("$(dtype) can be multiplied with a MvDeltaDistribution") do
            @fact length(methods(ForneyLab.prod!, [dtype, MvDeltaDistribution, Any])) --> greater_than_or_equal(1)
        end
    end

    # Matrixvariate PDF can be multiplied with MatrixDeltaDistribution
    for dtype in subtypes(MatrixvariateProbabilityDistribution)
        (dtype==MatrixDeltaDistribution) && continue
        context("$(dtype) can be multiplied with a MatrixDeltaDistribution") do
            @fact length(methods(ForneyLab.prod!, [dtype, MatrixDeltaDistribution, Any])) --> greater_than_or_equal(1)
        end
    end
end
