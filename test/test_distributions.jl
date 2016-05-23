# Probability distribution test
facts("ProbabilityDistribution unit tests") do
    for probdist_type in subtypes(ProbabilityDistribution)
        probdist_type.abstract && continue # Skip abstract types

        context("$(probdist_type) should have a default constructor and a == operator") do
            @fact probdist_type()==probdist_type() --> true
        end
        context("dimensions() should be implemented for $(probdist_type)") do
            @fact applicable(dimensions, probdist_type) --> true
            @fact applicable(dimensions, probdist_type()) --> true
        end
        context("$(probdist_type) should provide a isProper method") do
            @fact isProper(probdist_type()) --> true
        end

        if !(probdist_type <: AbstractDelta)
            context("$(probdist_type) should provide a vague! method") do
                @fact applicable(ForneyLab.vague!, probdist_type()) --> true
            end
        end

    end

    # Univariate PDF can be multiplied with Delta
    for dtype in subtypes(Univariate)
        (dtype==Delta) && continue
        context("$(dtype) can be multiplied with a Delta") do
            @fact length(methods(ForneyLab.prod!, [dtype, Delta, Any])) --> greater_than_or_equal(1)
        end
    end

    # Multivariate PDF can be multiplied with MvDelta
    for dtype in subtypes(Multivariate)
        (dtype==MvDelta) && continue
        context("$(dtype) can be multiplied with a MvDelta") do
            @fact length(methods(ForneyLab.prod!, [dtype, MvDelta, Any])) --> greater_than_or_equal(1)
        end
    end

    # Matrixvariate PDF can be multiplied with MatrixDelta
    for dtype in subtypes(MatrixVariate)
        (dtype==MatrixDelta) && continue
        context("$(dtype) can be multiplied with a MatrixDelta") do
            @fact length(methods(ForneyLab.prod!, [dtype, MatrixDelta, Any])) --> greater_than_or_equal(1)
        end
    end
end
