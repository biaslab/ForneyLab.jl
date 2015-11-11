# Probability distribution test
facts("General ProbabilityDistribution unit tests") do
    for probdist_type in [subtypes(UnivariateProbabilityDistribution); subtypes(MultivariateProbabilityDistribution)]
        context("$(probdist_type) should have a default constructor and a == operator") do
            @fact probdist_type()==probdist_type() --> true
        end
        context("$(probdist_type) should provide a isProper method") do
            @fact isProper(probdist_type()) --> true
        end
    end
end
