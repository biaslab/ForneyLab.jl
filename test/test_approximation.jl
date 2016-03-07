facts("Approximation related types") do
    context("Parameters of Approximation should be subtype-constrained") do
        @fact Approximation{GaussianDistribution, Laplace} --> Approximation{GaussianDistribution, Laplace}
        @fact_throws Approximation{GaussianDistribution, Int64}
        @fact_throws Approximation{Int64, Laplace}
    end
end