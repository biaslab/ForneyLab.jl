facts("Approximations") do
    context("Parameters of Approximation should be subtype-constrained") do
        @fact Approximation{Gaussian, Laplace} --> Approximation{Gaussian, Laplace}
        @fact_throws Approximation{Gaussian, Int64}
        @fact_throws Approximation{Int64, Laplace}
    end
end