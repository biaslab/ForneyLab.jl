#####################
# Unit tests
#####################

facts("Dirichlet unit tests") do
    context("construction") do
        @fact Dirichlet([2.0; 1.0]).alpha --> [2.0; 1.0]
        @fact typeof(Dirichlet([0.1; 0.8; 0.1])) --> Dirichlet{3}
    end

    context("required methods") do
        @fact vague(Dirichlet{4}) --> Dirichlet([tiny;tiny;tiny;tiny])
        @fact ForneyLab.vague!(Dirichlet([0.1; 0.9])) --> Dirichlet([tiny;tiny])
        @fact Dirichlet() --> Dirichlet()
        @fact (Dirichlet([1.0; 0.75]) == vague(Dirichlet{2})) --> false
        @fact isProper(Dirichlet()) --> true
        @fact isProper(Dirichlet([-0.2; 1.2])) --> false
    end

    context("prod!") do
        @fact Dirichlet([2.0; 2.0]) * Dirichlet([2.0; 2.0]) --> Dirichlet([3.0; 3.0])
        @fact_throws MethodError Dirichlet([0.5; 0.5]) * Dirichlet([1.0;1.0;1.0])
        @fact Dirichlet([1.0;1.0]) * MvDelta([0.5;0.5]) --> MvDelta([0.5;0.5])
        @fact MvDelta([0.5;0.5]) * Dirichlet([1.0;1.0]) --> MvDelta([0.5;0.5])
        @fact_throws DomainError MvDelta([0.5;1.5]) * Dirichlet([1.0;1.0])
        @fact_throws MethodError MvDelta([0.5;0.5;0.5]) * Dirichlet([1.0;1.0])
    end
end