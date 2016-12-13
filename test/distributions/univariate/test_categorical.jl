#####################
# Unit tests
#####################

facts("Categorical unit tests") do
    context("construction") do
        @fact Categorical().p --> [0.5; 0.5]
        @fact Categorical([0.1; 0.9]).p --> [0.1; 0.9]
        @fact typeof(Categorical([0.1; 0.8; 0.1])) --> Categorical{3}
        @fact pdf(Categorical([0.1; 0.8; 0.1]), 2) --> 0.8
        @fact pdf(Categorical([0.1; 0.8; 0.1]), 4) --> 0.0
    end

    context("required methods") do
        @fact vague(Categorical{4}) --> Categorical([0.25; 0.25; 0.25; 0.25])
        @fact mean(Categorical([0.1, 0.3, 0.6])) --> [0.1, 0.3, 0.6]
        @fact ForneyLab.vague!(Categorical([0.1; 0.9])) --> Categorical([0.5; 0.5])
        @fact Categorical() --> Categorical()
        @fact (Categorical([0.25; 0.75]) == vague(Categorical{2})) --> false
        @fact isProper(Categorical()) --> true
        @fact isProper(Categorical([0.2; 0.2])) --> false
    end

    context("prod!") do
        @fact Categorical([0.25; 0.75]) * Categorical([0.75; 0.25]) --> Categorical([0.5; 0.5])
        @fact_throws Categorical([0.0; 1.]) * Categorical([1.0; 0.])
        @fact_throws MethodError Categorical([0.5; 0.5]) * Categorical([0.25;0.25;0.5])
        @fact Categorical([0.25;0.25;0.5]) * Delta(1) --> Delta(1)
        @fact Delta(3) * Categorical([0.25;0.25;0.5]) --> Delta(3)
        @fact_throws DomainError Delta(4) * Categorical([0.25;0.25;0.5])
        @fact_throws MethodError Delta(2.0) * Categorical([0.25;0.25;0.5])
        @fact_throws Delta(2) * Categorical([0.25;0.0;0.75])
        @fact ForneyLab.prod!(Categorical([0.25; 0.75]), Delta(1), Categorical()) --> Categorical([1.0;0.0])
        @fact ForneyLab.prod!(Delta(2), Categorical([0.25; 0.75]), Categorical()) --> Categorical([0.0;1.0])
        @fact_throws DomainError ForneyLab.prod!(Delta(3), Categorical([0.25; 0.75]), Categorical())
    end
end