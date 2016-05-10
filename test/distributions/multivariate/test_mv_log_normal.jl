#####################
# Unit tests
#####################

facts("MvLogNormal unit tests") do
    context("MvLogNormal() should initiatize a multi-variate log-normal distribution") do
        dist = MvLogNormal(m=[2.0], S=0.5*eye(1))
        @fact dist.m --> [2.0]
        @fact dist.S --> reshape([0.5],1,1)
        @fact mean(dist) --> roughly(exp(dist.m + 0.5*diag(dist.S)))
        @fact var(dist) --> roughly(exp(2.0*dist.m + diag(dist.S)).*(exp(diag(dist.S)) - 1.0))
        d = MvLogNormal(m=[0.0], S=-1.0*eye(1))
        @fact isnan(mean(d)[1]) --> true
    end

    context("vague() should initialize a vague (almost uninformative) log-normal distribution") do
        dist = vague(MvLogNormal{2})
        @fact dist.m --> [0.0, 0.0]
        @fact dist.S --> Diagonal([huge, huge])
    end

    context("prod!() should calculate product with a MvDelta") do
        @fact MvLogNormal(m=[2.0], S=0.5*eye(1)) * MvDelta(ones(1)) --> MvDelta(ones(1))
        @fact MvDelta(ones(1)) * MvLogNormal(m=[2.0], S=0.5*eye(1)) --> MvDelta(ones(1))
        @fact_throws MethodError MvDelta(ones(2)) * MvLogNormal(m=[2.0], S=0.5*eye(1))
        @fact_throws DomainError MvDelta(-1.*ones(1)) * MvLogNormal(m=[2.0], S=0.5*eye(1))
        @fact ForneyLab.prod!(MvLogNormal(m=[2.0], S=0.5*eye(1)), MvDelta(ones(1)), MvLogNormal(m=[2.0], S=0.5*eye(1))) --> MvLogNormal(m=[0.0], S=tiny*eye(1))
        @fact_throws MethodError ForneyLab.prod!(MvLogNormal(m=[2.0], S=0.5*eye(1)), MvDelta(ones(2)), MvLogNormal(m=[2.0], S=0.5*eye(1)))
        @fact_throws DomainError ForneyLab.prod!(MvLogNormal(m=[2.0], S=0.5*eye(1)), MvDelta(-1.*ones(1)), MvLogNormal(m=[2.0], S=0.5*eye(1)))
    end
end