#####################
# Unit tests
#####################

facts("MvLogNormalDistribution unit tests") do
    context("MvLogNormalDistribution() should initiatize a multi-variate log-normal distribution") do
        dist = MvLogNormalDistribution(m=[2.0], S=0.5*eye(1))
        @fact dist.m --> [2.0]
        @fact dist.S --> reshape([0.5],1,1)
        @fact mean(dist) --> roughly(exp(dist.m + 0.5*diag(dist.S)))
        @fact var(dist) --> roughly(exp(2.0*dist.m + diag(dist.S)).*(exp(diag(dist.S)) - 1.0))
        d = MvLogNormalDistribution(m=[0.0], S=-1.0*eye(1))
        @fact isnan(mean(d)[1]) --> true
    end

    context("vague() should initialize a vague (almost uninformative) log-normal distribution") do
        dist = vague(MvLogNormalDistribution{2})
        @fact dist.m --> [0.0, 0.0]
        @fact dist.S --> Diagonal([huge, huge])
    end

    context("prod!() should calculate product with a MvDeltaDistribution") do
        @fact MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)) * MvDeltaDistribution(ones(1)) --> MvDeltaDistribution(ones(1))
        @fact MvDeltaDistribution(ones(1)) * MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)) --> MvDeltaDistribution(ones(1))
        @fact_throws MethodError MvDeltaDistribution(ones(2)) * MvLogNormalDistribution(m=[2.0], S=0.5*eye(1))
        @fact_throws DomainError MvDeltaDistribution(-1.*ones(1)) * MvLogNormalDistribution(m=[2.0], S=0.5*eye(1))
        @fact ForneyLab.prod!(MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)), MvDeltaDistribution(ones(1)), MvLogNormalDistribution(m=[2.0], S=0.5*eye(1))) --> MvLogNormalDistribution(m=[0.0], S=tiny*eye(1))
        @fact_throws MethodError ForneyLab.prod!(MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)), MvDeltaDistribution(ones(2)), MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)))
        @fact_throws DomainError ForneyLab.prod!(MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)), MvDeltaDistribution(-1.*ones(1)), MvLogNormalDistribution(m=[2.0], S=0.5*eye(1)))
    end
end