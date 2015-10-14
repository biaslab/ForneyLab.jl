#####################
# Unit tests
#####################

facts("StudentsTDistribution unit tests") do
    context("StudentsTDistribution() should initiatize a student's t-distribution") do
        dist = StudentsTDistribution()
        @fact dist.m => 0.0
        @fact dist.lambda => 1.0
        @fact dist.nu => huge
        @fact mean(StudentsTDistribution(m=0.0, lambda=1.0, nu=3.0)) => 0.0
        @fact var(StudentsTDistribution(m=0.0, lambda=1.0, nu=3.0)) => 3.0
        @fact isnan(mean(StudentsTDistribution(m=0.0, lambda=1.0, nu=1.0))) => true
        @fact isnan(var(StudentsTDistribution(m=0.0, lambda=1.0, nu=0.5))) => true
        @fact vague(StudentsTDistribution) => StudentsTDistribution(m=0.0, lambda=tiny, nu=huge)
    end
end
