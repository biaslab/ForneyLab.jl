#####################
# Unit tests
#####################

facts("StudentsTDistribution unit tests") do
    context("StudentsTDistribution() should initiatize a student's t-distribution") do
        dist = StudentsTDistribution()
        @fact dist.m => [0.0]
        @fact dist.W => reshape([1.0], 1, 1)
        @fact dist.nu => 1.0
        @fact mean(StudentsTDistribution(m=0.0, W=1.0, nu=3.0)) => [0.0]
        @fact var(StudentsTDistribution(m=0.0, W=1.0, nu=3.0)) => reshape([3.0], 1, 1)
        @fact isnan(mean(StudentsTDistribution(m=0.0, W=1.0, nu=0.5))[1]) => true
        @fact isnan(var(StudentsTDistribution(m=0.0, W=1.0, nu=1.0))[1,1]) => true
        @fact vague(StudentsTDistribution) => StudentsTDistribution(m=0.0, W=tiny(), nu=huge())
    end
end
