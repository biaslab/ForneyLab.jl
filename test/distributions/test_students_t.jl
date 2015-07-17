#####################
# Unit tests
#####################

facts("StudentsTDistribution unit tests") do
    context("StudentsTDistribution() should initiatize a student's t-distribution") do
        dist = StudentsTDistribution()
        @fact dist.m => [0.0]
        @fact dist.lambda => reshape([1.0], 1, 1)
        @fact dist.nu => 1.0
        @fact mean(StudentsTDistribution(0.0, 1.0, 3.0)) => [0.0]
        @fact var(StudentsTDistribution(0.0, 1.0, 3.0)) => reshape([3.0], 1, 1)
        @fact isnan(mean(StudentsTDistribution(0.0, 1.0, 0.5))[1]) => true
        @fact isnan(var(StudentsTDistribution(0.0, 1.0, 1.0))[1,1]) => true
        @fact vague(StudentsTDistribution) => StudentsTDistribution(0.0, tiny(), huge())
    end
end
