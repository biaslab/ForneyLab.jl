#####################
# Unit tests
#####################

facts("StudentsTDistribution unit tests") do
    context("StudentsTDistribution() should initiatize a student's t-distribution") do
        dist = StudentsTDistribution()
        @fact dist.m => [0.0]
        @fact dist.W => reshape([1.0], 1, 1)
        @fact dist.nu => 1.0
    end
end
