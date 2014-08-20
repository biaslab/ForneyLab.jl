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

    context("uninformative() should initialize an uninformative student's t-distribution") do
        dist = uninformative(StudentsTDistribution)
        @fact dist.m => [0.0]
        @fact dist.W => reshape([0.001], 1, 1)
        @fact dist.nu => 1000.0
    end
end
