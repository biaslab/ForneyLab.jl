#####################
# Unit tests
#####################

facts("Helper function unit tests") do
    context("ensureMatrix should convert an array with one element to a matrix type") do
        @fact typeof(ForneyLab.ensureMatrix([1.0])) => Array{Float64, 2} # Cast 1D to 2D array
        @fact ForneyLab.ensureMatrix([1.0]) => reshape([1.0], 1, 1)
        @fact ForneyLab.ensureMatrix(eye(2)) => eye(2)
    end

    context("isApproxEqual should work for scalars, vectors and matrices") do
        @fact isApproxEqual(1.0, 1.0+1e-15) => true
        @fact isApproxEqual(1.0, 1.0+1e-9) => false
        @fact isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-15) => true
        @fact isApproxEqual([1.0, 1.0], [1.0, 1.0]+1e-9) => false
        @fact isApproxEqual(eye(3,3), eye(3,3)+1e-15) => true
        @fact isApproxEqual(eye(3,3), eye(3,3)+1e-9) => false
    end

    context("getOrCreateMessage should assign a message to an interface if there is none and otherwise set a standard message") do
        node = TerminalNode(GaussianDistribution())
        @fact node.out.message => nothing
        getOrCreateMessage(node.out, GaussianDistribution)
        @fact typeof(node.out.message) => Message{GaussianDistribution}
        node2 = TerminalNode(2.0)
        @fact node2.out.message => nothing
        getOrCreateMessage(node2.out, Float64)
        @fact node2.out.message.payload => 1.0
        ForneyLab.updateNodeMessage!(node2, 1, Float64, nothing)
        @fact getOrCreateMessage(node2.out, Float64).payload => 2.0
    end

    context("getDuplicatedIds() should return values of duplicate entries") do
        @fact getDuplicatedIds([1, 2, 3]) => []
        @fact getDuplicatedIds([1, 2, 2, 3, 6, 6, 2, 2]) => [2, 6]
    end
end