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

    context("getOrCreateMessage should assign a message to an interface if there is none") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        @fact node.out.message => nothing
        @fact ForneyLab.getOrCreateMessage(node.out, GaussianDistribution) => Message(GaussianDistribution())
        @fact node.out.message => Message(GaussianDistribution())
    end

    context("getOrCreateMessage should return the present message is there is one and the type matches") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        node.out.message = Message(GaussianDistribution(m=5.0, V=1.0))
        @fact ForneyLab.getOrCreateMessage(node.out, GaussianDistribution) => Message(GaussianDistribution(m=5.0, V=1.0))
        @fact node.out.message => Message(GaussianDistribution(m=5.0, V=1.0))
    end

    context("getOrCreateMessage should assign a message to an interface if the type does not match") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        node.out.message = Message(GaussianDistribution(m=5.0, V=1.0))
        @fact ForneyLab.getOrCreateMessage(node.out, DeltaDistribution) => Message(DeltaDistribution())
        @fact node.out.message => Message(DeltaDistribution())
    end

    context("KLpq should calculate a numeric approximation to the KL divergence") do
        x = [0:0.1:2]
        p = vec([ones(11,1)-tiny(); zeros(10,1)+tiny()])
        q1 = vec([ones(11,1)-tiny(); zeros(10,1)+tiny()])
        q2 = vec([zeros(10,1)+tiny(); ones(11,1)-tiny()])
        @fact KLpq(x, p, q1) => 0.0
        @fact KLpq(x, p, q2) => 27.631021115875058
    end
end