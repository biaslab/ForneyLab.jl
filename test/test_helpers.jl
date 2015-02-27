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

    context("isRoundedPosDef should check positive difiniteness with robustness for numerical errors") do
        K = inv([4.0 3.0 2.0;
                 3.0 4.0 3.0;
                 2.0 3.0 4.0])
        @fact isposdef(K) => false
        @fact ForneyLab.isRoundedPosDef(K) => true
    end

    context("ensureMessage! should assign a message to an interface if there is none") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        @fact node.out.message => nothing
        @fact ForneyLab.ensureMessage!(node.out, GaussianDistribution) => Message(vague(GaussianDistribution))
        @fact node.out.message => Message(vague(GaussianDistribution))
    end

    context("ensureMessage! should return the present message is there is one and the type matches") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        node.out.message = Message(GaussianDistribution(m=5.0, V=1.0))
        @fact ForneyLab.ensureMessage!(node.out, GaussianDistribution) => Message(GaussianDistribution(m=5.0, V=1.0))
        @fact node.out.message => Message(GaussianDistribution(m=5.0, V=1.0))
    end

    context("ensureMessage! should assign a message to an interface if the type does not match") do
        node = TerminalNode(GaussianDistribution(m=5.0, V=1.0))
        node.out.message = Message(GaussianDistribution(m=5.0, V=1.0))
        @fact ForneyLab.ensureMessage!(node.out, DeltaDistribution) => Message(DeltaDistribution())
        @fact node.out.message => Message(DeltaDistribution())
    end

    context("truncate() should truncate a string to a specified length") do
        @fact ForneyLab.truncate("spsbrats", 9) => "spsbrats"
        @fact ForneyLab.truncate("spsbrats", 7) => "spsb..."
    end

    context("pad() should pad a string with spaces to a specified length") do
        @fact ForneyLab.pad("spsbrats", 9) => "spsbrats "
        @fact ForneyLab.pad("spsbrats", 7) => "spsb..."
    end
end