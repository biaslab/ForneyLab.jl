#####################
# Unit tests
#####################

facts("GainNode unit tests") do
    context("GainNode(gain=A) should initialize a GainNode with 2 interfaces") do
        FactorGraph()
        GainNode(gain=[1.0], id=:node)
        @fact length(n(:node).interfaces) --> 2
        @fact n(:node).i[:in] --> n(:node).interfaces[1]
        @fact n(:node).i[:out] --> n(:node).interfaces[2]
        @fact typeof(n(:node).gain) <: AbstractArray --> true
        @fact length(size(n(:node).gain)) --> 2 # A should always be a matrix
    end

    context("GainNode() should initialize a GainNode with 3 interfaces") do
        FactorGraph()
        GainNode(id=:node)
        @fact length(n(:node).interfaces) --> 3
        @fact n(:node).i[:in] --> n(:node).interfaces[1]
        @fact n(:node).i[:out] --> n(:node).interfaces[2]
        @fact n(:node).i[:gain] --> n(:node).interfaces[3]
    end

    context("GainNode should provide sumProductRule! for Delta{Float64}") do
        # Backward message
        validateOutboundMessage(GainNode(gain=2.0),
                                1,
                                [nothing, Message(Delta(3.0))],
                                Delta(1.5))
        validateOutboundMessage(GainNode(),
                                1,
                                [nothing, Message(Delta(3.0)), Message(Delta(2.0))],
                                Delta(1.5))
        # Forward message
        validateOutboundMessage(GainNode(gain=2.0),
                                2,
                                [Message(Delta(3.0)), nothing],
                                Delta(6.0))
        validateOutboundMessage(GainNode(),
                                2,
                                [Message(Delta(3.0)), nothing, Message(Delta(2.0))],
                                Delta(6.0))
    end

    context("GainNode should provide sumProductRule! for MvDelta{Float64}") do
        # Backward message
        A = [1.0 0.5; -0.5 2.0]
        validateOutboundMessage(GainNode(gain=A),
                                1,
                                [nothing, Message(MvDelta([30.0, 10.0]))],
                                MvDelta(pinv(A)*[30.0, 10.0]))
        validateOutboundMessage(GainNode(),
                                1,
                                [nothing, Message(MvDelta([30.0, 10.0])), Message(MatrixDelta(A))],
                                MvDelta(pinv(A)*[30.0, 10.0]))
        # Forward message
        validateOutboundMessage(GainNode(gain=A),
                                2,
                                [Message(MvDelta([30.0, 10.0])), nothing],
                                MvDelta(A*[30.0, 10.0]))
        validateOutboundMessage(GainNode(),
                                2,
                                [Message(MvDelta([30.0, 10.0])), nothing, Message(MatrixDelta(A))],
                                MvDelta(A*[30.0, 10.0]))
    end

    context("GainNode should provide sumProductRule! for Gaussian") do
        context("(m,V) parametrization") do
            # Backward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, V=5.0))],
                                    Gaussian(m=1.5, V=1.25))
            validateOutboundMessage(GainNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, V=5.0)), Message(Delta(2.0))],
                                    Gaussian(m=1.5, V=1.25))
            # Forward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    2,
                                    [Message(Gaussian(m=3.0, V=5.0)), nothing],
                                    Gaussian(m=6.0, V=20.0))
            validateOutboundMessage(GainNode(),
                                    2,
                                    [Message(Gaussian(m=3.0, V=5.0)), nothing, Message(Delta(2.0))],
                                    Gaussian(m=6.0, V=20.0))
        end
        context("(m,W) parametrization") do
            # Backward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, W=2.0))],
                                    Gaussian(m=1.5, W=8.0))
            validateOutboundMessage(GainNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, W=2.0)), Message(Delta(2.0))],
                                    Gaussian(m=1.5, W=8.0))
            # Forward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    2,
                                    [Message(Gaussian(m=3.0, W=2.0)), nothing],
                                    Gaussian(m=6.0, W=0.5))
            validateOutboundMessage(GainNode(),
                                    2,
                                    [Message(Gaussian(m=3.0, W=2.0)), nothing, Message(Delta(2.0))],
                                    Gaussian(m=6.0, W=0.5))
        end
        context("(xi,W) parametrization") do
            # Backward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    1,
                                    [nothing, Message(Gaussian(xi=6.0, W=2.0))],
                                    Gaussian(xi=12.0, W=8.0))
            validateOutboundMessage(GainNode(),
                                    1,
                                    [nothing, Message(Gaussian(xi=6.0, W=2.0)), Message(Delta(2.0))],
                                    Gaussian(xi=12.0, W=8.0))
            # Forward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    2,
                                    [Message(Gaussian(xi=6.0, W=2.0)), nothing],
                                    Gaussian(xi=3.0, W=0.5))
            validateOutboundMessage(GainNode(),
                                    2,
                                    [Message(Gaussian(xi=6.0, W=2.0)), nothing, Message(Delta(2.0))],
                                    Gaussian(xi=3.0, W=0.5))
        end
        context("Improper distributions") do
            # Backward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, V=-5.0))],
                                    Gaussian(m=1.5, V=-1.25))
            validateOutboundMessage(GainNode(),
                                    1,
                                    [nothing, Message(Gaussian(m=3.0, V=-5.0)), Message(Delta(2.0))],
                                    Gaussian(m=1.5, V=-1.25))
            # Forward message
            validateOutboundMessage(GainNode(gain=2.0),
                                    2,
                                    [Message(Gaussian(m=3.0, V=-5.0)), nothing],
                                    Gaussian(m=6.0, V=-20.0))
            validateOutboundMessage(GainNode(),
                                    2,
                                    [Message(Gaussian(m=3.0, V=-5.0)), nothing, Message(Delta(2.0))],
                                    Gaussian(m=6.0, V=-20.0))
        end
    end

    context("GainNode should provide sumProductRule! for MvGaussian") do
        # The following tests on the update rules correspond to nodes 3 and 4 from Table 4.1 in:
        # Korl, Sascha. “A Factor Graph Approach to Signal Modelling, System Identification and Filtering.” Hartung-Gorre, 2005.
        A = [   3.0 2.0 1.0;
                2.0 3.0 2.0;
                1.0 2.0 3.0]
        context("(m,V) parametrization") do
            mean = collect(1.0:3.0)
            variance = [4.0 3.0 2.0;
                        3.0 4.0 3.0;
                        2.0 3.0 4.0]
            # Backward message
            validateOutboundMessage(GainNode(gain=A),
                                    1,
                                    [nothing, Message(MvGaussian(m=mean, V=variance))],
                                    MvGaussian(m=inv(A) * mean, V=inv(A) * variance * inv(A)'))
            validateOutboundMessage(GainNode(),
                                    1,
                                    [nothing, Message(MvGaussian(m=mean, V=variance)), Message(MatrixDelta(A))],
                                    MvGaussian(m=inv(A) * mean, V=inv(A) * variance * inv(A)'))
            # Forward message
            validateOutboundMessage(GainNode(gain=A),
                                    2,
                                    [Message(MvGaussian(m=mean, V=variance)), nothing],
                                    MvGaussian(m=A * mean, V=A * variance * A'))
            validateOutboundMessage(GainNode(),
                                    2,
                                    [Message(MvGaussian(m=mean, V=variance)), nothing, Message(MatrixDelta(A))],
                                    MvGaussian(m=A * mean, V=A * variance * A'))
        end

        context("(m,W) parametrization") do
            mean = collect(1.0:3.0)
            precision = inv([   4.0 3.0 2.0;
                                3.0 4.0 3.0;
                                2.0 3.0 4.0])
            # Backward message
            dist = ForneyLab.sumProductRule!(
                                        GainNode(gain=A),
                                        Val{1},
                                        MvGaussian(m=zeros(3), V=eye(3)),
                                        nothing, Message(MvGaussian(m=mean, W=precision)))
            ensureParameters!(dist, (:m, :W))
            @fact dist.m --> roughly(inv(A) * mean)
            @fact dist.W --> roughly(A' * precision * A)

            dist = ForneyLab.sumProductRule!(
                                        GainNode(),
                                        Val{1},
                                        MvGaussian(m=zeros(3), V=eye(3)),
                                        nothing, Message(MvGaussian(m=mean, W=precision)), Message(MatrixDelta(A)))
            ensureParameters!(dist, (:m, :W))
            @fact dist.m --> roughly(inv(A) * mean)
            @fact dist.W --> roughly(A' * precision * A)

            # Forward message
            dist = ForneyLab.sumProductRule!(
                                        GainNode(gain=A),
                                        Val{2},
                                        MvGaussian(m=zeros(3), V=eye(3)),
                                        Message(MvGaussian(m=mean, W=precision)), nothing)
            ensureParameters!(dist, (:m, :W))
            @fact dist.m --> roughly(A * mean)
            @fact dist.W --> roughly(inv(A)' * precision * inv(A))

            dist = ForneyLab.sumProductRule!(
                                        GainNode(),
                                        Val{2},
                                        MvGaussian(m=zeros(3), V=eye(3)),
                                        Message(MvGaussian(m=mean, W=precision)), nothing, Message(MatrixDelta(A)))
            ensureParameters!(dist, (:m, :W))
            @fact dist.m --> roughly(A * mean)
            @fact dist.W --> roughly(inv(A)' * precision * inv(A))
        end
    end
end
