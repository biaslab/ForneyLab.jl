#####################
# Unit tests
#####################

facts("StudentsT unit tests") do
    context("StudentsT() should initiatize a student's t-distribution") do
        dist = StudentsT()
        @fact dist.m --> 0.0
        @fact dist.lambda --> 1.0
        @fact dist.nu --> huge
        @fact mean(StudentsT(m=0.0, lambda=1.0, nu=3.0)) --> 0.0
        @fact var(StudentsT(m=0.0, lambda=1.0, nu=3.0)) --> 3.0
        @fact isProper(StudentsT(m=0.0, lambda=1.0, nu=1.0)) --> false
        @fact isProper(StudentsT(m=0.0, lambda=1.0, nu=0.5)) --> false
        @fact vague(StudentsT) --> StudentsT(m=0.0, lambda=tiny, nu=4.0)
    end

    context("prod!() should yield correct result") do
        @fact Delta(0.5) * StudentsT(m=1.0, lambda=2.0, nu=4.0) --> Delta(0.5)
        @fact typeof(ForneyLab.prod!(Delta(0.5), StudentsT(m=1.0, lambda=2.0, nu=4.0), StudentsT())) --> StudentsT
        @fact mean(ForneyLab.prod!(StudentsT(m=1.0, lambda=2.0, nu=4.0), Delta(0.5), StudentsT())) --> roughly(0.5)
        @fact var(ForneyLab.prod!(StudentsT(m=1.0, lambda=2.0, nu=4.0), Delta(0.5), StudentsT())) --> less_than(1e-6)
    end

    context("Product of Gaussian and student's t-distribution should yield moment-matched Gaussian") do
        validateOutboundMessage(EqualityNode(),
                                3,
                                [Message(Gaussian()), Message(StudentsT(m=1.0, lambda=2.0, nu=4.0)), nothing],
                                Gaussian(m=0.5, W=2.0),
                                ForneyLab.sumProductRule!,
                                MomentMatching)
    end
end
