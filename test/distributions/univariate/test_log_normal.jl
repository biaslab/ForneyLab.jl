#####################
# Unit tests
#####################

facts("LogNormal unit tests") do
    context("LogNormal() should initiatize a log-normal distribution") do
        dist = LogNormal(m=2.0, s=0.5)
        @fact dist.m --> 2.0
        @fact dist.s --> 0.5
        @fact mean(dist) --> exp(dist.m + 0.5*dist.s)
        @fact var(dist) --> (exp(dist.s) - 1)*exp(2*dist.m + dist.s)
        @fact isProper(LogNormal(m=0.0, s=-1.0)) --> false
        @fact pdf(LogNormal(m=2.0, s=0.5^2), 1.0) --> roughly(0.00026766045152977074, atol=1e-8)
    end

    context("vague() should initialize a vague (almost uninformative) log-normal distribution") do
        dist = vague(LogNormal)
        @fact dist.m --> 0.0
        @fact dist.s --> huge
    end

    context("prod!() should yield correct result") do
        @fact Delta(0.5) * LogNormal(m=2.0, s=0.5) --> Delta(0.5)
        @fact_throws DomainError LogNormal(m=2.0, s=0.5) * Delta(-1.1)
        @fact typeof(ForneyLab.prod!(Delta(0.5), LogNormal(m=2.0, s=0.5), LogNormal())) --> LogNormal
        @fact ForneyLab.unsafeMean(ForneyLab.prod!(LogNormal(m=2.0, s=0.5), Delta(0.5), LogNormal())) --> roughly(0.5)
        @fact ForneyLab.unsafeVar(ForneyLab.prod!(LogNormal(m=2.0, s=0.5), Delta(0.5), LogNormal())) --> less_than(1e-6)
    end
end