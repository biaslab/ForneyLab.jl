#####################
# Unit tests
#####################

facts("GaussianMessage unit tests") do
    context("GaussianMessage() should initialize a Gaussian message") do
        @fact GaussianMessage().V => ones(1, 1)
        @fact GaussianMessage().m => [0.0]
        @fact GaussianMessage(m=[0.0], W=[1.0]).W => ones(1, 1)
        @fact GaussianMessage(xi=[0.0], W=[1.0]).W => ones(1, 1)
        @fact typeof(GaussianMessage(m=[1.0], V=[1.0]).V) => Array{Float64, 2} # cast single value to matrix
        @fact_throws GaussianMessage(V=[1.0], W=[1.0])
        @fact_throws GaussianMessage(m=[0.0], xi=[0.0])
        @fact_throws GaussianMessage(xi=[0.0])
    end
    context("Underdetermined GaussianMessage should be detected by isWellDefined()") do
        @fact isWellDefined(GaussianMessage()) => true
        @fact isWellDefined(GaussianMessage(m=[0.0], V=[1.0])) => true
        @fact isWellDefined(GaussianMessage(m=[0.0], W=[1.0])) => true
        @fact isWellDefined(GaussianMessage(xi=[0.0], V=[1.0])) => true
        @fact isWellDefined(GaussianMessage(xi=[0.0], W=[1.0])) => true
        @fact isWellDefined(GaussianMessage(m=[0.0], xi=[0.0], W=[1.0], V=[1.0])) => true

        msg = GaussianMessage(m=[0.0], V=[1.0])
        msg.m = nothing
        @fact isWellDefined(msg) => false

        msg = GaussianMessage(m=[0.0], V=[1.0])
        msg.V = nothing
        @fact isWellDefined(msg) => false

        msg = GaussianMessage(xi=[0.0], W=[1.0])
        msg.xi = nothing
        @fact isWellDefined(msg) => false

        msg = GaussianMessage(xi=[0.0], W=[1.0])
        msg.W = nothing
        @fact isWellDefined(msg) => false

        msg = GaussianMessage(m=[0.0], V=[1.0], W=[1.0])
        msg.m = nothing
        @fact isWellDefined(msg) => false

        msg = GaussianMessage(m=[0.0], xi=[0.0], V=[1.0])
        msg.V = nothing
        @fact isWellDefined(msg) => false
    end
    context("Conversions between valid parametrizations of a GaussianMessage should be consistent") do
        # Defined as (m,V)
        @fact isConsistent(ensureMWParametrization!(GaussianMessage(m=[0.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianMessage(m=[0.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianMessage(m=[0.0], V=[1.0]))) => true
        # Defined as (m,W)
        @fact isConsistent(ensureMVParametrization!(GaussianMessage(m=[0.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianMessage(m=[0.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianMessage(m=[0.0], W=[1.0]))) => true
        # Defined as (xi,V)
        @fact isConsistent(ensureMVParametrization!(GaussianMessage(xi=[2.0], V=[1.0]))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianMessage(xi=[2.0], V=[1.0]))) => true
        @fact isConsistent(ensureXiWParametrization!(GaussianMessage(xi=[2.0], V=[1.0]))) => true
        # Defined as (xi,W)
        @fact isConsistent(ensureMVParametrization!(GaussianMessage(xi=[2.0], W=[1.0]))) => true
        @fact isConsistent(ensureMWParametrization!(GaussianMessage(xi=[2.0], W=[1.0]))) => true
        @fact isConsistent(ensureXiVParametrization!(GaussianMessage(xi=[2.0], W=[1.0]))) => true
    end
    context("Inconsistent overdetermined GaussianMessage should be detected by isConsistent()") do
        @fact isConsistent(GaussianMessage(m=[0.0], xi=[1.0], W=[1.0])) => false
        @fact isConsistent(GaussianMessage(m=[0.0], V=[1.0], W=[2.0])) => false
        @fact isConsistent(GaussianMessage(m=[0.0], xi=[1.0], V=[1.0], W=[2.0])) => false
    end
end

facts("GeneralMessage unit tests") do
    context("GeneralMessage() should initiatize a message value") do
        @fact GeneralMessage(1.0).value => 1.0
        @fact GeneralMessage([1.0, 2.0]).value => [1.0, 2.0]
    end
end

facts("GammaMessage unit tests") do
    context("GammaMessage() should initiatize a message value") do
        msg = GammaMessage(a=2.0, b=0.5, inverted=true)
        @fact msg.inverted => true
        @fact msg.a => 2.0
        @fact msg.b => 0.5
    end
end