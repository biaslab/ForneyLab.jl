#####################
# Unit tests
#####################

facts("Interface unit tests") do
    # Test setMessage!, clearMessage!, message, ensureMessage!, handle
    n = MockNode()
    @fact_throws deepcopy(n.i[:out])
    @fact setMessage!(n.i[:out], Message(GaussianDistribution(m=3.0, V=2.0))) --> Message(GaussianDistribution(m=3.0, V=2.0))
    @fact typeof(n.i[:out].message) --> Message{GaussianDistribution}
    @fact message(n.i[:out]) --> n.i[:out].message
    @fact ForneyLab.ensureMessage!(n.i[:out], GaussianDistribution) --> Message(GaussianDistribution(m=3.0, V=2.0))
    @fact clearMessage!(n.i[:out]) --> nothing
    @fact message(n.i[:out]) --> nothing
    @fact ForneyLab.ensureMessage!(n.i[:out], GaussianDistribution) --> Message(vague(GaussianDistribution))
    @fact handle(n.interfaces[1]) --> :out
    @fact handle(EqualityNode().interfaces[1]) --> 1
end
