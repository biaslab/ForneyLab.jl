#####################
# Unit tests
#####################

facts("Interface unit tests") do
    # Test setMessage!, clearMessage!, message, ensureMessage!, handle
    FactorGraph()
    n = MockNode()
    @fact_throws deepcopy(n.i[:out])
    @fact setMessage!(n.i[:out], Message(Gaussian(m=3.0, V=2.0))) --> Message(Gaussian(m=3.0, V=2.0))
    @fact typeof(n.i[:out].message) --> Message{Gaussian}
    @fact ForneyLab.ensureMessage!(n.i[:out], Gaussian) --> Message(Gaussian(m=3.0, V=2.0))
    @fact clearMessage!(n.i[:out]) --> nothing
    @fact n.i[:out].message --> nothing
    @fact ForneyLab.ensureMessage!(n.i[:out], Gaussian) --> Message(vague(Gaussian))
    @fact handle(n.interfaces[1]) --> "out"
    @fact handle(EqualityNode().interfaces[1]) --> "1"
end
