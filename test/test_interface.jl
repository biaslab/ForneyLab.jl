#####################
# Unit tests
#####################

facts("Interface unit tests") do
    # Test setMessage!, clearMessage!, message, getOrCreateMessage!, name
    n = MockNode()
    @fact setMessage!(n.out, Message(GaussianDistribution(m=3.0, V=2.0))) => Message(GaussianDistribution(m=3.0, V=2.0))
    @fact typeof(n.out.message) => Message{GaussianDistribution}
    @fact message(n.out) => n.out.message
    @fact ForneyLab.getOrCreateMessage!(n.out, GaussianDistribution) => Message(GaussianDistribution(m=3.0, V=2.0))
    @fact clearMessage!(n.out) => nothing
    @fact message(n.out) => nothing
    @fact ForneyLab.getOrCreateMessage!(n.out, GaussianDistribution) => Message(GaussianDistribution())
    @fact name(n.interfaces[1]) => "out"
end
