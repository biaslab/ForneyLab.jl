@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Void, Message{PointMass}),
                :name          => SPBernoulliOutP)

@variationalRule(:node_type     => Bernoulli,
                 :outbound_type => Message{Bernoulli},
                 :inbound_types => (Void, ProbabilityDistribution),
                 :name          => VBBernoulliOut)