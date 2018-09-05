@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPBernoulliOutVP)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBBernoulliOut)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Beta},
                      :inbound_types => (ProbabilityDistribution, Nothing),
                      :name          => VBBernoulliIn1)