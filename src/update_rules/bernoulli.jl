@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Void, Message{PointMass}),
                :name          => SPBernoulliOutVP)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Void, ProbabilityDistribution),
                      :name          => VBBernoulliOut)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Beta},
                      :inbound_types => (ProbabilityDistribution, Void),
                      :name          => VBBernoulliIn1)