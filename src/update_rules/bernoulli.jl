@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPBernoulliOutVP)

@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Beta},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPBernoulliIn1PV)

@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{Beta}),
                :name          => SPBernoulliOutVB)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBBernoulliOut)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Beta},
                      :inbound_types => (ProbabilityDistribution, Nothing),
                      :name          => VBBernoulliIn1)
