@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPBernoulliOutNP)

@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Beta},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPBernoulliIn1PN)

@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{Beta}),
                :name          => SPBernoulliOutNB)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Nothing, Distribution),
                      :name          => VBBernoulliOut)

@naiveVariationalRule(:node_type     => Bernoulli,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Distribution, Nothing),
                      :name          => VBBernoulliIn1)
