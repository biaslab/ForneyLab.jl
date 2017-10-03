@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Message{PointMass}, Void),
                :name          => SPBernoulliPV)

@variationalRule(:node_type     => Bernoulli,
                 :outbound_type => Message{Bernoulli},
                 :outbound_id   => 2,
                 :name          => VBBernoulli2)