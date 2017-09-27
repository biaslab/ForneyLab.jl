@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Message{PointMass}, Void),
                :name          => SPBernoulliPV)