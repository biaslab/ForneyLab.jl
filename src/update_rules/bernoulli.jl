@sumProductRule(:node_type     => Bernoulli,
                :outbound_type => Message{Univariate{Bernoulli}},
                :inbound_types => (Void, Message{Univariate{PointMass}}),
                :name          => SPBernoulliOutP)

@variationalRule(:node_type     => Bernoulli,
                 :outbound_type => Message{Univariate{Bernoulli}},
                 :inbound_types => (Void, Univariate),
                 :name          => VBBernoulliOut)