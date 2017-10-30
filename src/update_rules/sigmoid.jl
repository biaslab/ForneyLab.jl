@sumProductRule(:node_type     => Sigmoid,
                :outbound_type => Message{Univariate{Bernoulli}},
                :inbound_types => (Void, Message{Univariate{Gaussian}}),
                :name          => SPSigmoidBinG)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{Univariate{Gaussian}},
                            :inbound_types => (Message{Univariate{Bernoulli}}, Message{Univariate{Gaussian}}),
                            :outbound_id   => 2,
                            :name          => EPSigmoidRealGB)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{Univariate{Gaussian}},
                            :inbound_types => (Message{Univariate{PointMass}}, Message{Univariate{Gaussian}}),
                            :outbound_id   => 2,
                            :name          => EPSigmoidRealGP)