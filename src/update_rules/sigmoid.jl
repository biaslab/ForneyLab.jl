@sumProductRule(:node_type     => Sigmoid,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPSigmoidBinVG)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Bernoulli}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPSigmoidRealGB)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Categorical}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPSigmoidRealGC)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{PointMass}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPSigmoidRealGP)