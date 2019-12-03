@sumProductRule(:node_type     => Probit,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPProbitOutNG)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Bernoulli}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1GB)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Categorical}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1GC)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{PointMass}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1GP)