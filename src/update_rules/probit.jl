@sumProductRule(:node_type     => Probit,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPProbitOutNG)

@sumProductRule(:node_type     => Probit,
                :outbound_type => Message{Function},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPProbitIn1PN)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Bernoulli}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1BG)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{Categorical}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1CG)

@expectationPropagationRule(:node_type     => Probit,
                            :outbound_type => Message{GaussianWeightedMeanPrecision},
                            :inbound_types => (Message{PointMass}, Message{Gaussian}),
                            :outbound_id   => 2,
                            :name          => EPProbitIn1PG)