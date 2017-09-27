@sumProductRule(:node_type     => Sigmoid,
                :outbound_type => Message{Bernoulli},
                :inbound_types => (Message{Gaussian}, Void),
                :name          => SPSigmoidGV)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{Gaussian},
                            :inbound_types => (Message{Gaussian}, Message{Bernoulli}),
                            :outbound_id   => 1,
                            :name          => EPSigmoidGB1)

@expectationPropagationRule(:node_type     => Sigmoid,
                            :outbound_type => Message{Gaussian},
                            :inbound_types => (Message{Gaussian}, Message{PointMass}),
                            :outbound_id   => 1,
                            :name          => EPSigmoidGP1)