@sumProductRule(:node_type     => GainEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}, Message{PointMass}),
                :name          => SPGainEqualityOutVGGP)

@sumProductRule(:node_type     => GainEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGainEqualityIn1GVGP)
