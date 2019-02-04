@sumProductRule(:node_type     => GainEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPGainEqualityOutVGG)

@sumProductRule(:node_type     => GainEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}),
                :name          => SPGainEqualityIn1GVG)

@sumProductRule(:node_type     => GainEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing),
                :name          => SPGainEqualityIn2GGV)
