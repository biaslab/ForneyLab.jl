@sumProductRule(:node_type     => UnitEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPUnitEqualityOutVGG)

@sumProductRule(:node_type     => UnitEquality,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}),
                :name          => SPUnitEqualityIn1GVG)
