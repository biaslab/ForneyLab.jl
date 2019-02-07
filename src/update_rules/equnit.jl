@sumProductRule(:node_type     => EqUnit,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPEqUnitOutVGG)

@sumProductRule(:node_type     => EqUnit,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{Gaussian}),
                :name          => SPEqUnitIn1GVG)
