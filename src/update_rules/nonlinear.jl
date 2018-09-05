@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearOutVG)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearIn1GV)