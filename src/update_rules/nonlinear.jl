@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{Gaussian}),
                :name          => SPNonlinearOutVG)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{Gaussian}, Void),
                :name          => SPNonlinearIn1GV)