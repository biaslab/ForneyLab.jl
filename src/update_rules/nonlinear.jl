@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPNonlinearOutVGP)
