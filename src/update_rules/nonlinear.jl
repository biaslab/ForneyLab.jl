@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}),
                :name          => SPNonlinearOutVG)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void),
                :name          => SPNonlinearIn1GV)