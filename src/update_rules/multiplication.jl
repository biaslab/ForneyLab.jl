@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Void),
                :name          => SPMultiplicationGPV)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPMultiplicationVPG)