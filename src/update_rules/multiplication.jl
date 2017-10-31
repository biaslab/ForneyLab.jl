@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPMultiplicationOutVGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPMultiplicationInGVP)