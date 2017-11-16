@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPMultiplicationOutVGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPMultiplicationOutVPG)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPMultiplicationIn1GVP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Void),
                :name          => SPMultiplicationIn2GPV)