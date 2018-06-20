@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPMultiplicationOutVGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Void, Message{PointMass}, Message{Gaussian}),
                :name          => SPMultiplicationOutVPG)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPMultiplicationOutVPP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPMultiplicationIn1GVP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPMultiplicationIn1PVP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Void),
                :name          => SPMultiplicationAGPV)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Void),
                :name          => SPMultiplicationAPPV)