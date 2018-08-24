@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPMultiplicationOutVGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{Gaussian}),
                :name          => SPMultiplicationOutVPG)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPMultiplicationOutVPP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPMultiplicationIn1GVP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPMultiplicationIn1PVP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Nothing),
                :name          => SPMultiplicationAGPV)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Nothing),
                :name          => SPMultiplicationAPPV)