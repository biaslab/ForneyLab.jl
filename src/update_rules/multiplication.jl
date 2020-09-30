@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPMultiplicationOutNGP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{Gaussian}),
                :name          => SPMultiplicationOutNPG)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPMultiplicationOutNPP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPMultiplicationIn1GNP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPMultiplicationIn1PNP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{GaussianWeightedMeanPrecision},
                :inbound_types => (Message{Gaussian}, Message{PointMass}, Nothing),
                :name          => SPMultiplicationAGPN)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{PointMass},
                :inbound_types => (Message{PointMass}, Message{PointMass}, Nothing),
                :name          => SPMultiplicationAPPN)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gamma},
                :inbound_types => (Nothing, Message{Union{Gamma, Wishart}}, Message{PointMass}),
                :name          => SPMultiplicationOutNΓP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gamma},
                :inbound_types => (Nothing, Message{PointMass}, Message{Union{Gamma, Wishart}}),
                :name          => SPMultiplicationOutNPΓ)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gamma},
                :inbound_types => (Message{Union{Gamma, Wishart}}, Nothing, Message{PointMass}),
                :name          => SPMultiplicationIn1ΓNP)

@sumProductRule(:node_type     => Multiplication,
                :outbound_type => Message{Gamma},
                :inbound_types => (Message{Union{Gamma, Wishart}}, Message{PointMass}, Nothing),
                :name          => SPMultiplicationAΓPN)
