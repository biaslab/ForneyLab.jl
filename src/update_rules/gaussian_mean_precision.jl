@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutVPP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMPVP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutVGP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMGVP)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecisionOut)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecisionM)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Union{Gamma, Wishart}},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                 :name          => VBGaussianMeanPrecisionW)