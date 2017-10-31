@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecisionOut)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecisionM)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{AbstractGamma},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                 :name          => VBGaussianMeanPrecisionW)