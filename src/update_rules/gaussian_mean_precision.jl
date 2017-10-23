@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecision1)

# TODO: how to handle Wishart
@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gamma},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMeanPrecision2)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                 :name          => VBGaussianMeanPrecision3)