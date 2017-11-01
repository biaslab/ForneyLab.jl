@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Bernoulli},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixtureZ)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixtureM1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{AbstractGamma},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixtureW1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMixtureM2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{AbstractGamma},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void),
                 :name          => VBGaussianMixtureW2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixtureOut)