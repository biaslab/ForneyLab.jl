@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixture1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gamma},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixture2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixture3)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gamma},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMixture4)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Bernoulli},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMixture5)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Void),
                 :name          => VBGaussianMixture6)