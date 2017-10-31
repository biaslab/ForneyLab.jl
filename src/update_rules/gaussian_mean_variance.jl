@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutGP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMGP)

@variationalRule(:node_type     => GaussianMeanVariance,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                 :name          => VBGaussianMeanVarianceM)

@variationalRule(:node_type     => GaussianMeanVariance,
                 :outbound_type => Message{Gaussian},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBGaussianMeanVarianceOut)