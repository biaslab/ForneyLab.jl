@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutVPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMPVP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutVGP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMGVP)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceM)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceOut)