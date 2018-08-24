@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutVPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMPVP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutVGP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMGVP)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{GaussianMeanVariance},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceM)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{GaussianMeanVariance},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceOut)