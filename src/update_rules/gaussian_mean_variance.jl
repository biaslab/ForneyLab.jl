@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutNPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMPNP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutNGP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanVarianceMGNP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Function},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing),
                :name          => SPGaussianMeanVarianceVGGN)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Function},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Nothing),
                :name          => SPGaussianMeanVarianceVPGN)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{SampleList}, Message{PointMass}),
                :name          => SPGaussianMeanVarianceOutNSP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{SampleList},Nothing,Message{PointMass}),
                :name          => SPGaussianMeanVarianceMSNP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian},Message{SampleList}),
                :name          => SPGaussianMeanVarianceOutNGS)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{Gaussian}, Nothing,Message{SampleList}),
                :name          => SPGaussianMeanVarianceMGNS)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{GaussianMeanVariance},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceM)

@naiveVariationalRule(:node_type     => GaussianMeanVariance,
                      :outbound_type => Message{GaussianMeanVariance},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGaussianMeanVarianceOut)
