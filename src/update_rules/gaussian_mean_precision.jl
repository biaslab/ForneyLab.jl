@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutNPP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{GaussianMeanPrecision},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMPNP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutNGP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMGNP)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{GaussianMeanPrecision},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGaussianMeanPrecisionOut)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{GaussianMeanPrecision},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBGaussianMeanPrecisionM)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{Union{Gamma, Wishart}},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBGaussianMeanPrecisionW)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{GaussianMeanVariance},
                           :inbound_types => (Nothing, Message{Gaussian}, ProbabilityDistribution),
                           :name          => SVBGaussianMeanPrecisionOutVGD)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{GaussianMeanVariance},
                           :inbound_types => (Message{Gaussian}, Nothing, ProbabilityDistribution),
                           :name          => SVBGaussianMeanPrecisionMGVD)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{Union{Gamma, Wishart}},
                           :inbound_types => (ProbabilityDistribution, Nothing),
                           :name          => SVBGaussianMeanPrecisionW)

@marginalRule(:node_type => GaussianMeanPrecision,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution),
              :name => MGaussianMeanPrecisionGGD)

@marginalRule(:node_type => GaussianMeanPrecision,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing), # Precision is marginalized out
              :name => MGaussianMeanPrecisionGGN)
