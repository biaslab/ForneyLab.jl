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
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBGaussianMeanPrecisionOut)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{GaussianMeanPrecision},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBGaussianMeanPrecisionM)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{Union{Gamma, Wishart}},
                      :inbound_types => (Distribution, Distribution, Nothing),
                      :name          => VBGaussianMeanPrecisionW)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{GaussianMeanVariance},
                           :inbound_types => (Nothing, Message{Gaussian}, Distribution),
                           :name          => SVBGaussianMeanPrecisionOutVGD)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{GaussianMeanVariance},
                           :inbound_types => (Message{Gaussian}, Nothing, Distribution),
                           :name          => SVBGaussianMeanPrecisionMGVD)

@structuredVariationalRule(:node_type     => GaussianMeanPrecision,
                           :outbound_type => Message{Union{Gamma, Wishart}},
                           :inbound_types => (Distribution, Nothing),
                           :name          => SVBGaussianMeanPrecisionW)

@marginalRule(:node_type => GaussianMeanPrecision,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, Distribution),
              :name => MGaussianMeanPrecisionGGD)

@marginalRule(:node_type => GaussianMeanPrecision,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing), # Precision is marginalized out
              :name => MGaussianMeanPrecisionGGN)
