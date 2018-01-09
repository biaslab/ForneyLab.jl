@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutVPP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Void, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMPVP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionOutVGP)

@sumProductRule(:node_type     => GaussianMeanPrecision,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPGaussianMeanPrecisionMGVP)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGaussianMeanPrecisionOut)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                      :name          => VBGaussianMeanPrecisionM)

@naiveVariationalRule(:node_type     => GaussianMeanPrecision,
                      :outbound_type => Message{Union{Gamma, Wishart}},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                      :name          => VBGaussianMeanPrecisionW)

# @structuredVariationalRule(:node_type     => GaussianMeanPrecision,
#                            :outbound_type => Message{Gaussian},
#                            :inbound_types => (Void, Message{Gaussian}, ProbabilityDistribution),
#                            :name          => SVBGaussianMeanPrecisionOutVGD)

# @structuredVariationalRule(:node_type     => GaussianMeanPrecision,
#                            :outbound_type => Message{Gaussian},
#                            :inbound_types => (Message{Gaussian}, Void, ProbabilityDistribution),
#                            :name          => SVBGaussianMeanPrecisionOutGVD)

# @structuredVariationalRule(:node_type     => GaussianMeanPrecision,
#                            :outbound_type => Message{Gamma},
#                            :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
#                            :name          => SVBGaussianMeanPrecisionW)