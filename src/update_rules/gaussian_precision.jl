@sumProductRule(:node_type     => Gaussian{Precision},
                :outbound_type => Message{Gaussian{Precision}},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianPrecisionOutNPP)

@sumProductRule(:node_type     => Gaussian{Precision},
                :outbound_type => Message{Gaussian{Precision}},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPGaussianPrecisionMPNP)

@sumProductRule(:node_type     => Gaussian{Precision},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianPrecisionOutNGP)

@sumProductRule(:node_type     => Gaussian{Precision},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPGaussianPrecisionMGNP)

@naiveVariationalRule(:node_type     => Gaussian{Precision},
                      :outbound_type => Message{Gaussian{Precision}},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBGaussianPrecisionOut)

@naiveVariationalRule(:node_type     => Gaussian{Precision},
                      :outbound_type => Message{Gaussian{Precision}},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBGaussianPrecisionM)

@naiveVariationalRule(:node_type     => Gaussian{Precision},
                      :outbound_type => Message{Union{Gamma, Wishart}},
                      :inbound_types => (Distribution, Distribution, Nothing),
                      :name          => VBGaussianPrecisionW)

@structuredVariationalRule(:node_type     => Gaussian{Precision},
                           :outbound_type => Message{Gaussian{Moments}},
                           :inbound_types => (Nothing, Message{Gaussian}, Distribution),
                           :name          => SVBGaussianPrecisionOutVGD)

@structuredVariationalRule(:node_type     => Gaussian{Precision},
                           :outbound_type => Message{Gaussian{Moments}},
                           :inbound_types => (Message{Gaussian}, Nothing, Distribution),
                           :name          => SVBGaussianPrecisionMGVD)

@structuredVariationalRule(:node_type     => Gaussian{Precision},
                           :outbound_type => Message{Union{Gamma, Wishart}},
                           :inbound_types => (Distribution, Nothing),
                           :name          => SVBGaussianPrecisionW)

@marginalRule(:node_type => Gaussian{Precision},
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, Distribution),
              :name => MGaussian{Precision}GGD)

@marginalRule(:node_type => Gaussian{Precision},
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing), # Precision is marginalized out
              :name => MGaussian{Precision}GGN)
