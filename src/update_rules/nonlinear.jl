@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPNonlinearOutVGP)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{Gaussian}, Void, Message{PointMass}),
                :name          => SPNonlinearIn1GVP)

@naiveVariationalRule(:node_type     => Nonlinear,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBNonlinearOut)

@naiveVariationalRule(:node_type     => Nonlinear,
                      :outbound_type => Message{Gaussian},
                      :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                      :name          => VBNonlinearIn1)

@naiveVariationalRule(:node_type     => Nonlinear,
                      :outbound_type => Message{Union{Gamma, Wishart}},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                      :name          => VBNonlinearW)

@structuredVariationalRule(:node_type     => Nonlinear,
                           :outbound_type => Message{Gaussian},
                           :inbound_types => (Void, Message{Gaussian}, ProbabilityDistribution),
                           :name          => SVBNonlinearOutVGD)

@structuredVariationalRule(:node_type     => Nonlinear,
                           :outbound_type => Message{Gaussian},
                           :inbound_types => (Message{Gaussian}, Void, ProbabilityDistribution),
                           :name          => SVBNonlinearIn1GVD)

@structuredVariationalRule(:node_type     => Nonlinear,
                           :outbound_type => Message{Union{Gamma, Wishart}},
                           :inbound_types => (ProbabilityDistribution, Void),
                           :name          => SVBNonlinearW)

@marginalRule(:node_type => Nonlinear,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution),
              :name => MNonlinearGGD)