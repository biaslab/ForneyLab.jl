@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Void, Message{Gaussian}, Message{PointMass}),
                :name          => SPNonlinearOutVGP)

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
