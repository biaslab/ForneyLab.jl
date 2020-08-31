@sumProductRule(:node_type     => Gamma,
                :outbound_type => Message{Gamma},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGammaOutNPP)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGammaOut)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBGammaA)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBGammaB)
