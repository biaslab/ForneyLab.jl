@sumProductRule(:node_type     => Gamma,
                :outbound_type => Message{Gamma},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGammaOutVPP)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGammaOut)