@sumProductRule(:node_type     => Gamma,
                :outbound_type => Message{Gamma},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPGammaOutVPP)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBGammaOut)