@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Beta},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPBetaOutVPP)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBBetaOut)