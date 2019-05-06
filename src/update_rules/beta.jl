@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Beta},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPBetaOutNPP)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBBetaOut)