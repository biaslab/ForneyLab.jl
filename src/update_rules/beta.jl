# @sumProductRule(:node_type     => Beta,
#                 :outbound_type => Message{Beta},
#                 :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
#                 :name          => SPBetaOutNPP)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Union{Beta,SampleList}},
                :inbound_types => (Nothing, Message, Message),
                :name          => SPBetaOutNPP)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing, Message),
                :name          => SPBetaA)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Message, Nothing),
                :name          => SPBetaB)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBBetaOut)
