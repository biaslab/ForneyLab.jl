@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Union{Beta,SampleList}},
                :inbound_types => (Nothing, Message, Message),
                :name          => SPBetaOutNMM)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Nothing, Message),
                :name          => SPBetaMNM)

@sumProductRule(:node_type     => Beta,
                :outbound_type => Message{Function},
                :inbound_types => (Message, Message, Nothing),
                :name          => SPBetaMMN)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Beta},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBBetaOut)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBBetaA)

@naiveVariationalRule(:node_type     => Beta,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBBetaB)
