@sumProductRule(:node_type     => SVI,
                :outbound_type => Message{SetProbDist},
                :inbound_types => (Nothing, Message{FactorNode}),
                :name          => SPSVIOutNM)

@sumProductRule(:node_type     => SVI,
                :outbound_type => Message{FactorNode},
                :inbound_types => (Message{FactorNode}, Nothing),
                :name          => SPSVIIn1MN)
