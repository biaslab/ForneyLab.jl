@variationalRule(:node_type     => LogNormal,
                 :outbound_type => Message{LogNormal},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBLogNormalOut)
