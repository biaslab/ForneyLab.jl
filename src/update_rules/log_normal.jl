@sumProductRule(:node_type     => LogNormal,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPLogNormalOutVPP)

@naiveVariationalRule(:node_type     => LogNormal,
                 	  :outbound_type => Message{LogNormal},
                 	  :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 	  :name          => VBLogNormalOut)
