@sumProductRule(:node_type     => LogNormal,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPLogNormalOutVPP)

@naiveVariationalRule(:node_type     => LogNormal,
                 	  :outbound_type => Message{LogNormal},
                 	  :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                 	  :name          => VBLogNormalOut)
