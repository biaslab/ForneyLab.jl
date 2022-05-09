@sumProductRule(:node_type     => LogNormal,
                :outbound_type => Message{LogNormal},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPLogNormalOutNPP)

@naiveVariationalRule(:node_type     => LogNormal,
                 	  :outbound_type => Message{LogNormal},
                 	  :inbound_types => (Nothing, Distribution, Distribution),
                 	  :name          => VBLogNormalOut)
