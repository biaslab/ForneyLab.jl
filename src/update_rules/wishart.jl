@sumProductRule(:node_type     => Wishart,
                :outbound_type => Message{Wishart},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPWishartOutNPP)

@naiveVariationalRule(:node_type     => Wishart,
                 	  :outbound_type => Message{Wishart},
                 	  :inbound_types => (Nothing, Distribution, Distribution),
                 	  :name          => VBWishartOut)