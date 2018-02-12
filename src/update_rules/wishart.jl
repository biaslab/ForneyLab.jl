@sumProductRule(:node_type     => Wishart,
                :outbound_type => Message{Wishart},
                :inbound_types => (Void, Message{PointMass}, Message{PointMass}),
                :name          => SPWishartOutVPP)

@naiveVariationalRule(:node_type     => Wishart,
                 	  :outbound_type => Message{Wishart},
                 	  :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 	  :name          => VBWishartOut)