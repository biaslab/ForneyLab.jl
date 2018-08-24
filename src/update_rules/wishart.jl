@sumProductRule(:node_type     => Wishart,
                :outbound_type => Message{Wishart},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPWishartOutVPP)

@naiveVariationalRule(:node_type     => Wishart,
                 	  :outbound_type => Message{Wishart},
                 	  :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                 	  :name          => VBWishartOut)