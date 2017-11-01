@variationalRule(:node_type     => Wishart,
                 :outbound_type => Message{Wishart},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBWishartOut)