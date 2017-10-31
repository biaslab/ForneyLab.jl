@variationalRule(:node_type     => Wishart,
                 :outbound_type => Message{AbstractGamma},
                 :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                 :name          => VBWishartOut)