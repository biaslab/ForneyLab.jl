@sumProductRule(:node_type     => Dirichlet,
                :outbound_type => Message{Dirichlet},
                :inbound_types => (Void, Message{PointMass}),
                :name          => SPDirichletOutVP)

@variationalRule(:node_type     => Dirichlet,
                 :outbound_type => Message{Dirichlet},
                 :inbound_types => (Void, ProbabilityDistribution),
                 :name          => VBDirichletOut)