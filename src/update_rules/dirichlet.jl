@sumProductRule(:node_type     => Dirichlet,
                :outbound_type => Message{Dirichlet},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPDirichletOutNP)

@naiveVariationalRule(:node_type     => Dirichlet,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBDirichletOut)