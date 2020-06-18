@sumProductRule(:node_type     => Dirichlet,
                :outbound_type => Message{Dirichlet},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPDirichletOutNP)

@naiveVariationalRule(:node_type     => Dirichlet,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBDirichletOut)

@naiveVariationalRule(:node_type     => Dirichlet,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, Nothing),
                      :name          => VBDirichletIn1)
