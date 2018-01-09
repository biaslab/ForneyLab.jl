@sumProductRule(:node_type     => Dirichlet,
                :outbound_type => Message{Dirichlet},
                :inbound_types => (Void, Message{PointMass}),
                :name          => SPDirichletOutVP)

@naiveVariationalRule(:node_type     => Dirichlet,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (Void, ProbabilityDistribution),
                      :name          => VBDirichletOut)