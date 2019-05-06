@sumProductRule(:node_type     => Categorical,
                :outbound_type => Message{Categorical},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPCategoricalOutNP)

@naiveVariationalRule(:node_type     => Categorical,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Nothing, ProbabilityDistribution),
                      :name          => VBCategoricalOut)

@naiveVariationalRule(:node_type     => Categorical,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (ProbabilityDistribution, Nothing),
                      :name          => VBCategoricalIn1)