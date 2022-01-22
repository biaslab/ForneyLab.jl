@sumProductRule(:node_type     => Categorical,
                :outbound_type => Message{Categorical},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPCategoricalOutNP)

@sumProductRule(:node_type     => Categorical,
                :outbound_type => Message{Dirichlet},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPCategoricalIn1PN)

@naiveVariationalRule(:node_type     => Categorical,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Nothing, Distribution),
                      :name          => VBCategoricalOut)

@naiveVariationalRule(:node_type     => Categorical,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (Distribution, Nothing),
                      :name          => VBCategoricalIn1)