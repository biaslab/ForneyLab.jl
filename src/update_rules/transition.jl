@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Nothing, Message{Categorical}, Message{PointMass}),
                :name          => SPTransitionOutVCP)

@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Message{Categorical}, Nothing, Message{PointMass}),
                :name          => SPTransitionIn1CVP)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBTransitionOut)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBTransitionIn1)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBTransitionA)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Categorical},
                           :inbound_types => (Nothing, Message{Categorical}, ProbabilityDistribution),
                           :name          => SVBTransitionOutVCD)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Categorical},
                           :inbound_types => (Message{Categorical}, Nothing, ProbabilityDistribution),
                           :name          => SVBTransitionIn1CVD)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Dirichlet},
                           :inbound_types => (ProbabilityDistribution, Nothing),
                           :name          => SVBTransitionADV)

@marginalRule(:node_type => Transition,
              :inbound_types => (Message{Categorical}, Message{Categorical}, ProbabilityDistribution),
              :name => MTransitionCCD)