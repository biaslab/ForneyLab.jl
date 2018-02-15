@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Void, Message{Categorical}, Message{PointMass}),
                :name          => SPTransitionOutVCP)

@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Message{Categorical}, Void, Message{PointMass}),
                :name          => SPTransitionIn1CVP)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Void, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBTransitionOut)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (ProbabilityDistribution, Void, ProbabilityDistribution),
                      :name          => VBTransitionIn1)

@naiveVariationalRule(:node_type     => Transition,
                      :outbound_type => Message{Dirichlet},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Void),
                      :name          => VBTransitionA)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Categorical},
                           :inbound_types => (Void, Message{Categorical}, ProbabilityDistribution),
                           :name          => SVBTransitionOutVCD)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Categorical},
                           :inbound_types => (Message{Categorical}, Void, ProbabilityDistribution),
                           :name          => SVBTransitionIn1CVD)

@structuredVariationalRule(:node_type     => Transition,
                           :outbound_type => Message{Dirichlet},
                           :inbound_types => (ProbabilityDistribution, Void),
                           :name          => SVBTransitionADV)

@marginalRule(:node_type => Transition,
              :inbound_types => (Message{Categorical}, Message{Categorical}, ProbabilityDistribution),
              :name => MTransitionCCD)