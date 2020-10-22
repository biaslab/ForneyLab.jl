@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPTransitionOutNPP)

@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPTransitionIn1PNP)

@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Nothing, Message{Categorical}, Message{PointMass}),
                :name          => SPTransitionOutNCP)

@sumProductRule(:node_type     => Transition,
                :outbound_type => Message{Categorical},
                :inbound_types => (Message{Categorical}, Nothing, Message{PointMass}),
                :name          => SPTransitionIn1CNP)

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

@marginalRule(:node_type     => Transition,
              :inbound_types => (Message{Categorical}, Message{Categorical}, ProbabilityDistribution),
              :name          => MTransitionCCD)

@marginalRule(:node_type     => Transition,
              :inbound_types => (Message{Categorical}, Message{Categorical}, Nothing), # "Nothing" indicates marginalization
              :name          => MTransitionCCN)