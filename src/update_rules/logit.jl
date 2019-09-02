@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBLogitOut)

@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{GaussianWeightedMeanPrecision},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBLogitIn1)

@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBLogitXi)