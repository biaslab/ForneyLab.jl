@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Nothing, ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBSoftmaxOut)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{GaussianWeightedMeanPrecision},
                      :inbound_types => (ProbabilityDistribution, Nothing, ProbabilityDistribution, ProbabilityDistribution),
                      :name          => VBSoftmaxIn1)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, Nothing, ProbabilityDistribution),
                      :name          => VBSoftmaxXi)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Function},
                      :inbound_types => (ProbabilityDistribution, ProbabilityDistribution, ProbabilityDistribution, Nothing),
                      :name          => VBSoftmaxA)