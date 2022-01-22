@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Categorical},
                      :inbound_types => (Nothing, Distribution, Distribution, Distribution),
                      :name          => VBSoftmaxOut)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{GaussianWeightedMeanPrecision},
                      :inbound_types => (Distribution, Nothing, Distribution, Distribution),
                      :name          => VBSoftmaxIn1)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Distribution, Nothing, Distribution),
                      :name          => VBSoftmaxXi)

@naiveVariationalRule(:node_type     => Softmax,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Distribution, Distribution, Nothing),
                      :name          => VBSoftmaxA)