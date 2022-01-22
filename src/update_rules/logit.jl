@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{Bernoulli},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBLogitOut)

@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{GaussianWeightedMeanPrecision},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBLogitIn1)

@naiveVariationalRule(:node_type     => Logit,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Distribution, Nothing),
                      :name          => VBLogitXi)