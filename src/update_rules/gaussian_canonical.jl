@sumProductRule(:node_type     => Gaussian{Canonical},
                :outbound_type => Message{Gaussian{Canonical}},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianCanonicalOutNPP)

@naiveVariationalRule(:node_type     => Gaussian{Canonical},
                      :outbound_type => Message{Gaussian{Canonical}},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBGaussianCanonicalOut)