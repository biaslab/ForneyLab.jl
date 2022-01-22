@sumProductRule(:node_type     => Gamma,
                :outbound_type => Message{Gamma},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGammaOutNPP)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBGammaOut)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Function},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBGammaA)

@naiveVariationalRule(:node_type     => Gamma,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Distribution, Distribution, Nothing),
                      :name          => VBGammaB)
