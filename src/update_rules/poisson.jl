# Register sumProductRule: SPPoissonOutNP, SPPoissonLPN
@sumProductRule(:node_type     => Poisson,
                :outbound_type => Message{Poisson},
                :inbound_types => (Nothing, Message{PointMass}),
                :name          => SPPoissonOutNP)

@sumProductRule(:node_type     => Poisson,
                :outbound_type => Message{Gamma},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPPoissonLPN)

# Register naiveVariationalRule: VBPoissonOut, VBPoissonL

@naiveVariationalRule(:node_type     => Poisson,
                      :outbound_type => Message{Poisson},
                      :inbound_types => (Nothing, Distribution),
                      :name          => VBPoissonOut)

@naiveVariationalRule(:node_type     => Poisson,
                      :outbound_type => Message{Gamma},
                      :inbound_types => (Distribution, Nothing),
                      :name          => VBPoissonL)
