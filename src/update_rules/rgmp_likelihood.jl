@sumProductRule(:node_type     => RGMP_likelihood,
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{PointMass}, Nothing),
                :name          => SPRGMP_likelihoodIn1PN)
