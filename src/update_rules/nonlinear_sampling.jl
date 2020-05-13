@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{FactorFunction}, Nothing, Message{Gaussian}),
                :name          => SPNonlinearSIn1MNG)

@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{FactorFunction}, Message{Gaussian}, Nothing),
                :name          => SPNonlinearSIn2MGN)

@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPNonlinearSOutNGG)

@marginalRule(:node_type     => Nonlinear{Sampling},
              :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}), # "Nothing" indicates marginalization over outbound edge
              :name          => MNonlinearSOutNGG)

@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Nothing),
                :name          => SPNonlinearSIn1MN)

@sumProductRule(:node_type     => Nonlinear{Sampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message),
                :name          => SPNonlinearSOutNM)
