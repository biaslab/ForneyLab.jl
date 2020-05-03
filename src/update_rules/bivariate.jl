@sumProductRule(:node_type     => Bivariate{Sampling},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{FactorFunction}, Nothing, Message{Gaussian}),
                :name          => SPBivariateSIn1MNG)

@sumProductRule(:node_type     => Bivariate{Sampling},
                :outbound_type => Message{Gaussian},
                :inbound_types => (Message{FactorFunction}, Message{Gaussian}, Nothing),
                :name          => SPBivariateSIn2MGN)

@sumProductRule(:node_type     => Bivariate{Sampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPBivariateSOutNGG)

@marginalRule(:node_type     => Bivariate{Sampling},
              :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}), # "Nothing" indicates marginalization over outbound edge
              :name          => MBivariateSOutNGG)
