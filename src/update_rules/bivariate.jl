@sumProductRule(:node_type     => Bivariate{Laplace},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Nothing, Message{Gaussian}),
                :name          => SPBivariateLIn1MNG)

@sumProductRule(:node_type     => Bivariate{Laplace},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Message{Gaussian}, Nothing),
                :name          => SPBivariateLIn2MGN)

@sumProductRule(:node_type     => Bivariate{Laplace},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}, Message{Gaussian}),
                :name          => SPBivariateLOutNGG)
