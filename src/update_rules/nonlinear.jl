@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearUTOutNG)

@sumProductRule(:node_type     => Nonlinear{Unscented},
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearUTIn1GG)

@sumProductRule(:node_type     => Nonlinear{ImportanceSampling},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Nothing),
                :name          => SPNonlinearISInMN)

@sumProductRule(:node_type     => Nonlinear{ImportanceSampling},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearISOutNG)

@sumProductRule(:node_type     => Nonlinear{Laplace},
                :outbound_type => Message{Function},
                :inbound_types => (Message{FactorFunction}, Nothing),
                :name          => SPNonlinearLInMN)

@sumProductRule(:node_type     => Nonlinear{Laplace},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearLOutNG)
