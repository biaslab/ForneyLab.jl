@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}),
                :name          => SPNonlinearOutNG)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing),
                :name          => SPNonlinearIn1GG)

@marginalRule(:node_type => Nonlinear,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}),
              :name => MNonlinearGG)