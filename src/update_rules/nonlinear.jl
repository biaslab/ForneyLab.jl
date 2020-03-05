@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPNonlinearOutNGP)

@sumProductRule(:node_type     => Nonlinear,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPNonlinearIn1GGP)

@marginalRule(:node_type => Nonlinear,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution),
              :name => MNonlinearGGD)