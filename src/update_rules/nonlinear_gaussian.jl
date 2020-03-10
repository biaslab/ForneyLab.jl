@sumProductRule(:node_type     => NonlinearGaussian,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPNonlinearGaussianOutNGP)

@sumProductRule(:node_type     => NonlinearGaussian,
                :outbound_type => Message{GaussianMeanVariance},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPNonlinearGaussianIn1GGP)

@marginalRule(:node_type => NonlinearGaussian,
              :inbound_types => (Message{Gaussian}, Message{Gaussian}, ProbabilityDistribution),
              :name => MNonlinearGaussianGGD)