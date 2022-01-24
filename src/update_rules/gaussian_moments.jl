@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPGaussianMomentsOutNPP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Message{PointMass}, Nothing, Message{PointMass}),
                :name          => SPGaussianMomentsMPNP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Nothing, Message{Gaussian}, Message{PointMass}),
                :name          => SPGaussianMomentsOutNGP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Gaussian{Moments}},
                :inbound_types => (Message{Gaussian}, Nothing, Message{PointMass}),
                :name          => SPGaussianMomentsMGNP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Function},
                :inbound_types => (Message{Gaussian}, Message{Gaussian}, Nothing),
                :name          => SPGaussianMomentsVGGN)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{Function},
                :inbound_types => (Message{PointMass}, Message{Gaussian}, Nothing),
                :name          => SPGaussianMomentsVPGN)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{SampleList}, Message{PointMass}),
                :name          => SPGaussianMomentsOutNSP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{SampleList},Nothing,Message{PointMass}),
                :name          => SPGaussianMomentsMSNP)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{Gaussian},Message{SampleList}),
                :name          => SPGaussianMomentsOutNGS)

@sumProductRule(:node_type     => Gaussian{Moments},
                :outbound_type => Message{SampleList},
                :inbound_types => (Message{Gaussian}, Nothing,Message{SampleList}),
                :name          => SPGaussianMomentsMGNS)

@naiveVariationalRule(:node_type     => Gaussian{Moments},
                      :outbound_type => Message{Gaussian{Moments}},
                      :inbound_types => (Distribution, Nothing, Distribution),
                      :name          => VBGaussianMomentsM)

@naiveVariationalRule(:node_type     => Gaussian{Moments},
                      :outbound_type => Message{Gaussian{Moments}},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBGaussianMomentsOut)
