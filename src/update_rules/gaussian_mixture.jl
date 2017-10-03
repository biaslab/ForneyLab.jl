@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :outbound_id   => 1,
                 :name          => VBGaussianMixture1)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gamma},
                 :outbound_id   => 2,
                 :name          => VBGaussianMixture2)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :outbound_id   => 3,
                 :name          => VBGaussianMixture3)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gamma},
                 :outbound_id   => 4,
                 :name          => VBGaussianMixture4)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Bernoulli},
                 :outbound_id   => 5,
                 :name          => VBGaussianMixture5)

@variationalRule(:node_type     => GaussianMixture,
                 :outbound_type => Message{Gaussian},
                 :outbound_id   => 6,
                 :name          => VBGaussianMixture6)