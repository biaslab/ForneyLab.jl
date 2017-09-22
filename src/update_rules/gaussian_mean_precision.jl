@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :outbound_id   => 1,
                 :name          => VBGaussianMeanPrecision1)

# TODO: how to handle Wishart
@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gamma},
                 :outbound_id   => 2,
                 :name          => VBGaussianMeanPrecision2)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Gaussian},
                 :outbound_id   => 3,
                 :name          => VBGaussianMeanPrecision3)