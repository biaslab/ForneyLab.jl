@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Void, Univariate, Univariate),
                 :name          => VBGaussianMeanPrecisionOut)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Univariate, Void, Univariate),
                 :name          => VBGaussianMeanPrecisionM)

@variationalRule(:node_type     => GaussianMeanPrecision,
                 :outbound_type => Message{Univariate{Gamma}},
                 :inbound_types => (Univariate, Univariate, Void),
                 :name          => VBGaussianMeanPrecisionW)