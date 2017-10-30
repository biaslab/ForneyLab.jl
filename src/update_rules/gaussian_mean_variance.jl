@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{PointMass}}, Message{Univariate{PointMass}}),
                :name          => SPGaussianMeanVarianceOutPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{PointMass}}, Void, Message{Univariate{PointMass}}),
                :name          => SPGaussianMeanVarianceMPP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{Gaussian}}, Message{Univariate{PointMass}}),
                :name          => SPGaussianMeanVarianceOutGP)

@sumProductRule(:node_type     => GaussianMeanVariance,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{Gaussian}}, Void, Message{Univariate{PointMass}}),
                :name          => SPGaussianMeanVarianceMGP)

@variationalRule(:node_type     => GaussianMeanVariance,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Univariate, Void, Univariate),
                 :name          => VBGaussianMeanVarianceM)

@variationalRule(:node_type     => GaussianMeanVariance,
                 :outbound_type => Message{Univariate{Gaussian}},
                 :inbound_types => (Void, Univariate, Univariate),
                 :name          => VBGaussianMeanVarianceOut)