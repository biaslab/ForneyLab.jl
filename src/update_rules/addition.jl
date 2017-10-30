@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{Gaussian}}, Message{Univariate{Gaussian}}),
                :name          => SPAdditionOutGG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{Gaussian}}, Message{Univariate{Gaussian}}, Void),
                :name          => SPAdditionIn2GG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{Gaussian}}, Void, Message{Univariate{Gaussian}}),
                :name          => SPAdditionIn1GG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{Gaussian}}, Message{Univariate{PointMass}}),
                :name          => SPAdditionOutGP)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Void, Message{Univariate{PointMass}}, Message{Univariate{Gaussian}}),
                :name          => SPAdditionOutPG)

@sumProductRule(:node_type     => Addition,
                :outbound_type => Message{Univariate{Gaussian}},
                :inbound_types => (Message{Univariate{PointMass}}, Void, Message{Univariate{Gaussian}}),
                :name          => SPAdditionIn1PG)