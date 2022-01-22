@sumProductRule(:node_type     => SampleList,
                :outbound_type => Message{SampleList},
                :inbound_types => (Nothing, Message{PointMass}, Message{PointMass}),
                :name          => SPSampleListOutNPP)

@naiveVariationalRule(:node_type     => SampleList,
                      :outbound_type => Message{SampleList},
                      :inbound_types => (Nothing, Distribution, Distribution),
                      :name          => VBSampleListOut)
