# This file holds the collapsedEP implementation for all supported expectation types

function collapsedEP{n_factors}(site::EPSite,
                                ::Type{PartitionedDistribution{Gaussian,n_factors}},
                                algo::ExpectationPropagation)
    # This function gets called when a ScheduleEntry{CollapsedExpectationRule} is executed
    # TODO
    println("Collapsed EP here...")
end