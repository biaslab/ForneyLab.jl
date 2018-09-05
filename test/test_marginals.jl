module MarginalsTest

using Test
using ForneyLab

import ForneyLab: generateId, addNode!, associate!, Product

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Vector{Interface}(undef, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = associate!(Interface(self), vars[idx])
        end

        return self
    end
end

@testset "MarginalScheduleEntry" begin
    FactorGraph()
    x = Variable()
    nd = MockNode([x])

    # Marginal schedule should be constructable by hand
    marginal_schedule = [MarginalScheduleEntry(x, [nd.i[1]], Nothing)]
    @test isa(marginal_schedule, MarginalSchedule)
end

@testset "marginalSchedule" begin
    # Construct
    #
    #    a    b    
    # [1]--[2]--

    FactorGraph()

    a = Variable()
    b = Variable()

    n1 = MockNode([a])
    n2 = MockNode([a, b])

    marginal_schedule = marginalSchedule([a, b])

    @test length(marginal_schedule) == 2
    @test marginal_schedule[1].target == a
    @test marginal_schedule[1].interfaces[1] == n1.i[1]
    @test marginal_schedule[1].interfaces[2] == n2.i[1]
    @test marginal_schedule[1].marginal_update_rule == Product

    @test marginal_schedule[2].target == b
    @test marginal_schedule[2].interfaces[1] == n2.i[2]
    @test marginal_schedule[2].marginal_update_rule == Nothing
end

end # module
