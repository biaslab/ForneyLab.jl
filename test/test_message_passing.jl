module MessagePassingTest

using Base.Test
import ForneyLab: ScheduleEntry, Schedule, summaryPropagationSchedule, Interface, generateId, addNode!, currentGraph, FactorGraph, Variable, FactorNode

# Integration helper
type MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(n_interfaces::Int; id=generateId(MockNode))
        self = new(id, Array(Interface, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = Interface(self)
        end

        return self
    end
end

@testset "ScheduleEntry" begin
    FactorGraph()
    nd = MockNode(2)

    # Schedule should be constructable by hand
    schedule = [ScheduleEntry(nd.i[1], Void), ScheduleEntry(nd.i[2], Void)]
    @test isa(schedule, Schedule)
end

@testset "summaryPropagationSchedule" begin
    # Construct
    #
    #    a    b    c
    # [1]--[2]--[3]--[4]
    #            |d
    #           [5]

    @test true == false
end

end # module