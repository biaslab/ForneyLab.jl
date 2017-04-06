module MessagePassingTest

using Base.Test
import ForneyLab: ScheduleEntry, Schedule, summaryPropagationSchedule, Interface, generateId, addNode!, currentGraph, FactorGraph, Variable, FactorNode, associate!

# Integration helper
type MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Array(Interface, n_interfaces), Dict{Int,Interface}())
        addNode!(currentGraph(), self)

        for idx = 1:n_interfaces
            self.i[idx] = self.interfaces[idx] = associate!(Interface(self), vars[idx])
        end

        return self
    end
end

@testset "ScheduleEntry" begin
    FactorGraph()
    nd = MockNode([Variable(), Variable()])

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

    FactorGraph()

    a = Variable()
    b = Variable()
    c = Variable()
    d = Variable()

    MockNode([a])
    MockNode([a, b])
    MockNode([b, c, d])
    MockNode([c])
    MockNode([d])

    schedule = summaryPropagationSchedule([d])
    println(schedule)

    @test true == false
end

end # module