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

    n1 = MockNode([a])
    n2 = MockNode([a, b])
    n3 = MockNode([b, c, d])
    n4 = MockNode([c])
    n5 = MockNode([d])

    schedule = summaryPropagationSchedule(d)

    @test schedule == [ ScheduleEntry(n1.i[1], Void),
                        ScheduleEntry(n2.i[2], Void),
                        ScheduleEntry(n4.i[1], Void),
                        ScheduleEntry(n3.i[3], Void),
                        ScheduleEntry(n5.i[1], Void)]
end

end # module