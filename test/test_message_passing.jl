module MessagePassingTest

using Base.Test
using ForneyLab
import ForneyLab: generateId, addNode!, associate!, summaryPropagationSchedule, matches, flatten

@testset "Message" begin
    msg = Message(Univariate, Gaussian, m=0.0, v=1.0)
    @test isa(msg, Message{Gaussian})
    @test isa(msg, Message{Gaussian, Univariate})
    @test !isa(msg, Message{PointMass})
    @test !isa(msg, Message{Gaussian, Multivariate})
    @test !isa(msg, Message{Gaussian, MatrixVariate})
    @test msg.dist == ProbabilityDistribution(Univariate, Gaussian, m=0.0, v=1.0)

    @test Message(Univariate, PointMass, m=0.0) == Message(Univariate, PointMass, m=0.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, Gaussian, m=0.0, v=1.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, PointMass, m=1.0)
end

@testset "matches" begin
    @test matches(Message{Gaussian, Univariate}, Message{Gaussian, Univariate})
    @test !matches(Message{Gaussian, Univariate}, Message{Gaussian, Multivariate})
    @test !matches(Message{Gaussian, Univariate}, Message{PointMass, Univariate})
    @test matches(Message{Gaussian, Univariate}, Message{Gaussian})
    @test matches(Message{Gaussian}, Message{Gaussian})
    @test !matches(Void, Message{Gaussian})
    @test matches(Message{Gamma, Univariate}, Message{Union{Gamma, Wishart}, Univariate})
    @test matches(Message{Gamma}, Message{Union{Gamma, Wishart}})
end

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Array{Interface}(n_interfaces), Dict{Int,Interface}())
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

    @test length(schedule) == 5
    @test ScheduleEntry(n1.i[1], Void) in schedule
    @test ScheduleEntry(n2.i[2], Void) in schedule
    @test ScheduleEntry(n4.i[1], Void) in schedule
    @test ScheduleEntry(n3.i[3], Void) in schedule
    @test ScheduleEntry(n5.i[1], Void) in schedule
end

@testset "flatten" begin
    FactorGraph()

    a = Variable()
    b = Variable()
    c = Variable()

    n1 = MockNode([a])
    n2 = MockNode([a, b])
    n3 = MockNode([b, c])
    n4 = MockNode([Variable()])

    schedule = [ScheduleEntry(n1.i[1], Void),
                ScheduleEntry(n2.i[2], Void),
                ScheduleEntry(n3.i[2], Void)]
    schedule[2].internal_schedule = [   ScheduleEntry(n1.i[1], Void),
                                        ScheduleEntry(n4.i[1], Void)]

    flat_schedule = flatten(schedule)

    @test length(flat_schedule) == 3
    @test flat_schedule[1] == schedule[1]
    @test flat_schedule[2] == schedule[2].internal_schedule[2]
    @test flat_schedule[3] == schedule[3]
end


end # module