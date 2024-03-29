module MessagePassingTest

using Test
using ForneyLab
using ForneyLab: generateId, addNode!, associate!, summaryPropagationSchedule, <<, matches, flatten, Cluster

@testset "Message" begin
    msg = Message(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test isa(msg, Message{Gaussian{Moments}})
    @test isa(msg, Message{Gaussian{Moments}, Univariate})
    @test !isa(msg, Message{PointMass})
    @test !isa(msg, Message{Gaussian{Precision}})
    @test !isa(msg, Message{Gaussian{Moments}, Multivariate})
    @test !isa(msg, Message{Gaussian{Moments}, MatrixVariate})
    @test msg.dist == Distribution(Univariate, Gaussian{Moments}, m=0.0, v=1.0)

    @test Message(Univariate, PointMass, m=0.0) == Message(Univariate, PointMass, m=0.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
    @test Message(Univariate, PointMass, m=0.0) != Message(Univariate, PointMass, m=1.0)
end

@testset "<<" begin
    @test <<(Message{Gaussian}, Message)
    @test <<(Message{Gaussian, Univariate}, Message)
    @test <<(Message{Gaussian, Univariate}, Message{Gaussian, Univariate})
    @test <<(Message{Gaussian{Moments}, Univariate}, Message{Gaussian, Univariate})
    @test !<<(Message{Gaussian{Moments}, Univariate}, Message{Gaussian, Multivariate})
    @test !<<(Message{Gaussian{Moments}, Univariate}, Message{PointMass, Univariate})
    @test <<(Message{Gaussian{Moments}, Univariate}, Message{Gaussian})
    @test <<(Message{Gaussian}, Message{Gaussian})
    @test <<(Message{Gaussian{Moments}}, Message{Gaussian})
    @test !<<(Nothing, Message{Gaussian})
    @test !<<(Nothing, Message{Gaussian{Moments}})
    @test <<(Message{Gamma, Univariate}, Message{Union{Gamma, Wishart}, Univariate})
    @test <<(Message{Gamma}, Message{Union{Gamma, Wishart}})
    @test !<<(Message{Gaussian{Moments}, Univariate}, Distribution{Univariate, Gaussian{Moments}})
    @test !<<(Distribution{Univariate, Gaussian{Moments}}, Message{Gaussian{Moments}, Univariate})
end

@testset "Message constructor alias" begin
    @test M(Univariate, Gaussian{Moments}, m=0.0, v=1.0) == Message(Univariate, Gaussian{Moments}, m=0.0, v=1.0)
end

@testset "matches" begin
    @test matches(Message{Gaussian}, Message{FactorFunction})
    @test matches(Message{FactorFunction}, Message{Gaussian})
    @test !matches(Nothing, Message{Gaussian})
end

# Integration helper
mutable struct MockNode <: FactorNode
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Int,Interface}

    function MockNode(vars::Vector{Variable}; id=generateId(MockNode))
        n_interfaces = length(vars)
        self = new(id, Array{Interface}(undef, n_interfaces), Dict{Int,Interface}())
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
    schedule = [ScheduleEntry(nd.i[1], Nothing), ScheduleEntry(nd.i[2], Nothing)]
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

    cl = Cluster(n1, [a.edges; b.edges])

    schedule = summaryPropagationSchedule(Variable[d], Cluster[cl])

    @test length(schedule) == 7
    @test ScheduleEntry(n1.i[1], Nothing) in schedule
    @test ScheduleEntry(n2.i[2], Nothing) in schedule
    @test ScheduleEntry(n4.i[1], Nothing) in schedule
    @test ScheduleEntry(n3.i[3], Nothing) in schedule
    @test ScheduleEntry(n5.i[1], Nothing) in schedule
    @test ScheduleEntry(n3.i[1], Nothing) in schedule
    @test ScheduleEntry(n2.i[1], Nothing) in schedule
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

    schedule = [ScheduleEntry(n1.i[1], Nothing),
                ScheduleEntry(n2.i[2], Nothing),
                ScheduleEntry(n3.i[2], Nothing)]
    schedule[2].internal_schedule = [   ScheduleEntry(n1.i[1], Nothing),
                                        ScheduleEntry(n4.i[1], Nothing)]

    flat_schedule = flatten(schedule)

    @test length(flat_schedule) == 3
    @test flat_schedule[1] == schedule[1]
    @test flat_schedule[2] == schedule[2].internal_schedule[2]
    @test flat_schedule[3] == schedule[3]
end


end # module