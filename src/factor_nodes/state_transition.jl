export StateTransition

type StateTransition <: SoftFactor
    id::Symbol
    interfaces::Vector{Interface}
    i::Dict{Symbol, Interface}
    inner_graph::FactorGraph
    terminals::Vector{Terminal}

    function StateTransition(x_prev::Variable, x::Variable, y::Variable; id=generateId(StateTransition))
        self = new(id, Array(Interface, 3), Dict{Symbol,Interface}())
        outer_graph = currentGraph()
        addNode!(outer_graph, self)
        self.i[:x_prev] = self.interfaces[1] = associate!(Interface(self), x_prev)
        self.i[:x] = self.interfaces[2] = associate!(Interface(self), x)
        self.i[:y] = self.interfaces[3] = associate!(Interface(self), y)

        # Build internal graph
        self.inner_graph = FactorGraph()
        self.terminals = Terminal[]
        let
            local x_prev = Variable(id=:x_prev)
            push!(self.terminals, Terminal(x_prev, id=:x_prev))
            local x = Variable(id=:x)
            push!(self.terminals, Terminal(x, id=:x))
            local y = Variable(id=:y)
            push!(self.terminals, Terminal(y, id=:y))

            Gaussian(x, x_prev, y)
        end
        setCurrentGraph(outer_graph)

        return self
    end
end