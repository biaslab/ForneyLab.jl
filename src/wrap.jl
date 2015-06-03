export Wrap, clearWraps, wraps

# Wrap type and functions
type Wrap <: AbstractWrap
    from::TerminalNode
    to::TerminalNode

    function Wrap(from::TerminalNode, to::TerminalNode)
        (from in values(current_graph.n) && to in values(current_graph.n)) || error("From and to node should belong to the current graph.")
        !is(from, to) || error("Cannot create wrap: from and to must be different nodes")
        wrap = new(from, to)
        push!(current_graph.wraps, wrap)
        return wrap
    end
end

wraps(g::FactorGraph=current_graph) = g.wraps

function clearWraps(graph::FactorGraph=current_graph)
    graph.wraps = Array(Wrap, 0)
end