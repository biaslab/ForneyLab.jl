export Wrap, wrap, wraps

# Wrap type and functions
type Wrap <: AbstractWrap
    id::Symbol
    source::TerminalNode
    sink::TerminalNode

    function Wrap(source::TerminalNode, sink::TerminalNode; id=symbol("$(source.id)_$(sink.id)"))
        hasNode(current_graph, source) || error("The source node does not belong to the current graph")
        hasNode(current_graph, sink) || error("The sink node does not belong to the current graph")
        !(sink in [wr.sink for wr in wraps(current_graph)]) || error("TerminalNode $(sink) already is a sink in another wrap")
        !is(source, sink) || error("Cannot create wrap: source and sink must be different nodes")
        !haskey(current_graph.wraps, id) || error("The wrap id $(id) already exists in the current graph. Consider specifying an explicit id.")

        wrap = new(id, source, sink)
        current_graph.wraps[id] = wrap
        return wrap
    end
end

show(io::IO, wrap::Wrap) = println(io, "Wrap with id $(wrap.id) from source $(wrap.source) to sink $(wrap.sink).")

wrap(id::Symbol, g::FactorGraph=current_graph) = g.wraps[id]

wraps(g::FactorGraph=current_graph) = Set{Wrap}(values(g.wraps))

function wraps(nd::TerminalNode, g::FactorGraph=current_graph)
    ws = Set{Wrap}()
    for w in values(g.wraps)
        (w.sink == nd) && push!(ws, w)
        (w.source == nd) && push!(ws, w)
    end
    return ws
end

hasWrap(graph::FactorGraph, wr::Wrap) = (haskey(graph.wraps, wr.id) && is(graph.wraps[wr.id], wr))

function Base.delete!(graph::FactorGraph, wr::Wrap)
    hasWrap(graph, wr) || error("Graph does not contain wrap")
    !graph.locked || error("Cannot delete node from locked graph")

    delete!(graph.wraps, wr.id)
    return graph
end