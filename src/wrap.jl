export Wrap, clearWraps, wrap, wraps

# Wrap type and functions
type Wrap <: AbstractWrap
    id::Symbol
    source::TerminalNode
    sink::TerminalNode

    function Wrap(source::TerminalNode, sink::TerminalNode; id=symbol("$(source.id)_$(sink.id)"))
        (source in values(current_graph.n) && sink in values(current_graph.n)) || error("The source and sink nodes should belong to the current graph")
        !(sink in [wr.sink for wr in wraps(current_graph)]) || error("TerminalNode $(sink) already is a sink in another wrap")
        !is(source, sink) || error("Cannot create wrap: source and sink must be different nodes")
        (!haskey(current_graph.wraps, id)) || error("The wrap id $(id) already exists in the current graph. Consider specifying an explicit id.")

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

function clearWraps(graph::FactorGraph=current_graph)
    graph.wraps = Dict{Symbol, Wrap}()
end