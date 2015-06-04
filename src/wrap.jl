export Wrap, clearWraps, wrap, wraps

# Wrap type and functions
type Wrap <: AbstractWrap
    id::Symbol
    from::TerminalNode
    to::TerminalNode

    function Wrap(from::TerminalNode, to::TerminalNode; id=symbol("$(from.id)_$(to.id)"))
        (from in values(current_graph.n) && to in values(current_graph.n)) || error("From and to node should belong to the current graph")
        !(to in [wr.to for wr in wraps(current_graph)]) || error("TerminalNode $(to) cannot receive multiple wraps")
        !is(from, to) || error("Cannot create wrap: from and to must be different nodes")
        (!haskey(current_graph.wraps, id)) || error("The wrap id $(id) already exists in the current graph. Consider specifying an explicit id.")

        wrap = new(id, from, to)
        current_graph.wraps[id] = wrap
        return wrap
    end
end

show(io::IO, wrap::Wrap) = println(io, "Wrap with id $(wrap.id) from $(wrap.from) to $(wrap.to).")

wrap(id::Symbol, g::FactorGraph=current_graph) = g.wraps[id]

wraps(g::FactorGraph=current_graph) = Set{Wrap}(values(g.wraps))
function wraps(nd::TerminalNode, g::FactorGraph=current_graph)
    ws = Set{Wrap}()
    for w in values(g.wraps)
        (w.to == nd) && push!(ws, w)
        (w.from == nd) && push!(ws, w)
    end
    return ws
end

function clearWraps(graph::FactorGraph=current_graph)
    graph.wraps = Dict{Symbol, Wrap}()
end