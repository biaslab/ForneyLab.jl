export Wrap, wrap, wraps, hasWrap

# Wrap type and functions
type Wrap <: AbstractWrap
    id::Symbol
    tail::TerminalNode
    head::TerminalNode
    tail_buffer::Vector{ProbabilityDistribution}
    head_buffer::Vector{ProbabilityDistribution}
    
    function Wrap(tail::TerminalNode, head::TerminalNode; id=symbol("$(tail.id)_$(head.id)"), block_size::Int64=0)
        current_graph = currentGraph()
        hasNode(current_graph, tail) || error("The tail node does not belong to the current graph")
        hasNode(current_graph, head) || error("The head node does not belong to the current graph")
        !(head in [wr.head for wr in wraps(current_graph)]) || error("TerminalNode $(head) already is a head in another wrap")
        !is(tail, head) || error("Cannot create wrap: tail and head must be different nodes")
        !haskey(current_graph.wraps, id) || error("The wrap id $(id) already exists in the current graph. Consider specifying an explicit id.")

        if !isdefined(current_graph, :block_size)
            if block_size >= 1 
                # Test whether there are wraps without defined block_size
                do_wraps_exist = isempty(wraps(current_graph)) || error("The graph contains wraps both with defined and undefined block_size")
                current_graph.block_size = block_size
                current_graph.current_section = 1
                tail_buffer = Vector{ProbabilityDistribution}(block_size)
                head_buffer = Vector{ProbabilityDistribution}(block_size)
            else
                tail_buffer = Vector{ProbabilityDistribution}(1)
                head_buffer = Vector{ProbabilityDistribution}(1)
            end
        else
            if current_graph.block_size != block_size
                error("The graph contains two wraps with different block sizes.")
            else
                tail_buffer = Vector{ProbabilityDistribution}(block_size)
                head_buffer = Vector{ProbabilityDistribution}(block_size)
                current_graph.current_section = 1
            end
        end

        if isdefined(current_graph, :block_size)
            for (component, write_buffer) in current_graph.write_buffers
                length(write_buffer) == block_size || error("Write buffer attached to $component has length $(length(write_buffer)), which is not equal to the block_size: $block_size.") 
            end
        end
        
        wrap = new(id, tail, head, tail_buffer, head_buffer)
        current_graph.wraps[id] = wrap
        return wrap
    end
end

show(io::IO, wrap::Wrap) = println(io, "Wrap with id $(wrap.id) from tail $(wrap.tail) to head $(wrap.head).")

wrap(id::Symbol, g::FactorGraph=current_graph) = g.wraps[id]

wraps(g::FactorGraph=current_graph) = Set{Wrap}(values(g.wraps))

function wraps(nd::TerminalNode, g::FactorGraph=current_graph)
    hasNode(g, nd) ||  error("The specified node does not belong to the specified or current graph")
    ws = Set{Wrap}()
    for w in values(g.wraps)
        (w.head == nd) && push!(ws, w)
        (w.tail == nd) && push!(ws, w)
    end
    return ws
end

hasWrap(graph::FactorGraph, wr::Wrap) = (haskey(graph.wraps, wr.id) && is(graph.wraps[wr.id], wr))
