export children

# ChildrenVisitor and related methods
# Used by the children(...) algorithm

type ChildrenVisitor{V} <: AbstractGraphVisitor
    vertices::Vector{V}
    breakers::Set{V}
    restrict_to::Set{V}
    allow_cycles::Bool

    function ChildrenVisitor(n::Int; allow_cycles::Bool=false, breakers::Set{V}=Set{V}(), restrict_to::Set{V}=Set{V}())
        vs = Array(Int, 0)
        sizehint!(vs, n)
        new(vs, breakers, restrict_to, allow_cycles)
    end
end

function Graphs.examine_neighbor!{V}(visitor::ChildrenVisitor{V}, u::V, v::V, vcolor::Int, ecolor::Int)
    if visitor.allow_cycles == false && vcolor == 1 && ecolor == 0
        throw(ArgumentError("The input graph contains a loop around $(v)."))
    end
end

function Graphs.discover_vertex!{V}(visitor::ChildrenVisitor{V}, v::V)
    (v in visitor.breakers) && return false # stop at breakers
    if length(visitor.restrict_to) > 0
        (v in visitor.restrict_to) || return false
    end

    return true
end

Graphs.close_vertex!{V}(visitor::ChildrenVisitor{V}, v::V) = push!(visitor.vertices, v)


"""
children(vertices, graph; allow_cycles=false, breakers=[], restrict_to=[])

Return a vector consisting of `vertices` and all their children in `graph`.
`v` is a child of `u` iff there exists a path from `u` to `v`.
The resulting array is sorted in reverse topological order.
For each directed edge `(u,v)`, `v` (child of `u`) appears before `u`.

Optional keyword arguments:

- `allow_cycles`: set to true to accept cycles.
- `breakers`: a Set of vertices on which the search will terminate.
- `restrict_to`: a Set of vertices to restrict the search to.

This function can be used to generate message passing schedules
if `graph` is a dependency graph.
"""
function children{V}(   vertices::Vector{V},
                        graph::AbstractGraph{V};
                        kwargs...)
    @graph_requires graph vertex_list incidence_list vertex_map

    cmap = zeros(Int, num_vertices(graph))
    visitor = ChildrenVisitor{V}(num_vertices(graph); kwargs...)

    for s in vertices
        if cmap[vertex_index(s, graph)] == 0
            traverse_graph(graph, DepthFirst(), s, visitor, vertexcolormap=cmap)
        end
    end

    visitor.vertices
end