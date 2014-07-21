# Functions for visualizing graphs
export graph2dot, graphviz, showDot

function getEdges(nodes::Array{Node,1})
    # Return a set of all Edges connected to nodes
    edges = Set()
    for node in nodes
        for interface in node.interfaces
            (interface.edge==nothing) || push!(edges, interface.edge)
        end
    end
    return edges
end

function graph2dot(nodes::Array{Node,1})
    # Return a string representing the graph that connects the nodes in DOT format for visualization.
    # http://en.wikipedia.org/wiki/DOT_(graph_description_language)
    node_type_symbols = {   AdditionNode => "+",
                            EqualityNode => "="}
    edges = getEdges(nodes)
    
    dot = "digraph G{splines=true;sep=\"+25,25\";overlap=scalexy;nodesep=1.6;compound=true;\n"
    dot *= "\tnode [shape=box, width=1.0, height=1.0, fontsize=9];\n"
    dot *= "\tedge [fontsize=8, arrowhead=onormal];\n"
    for node in nodes
        if typeof(node)==TerminalNode
            dot *= "\t$(object_id(node)) [label=\"$(node.name)\", style=filled, width=0.75, height=0.75]\n"
        else
            if haskey(node_type_symbols, typeof(node))
                dot *= "\t$(object_id(node)) [label=\"$(node_type_symbols[typeof(node)])\\n$(node.name)\"]\n"
            else
                dot *= "\t$(object_id(node)) [label=\"$(typeof(node))\\n$(node.name)\"]\n"
            end
        end
    end
    
    for edge in edges
        tail_id = findfirst(edge.tail.node.interfaces, edge.tail)
        tail_label = "$tail_id $(getName(edge.tail))"
        head_id = findfirst(edge.head.node.interfaces, edge.head)
        head_label = "$head_id $(getName(edge.head))"
        label =  string("FW: ", (edge.tail.message!=nothing) ? "&#9679;" : "&#9675;", " $(edge.tail.message_value_type)\n")
        label *= string("BW: ", (edge.head.message!=nothing) ? "&#9679;" : "&#9675;", " $(edge.head.message_value_type)")
        dot *= "\t$(object_id(edge.tail.node)) -> $(object_id(edge.head.node)) " 
        dot *= "[taillabel=\"$(tail_label)\", headlabel=\"$(head_label)\", label=\"$(label)\"]\n"
    end
    
    dot *= "}";
    
    return dot
end

function graph2dot(composite_node::CompositeNode)
    # Return graph2dot(nodes) where nodes are the internal nodes of composite_node
    nodes = Node[]
    for field in names(composite_node)
        if typeof(getfield(composite_node, field))<:Node
            push!(nodes, getfield(composite_node, field))
        end
    end
    (length(nodes) > 0) || error("CompositeNode does not contain any internal nodes.")

    return graph2dot(nodes)
end
graph2dot(seed_node::Node) = graph2dot(getAllNodes(seed_node, open_composites=false))

function showDot(dot_graph::String, output_file::String="")
    # Show DOT graph
    if output_file!=""
        stdin, proc = writesto(`dot -Tx11 -ofile $(output_file)`)
    else
        stdin, proc = writesto(`dot -Tx11`)
    end
    write(stdin, dot_graph)
    close(stdin)
end

# graphviz() generates a DOT graph and shows it 
graphviz(n::Union(Node, Array{Node,1}); output_file="") = showDot(graph2dot(n), output_file)