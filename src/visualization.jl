# Functions for visualizing graphs
export graphViz, graphPdf

####################################################
# graphViz methods
####################################################

function graphViz(dot_graph::String; external_viewer::Bool=false)
    # Show a DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    if external_viewer
        viewDotExternal(dot_graph)
    else
        try
            # For iJulia notebook
            display("image/svg+xml", dot2svg(dot_graph))
        catch
            viewDotExternal(dot_graph)
        end
    end
end

graphViz(factor_graph::FactorGraph; args...) = graphViz(genDot(getNodes(factor_graph, open_composites=false), getEdges(factor_graph)); args...)
graphViz(; args...) = graphViz(getCurrentGraph(); args...)

graphViz(subgraph::Subgraph; args...) = graphViz(genDot(subgraph.nodes, subgraph.internal_edges, subgraph.external_edges); args...)

function graphViz(composite_node::CompositeNode; args...)
    # Visualize the internals of composite_node
    nodes = Set{Node}()
    for field in names(composite_node)
        if typeof(getfield(composite_node, field)) <: Node
            push!(nodes, getfield(composite_node, field))
        end
    end
    (length(nodes) > 0) || error("CompositeNode does not contain any internal nodes.")

    return graphViz(genDot(nodes, getEdges(nodes, include_external=false)); args...)
end

graphViz(nodes::Set{Node}; args...) = graphViz(genDot(nodes, getEdges(nodes, include_external=false)); args...)
graphViz(nodes::Vector{Node}; args...) = graphViz(Set(nodes); args...)

####################################################
# graphPdf methods
####################################################

function graphPdf(dot_graph::String, filename::String)
    # Generates a DOT graph and writes it to a pdf file
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    open(`dot -Tpdf -o$(filename)`, "w", STDOUT) do io
        println(io, dot_graph)
    end
end

graphPdf(factor_graph::FactorGraph, filename::String) = graphPdf(genDot(getNodes(factor_graph, open_composites=false), getEdges(factor_graph)), filename)
graphPdf(filename::String) = graphPdf(getCurrentGraph(), filename)

graphPdf(subgraph::Subgraph, filename::String) = graphPdf(genDot(subgraph.nodes, subgraph.internal_edges, subgraph.external_edges), filename)

function graphPdf(composite_node::CompositeNode, filename::String)
    # Visualize the internals of composite_node
    nodes = Set{Node}()
    for field in names(composite_node)
        if typeof(getfield(composite_node, field)) <: Node
            push!(nodes, getfield(composite_node, field))
        end
    end
    (length(nodes) > 0) || error("CompositeNode does not contain any internal nodes.")

    return graphPdf(genDot(nodes, getEdges(nodes, include_external=false)), filename)
end

graphPdf(nodes::Set{Node}, filename::String) = graphPdf(genDot(nodes, getEdges(nodes, include_external=false)), filename)
graphPdf(nodes::Vector{Node}, filename::String) = graphPdf(Set(nodes), filename)


####################################################
# Internal functions
####################################################

function genDot(nodes::Set{Node}, edges::Set{Edge}, external_edges::Set{Edge}=Set{Edge}())
    # Return a string representing the graph in DOT format
    # External edges are edges of which only the head or tail is in the nodes set
    # http://en.wikipedia.org/wiki/DOT_(graph_description_language)
    node_type_symbols = {   AdditionNode => "+",
                            EqualityNode => "="}
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
        dot *= edgeDot(edge)
    end

    if length(external_edges) > 0
        # Add nodes connected to the external edges
        external_nodes = Set{Node}()
        for external_edge in external_edges
            for node in [external_edge.tail.node, external_edge.head.node]
                if !(node in nodes) push!(external_nodes, node) end
            end
        end
        for node in external_nodes
            dot *= "\t$(object_id(node)) [label=\"\", shape=none,  width=0.15, height=0.15]\n"
        end

        # Add the external edges themselves
        for edge in external_edges
            dot *= edgeDot(edge, is_external_edge=true)
        end
    end

    dot *= "}";
    
    return dot
end

function edgeDot(edge::Edge; is_external_edge=false)
    # Generate DOT code for an edge
    tail_id = findfirst(edge.tail.node.interfaces, edge.tail)
    tail_label = "$tail_id $(getName(edge.tail))"
    head_id = findfirst(edge.head.node.interfaces, edge.head)
    head_label = "$head_id $(getName(edge.head))"
    dot = "\t$(object_id(edge.tail.node)) -> $(object_id(edge.head.node)) " 
    if is_external_edge
        dot *= "[taillabel=\"$(tail_label)\", headlabel=\"$(head_label)\", style=\"dashed\" color=\"red\"]\n"
    else
        label =  string("FW: ", (edge.tail.message!=nothing) ? "&#9679;" : "&#9675;", " $(edge.tail.message_payload_type)\n")
        label *= string("BW: ", (edge.head.message!=nothing) ? "&#9679;" : "&#9675;", " $(edge.head.message_payload_type)\n")
        dot *= "[taillabel=\"$(tail_label)\", headlabel=\"$(head_label)\", label=\"$(label)\" color=\"black\"]\n"
    end
end

function dot2svg(dot_graph::String)
    # Generate SVG image from DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    stdout, stdin, proc = readandwrite(`dot -Tsvg`)
    write(stdin, dot_graph)
    close(stdin)
    return readall(stdout)
end

function validateGraphVizInstalled()
    # Check if GraphViz is installed
    try
        (readall(`dot -?`)[1:10] == "Usage: dot") || error()
    catch
        error("GraphViz is not installed correctly. Make sure GraphViz is installed. If you are on Windows, manually add the path to GraphViz to your path variable. You should be able to run 'dot' from the command line.")
    end
end

viewDotExternal(dot_graph::String) = (@windows? viewDotExternalImage(dot_graph::String) : viewDotExternalInteractive(dot_graph::String))

function viewDotExternalInteractive(dot_graph::String)
    # View a DOT graph in interactive viewer
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    open(`dot -Tx11`, "w", STDOUT) do io
        println(io, dot_graph)
    end
end

function viewDotExternalImage(dot_graph::String)
    # Write the image to a file and open it with the default image viewer
    svg = dot2svg(dot_graph)
    filename = tempname()*".svg"
    open(filename, "w") do f
        write(f, svg)
    end
    viewFile(filename)
end