# Functions for visualizing graphs
export draw, drawPdf, drawPng

####################################################
# draw methods
####################################################

function graphviz(dot_graph::AbstractString; external_viewer::Symbol=:None)
    # Show a DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    if external_viewer == :iTerm
        viewDotIniTerm(dot_graph)
    elseif external_viewer == :Default
        viewDotExternal(dot_graph)
    else
        try
            # For iJulia notebook
            #display("text/html", "<img src=\"data:image/gif;base64,"*base64(dot2gif(dot_graph))*"\" />")
            display("text/html", dot2svg(dot_graph))
        catch
            viewDotExternal(dot_graph)
        end
    end
end

draw(factor_graph::FactorGraph; args...) = graphviz(genDot(nodes(factor_graph), edges(factor_graph), wraps=wraps(factor_graph)); args...)
draw(; args...) = draw(currentGraph(); args...)
# draw(composite_node::CompositeNode; args...) = draw(composite_node.internal_graph; args...)

function drawPng(factor_graph::FactorGraph, filename::AbstractString)
    dot_graph = genDot(nodes(factor_graph), edges(factor_graph), wraps=wraps(factor_graph))
    png_graph = dot2png(dot_graph)
    png_file = open(filename, "w")
    write(png_file, png_graph)
end
drawPng(filename::AbstractString) = drawPng(currentGraph(), filename)

draw(nodes::Set{Node}; args...) = graphviz(genDot(nodes, edges(nodes)); args...)
draw(nodes::Vector{Node}; args...) = draw(Set(nodes); args...)

####################################################
# writePdf methods
####################################################

function dot2pdf(dot_graph::AbstractString, filename::AbstractString)
    # Generates a DOT graph and writes it to a pdf file
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    open(`dot -Tpdf -o$(filename)`, "w", STDOUT) do io
        println(io, dot_graph)
    end
end

drawPdf(factor_graph::FactorGraph, filename::AbstractString) = dot2pdf(genDot(nodes(factor_graph), edges(factor_graph)), filename)
drawPdf(filename::AbstractString) = drawPdf(currentGraph(), filename)
# drawPdf(composite_node::CompositeNode, filename::AbstractString) = drawPdf(composite_node.internal_graph, filename)

drawPdf(nodes::Set{Node}, filename::AbstractString) = dot2pdf(genDot(nodes, edges(nodes)), filename)
drawPdf(nodes::Vector{Node}, filename::AbstractString) = drawPdf(Set(nodes), filename)


####################################################
# Internal functions
####################################################

function genDot(nodes::Set{Node}, edges::Set{Edge}; external_edges::Set{Edge}=Set{Edge}(), wraps::Set{Wrap}=Set{Wrap}())
    # Return a string representing the graph in DOT format
    # External edges are edges of which only the head or tail is in the nodes set
    # http://en.wikipedia.org/wiki/DOT_(graph_description_language)
    node_type_symbols = Dict{DataType, AbstractString}(
                            AdditionNode => "+",
                            EqualityNode => "=",
                            ExponentialNode => "exp",
                            GainNode => "GainNode",
                            GainAdditionNode => "GainAdditionNode",
                            GainEqualityNode => "GainEqualityNode",
                            SigmoidNode => "\u03C3"
                        )

    dot = "digraph G{splines=true;sep=\"+25,25\";overlap=scalexy;nodesep=1.6;compound=true;\n"
    dot *= "\tnode [shape=box, width=1.0, height=1.0, fontsize=9];\n"
    dot *= "\tedge [fontsize=8, arrowhead=onormal];\n"
    for node in nodes
        if typeof(node)==TerminalNode
            dot *= "\t$(object_id(node)) [label=\"$(node.id)\", style=filled, width=0.75, height=0.75]\n"
        elseif typeof(node) <: GaussianNode
            dot *= "\t$(object_id(node)) [label=\"N\\n$(node.id)\"]\n"
        else
            if haskey(node_type_symbols, typeof(node))
                dot *= "\t$(object_id(node)) [label=\"$(node_type_symbols[typeof(node)])\\n$(node.id)\"]\n"
            else
                dot *= "\t$(object_id(node)) [label=\"$(typeof(node))\\n$(node.id)\"]\n"
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

    if !isempty(wraps)
        # Add edges to visualize time wraps
        for wrap in wraps
            dot *= "\t$(object_id(wrap.tail)) -> $(object_id(wrap.head)) [style=\"dotted\" color=\"green\"]\n"
        end
    end

    dot *= "}";

    return dot
end

function edgeDot(edge::Edge; is_external_edge=false)
    # Generate DOT code for an edge
    tail_id = findfirst(edge.tail.node.interfaces, edge.tail)
    tail_label = "$tail_id $(handle(edge.tail))"
    head_id = findfirst(edge.head.node.interfaces, edge.head)
    head_label = "$head_id $(handle(edge.head))"
    dot = "\t$(object_id(edge.tail.node)) -> $(object_id(edge.head.node)) "
    if is_external_edge
        dot *= "[taillabel=\"$(tail_label)\", headlabel=\"$(head_label)\", style=\"dashed\" color=\"red\"]\n"
    else
        label = ""
        label *= (typeof(edge.tail.message) <: Message) ? "FW: $(edge.tail.message.payload)" : ""
        label *= (typeof(edge.head.message) <: Message) ? "BW: $(edge.head.message.payload)" : ""
        label *= (typeof(edge.marginal) <: ProbabilityDistribution) ? "Marginal: $(edge.marginal)" : ""
        dot *= "[taillabel=\"$(tail_label)\", headlabel=\"$(head_label)\", label=\"$(label)\" color=\"black\"]"
    end
end

function dot2gif(dot_graph::AbstractString)
    # Generate SVG image from DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    stdout, stdin, proc = readandwrite(`dot -Tgif`)
    write(stdin, dot_graph)
    close(stdin)
    return readstring(stdout)
end

function dot2svg(dot_graph::AbstractString)
    # Generate SVG image from DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    stdout, stdin, proc = readandwrite(`dot -Tsvg`)
    write(stdin, dot_graph)
    close(stdin)
    return readstring(stdout)
end

function dot2png(dot_graph::AbstractString)
    # Generate PNG image from DOT graph
    validateGraphVizInstalled()
    stdout, stdin, proc = readandwrite(`dot -Tpng`)
    write(stdin, dot_graph)
    close(stdin)
    return readstring(stdout)
end

function validateGraphVizInstalled()
    # Check if GraphViz is installed
    try
        (readstring(`dot -?`)[1:10] == "Usage: dot") || error()
    catch
        error("GraphViz is not installed correctly. Make sure GraphViz is installed. If you are on Windows, manually add the path to GraphViz to your path variable. You should be able to run 'dot' from the command line.")
    end
end

viewDotExternal(dot_graph::AbstractString) = (is_linux() ? viewDotExternalInteractive(dot_graph::AbstractString) : viewDotExternalImage(dot_graph::AbstractString))

function viewDotExternalInteractive(dot_graph::AbstractString)
    # View a DOT graph in interactive viewer
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    open(`dot -Tx11`, "w", STDOUT) do io
        println(io, dot_graph)
    end
end

function viewDotExternalImage(dot_graph::AbstractString)
    # Write the image to a file and open it with the default image viewer
    svg = dot2svg(dot_graph)
    filename = tempname()*".svg"
    open(filename, "w") do f
        write(f, svg)
    end
    viewFile(filename)
end

# Based on imgcat script provided by iTerm developers (working in iTerm v.2.9.x)
function viewDotIniTerm(dot_graph::AbstractString)
    png = dot2png(dot_graph)
    base64png = base64encode(png)
    sequence = ("$(Char(0x1b))]1337;File=size=$(length(png));inline=1:$base64png\a\n")
    print(sequence)
end
