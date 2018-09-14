# Functions for visualizing graphs
export draw, drawPdf, drawPng

####################################################
# draw methods
####################################################

function graphviz(dot_graph::AbstractString; external_viewer::Symbol=:None)
    # Show a DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    if external_viewer == :iterm
        viewDotIniTerm(dot_graph)
    elseif external_viewer == :default
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

draw(factor_graph::FactorGraph; schedule=ScheduleEntry[], args...) = graphviz(genDot(nodes(factor_graph), edges(factor_graph), schedule=schedule); args...)
draw(; args...) = draw(currentGraph(); args...)

function drawPng(factor_graph::FactorGraph, filename::AbstractString)
    dot_graph = genDot(nodes(factor_graph), edges(factor_graph))
    png_graph = dot2png(dot_graph)
    png_file = open(filename, "w")
    write(png_file, png_graph)
end
drawPng(filename::AbstractString) = drawPng(currentGraph(), filename)

draw(nodes::Set{FactorNode}; args...) = graphviz(genDot(nodes, edges(nodes)); args...)
draw(nodes::Vector{FactorNode}; args...) = draw(Set(nodes); args...)

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

drawPdf(nodes::Set{FactorNode}, filename::AbstractString) = dot2pdf(genDot(nodes, edges(nodes)), filename)
drawPdf(nodes::Vector{FactorNode}, filename::AbstractString) = drawPdf(Set(nodes), filename)

####################################################
# Internal functions
####################################################

function genDot(nodeset::Set{FactorNode}, edgeset::Set{Edge}; schedule::Schedule=ScheduleEntry[], external_edges::Set{Edge}=Set{Edge}())
    # Return a string representing the graph in DOT format
    # http://en.wikipedia.org/wiki/DOT_(graph_description_language)

    dot = "graph G{splines=true;sep=\"+25,25\";overlap=scalexy;nodesep=1.6;compound=true;\n"
    dot *= "\tnode [shape=box, width=1.0, height=1.0, fontsize=9];\n"
    dot *= "\tedge [fontsize=8, arrowhead=onormal];\n"

    # Draw nodes
    for node in nodeset
        if isa(node, Clamp)
            dot *= "\t$(objectid(node)) [label=\"$(node.id)\", style=filled, width=0.75, height=0.75]\n"
        elseif isa(node, Terminal)
            dot *= "\t$(objectid(node)) [label=\"Terminal $(node.id)\", style=filled, width=0.75, height=0.75]\n"
        else
            dot *= "\t$(objectid(node)) [label=\"$(slug(typeof(node)))\\n$(node.id)\"]\n"
        end
    end

    # Build dictionary for message labels
    msg_labels = Dict{Interface, String}()
    for (i, entry) in enumerate(condense(schedule))
        if entry.msg_update_rule <: SumProductRule
            str = "($i)"
        elseif entry.msg_update_rule <: Union{NaiveVariationalRule, StructuredVariationalRule}
            str = "(($i))"
        elseif entry.msg_update_rule <: ExpectationPropagationRule
            str = "[$i]"
        else
            str = "?$i?"
        end
        msg_labels[entry.interface] = str
    end

    for edge in edgeset
        dot *= edgeDot(edge, msg_labels=msg_labels)
    end

    # Draw external edges
    if length(external_edges) > 0
        # Add invisible nodes to external edges
        external_nodes = setdiff(nodes(external_edges), nodeset)
        for node in external_nodes
            dot *= "\t$(objectid(node)) [label=\"\", shape=none,  width=0.15, height=0.15]\n"
        end

        # Add the external edges themselves
        for edge in external_edges
            dot *= edgeDot(edge, is_external_edge=true)
        end
    end

    dot *= "}";

    return dot
end

function edgeDot(edge::Edge; msg_labels=Dict{Interface, String}(), is_external_edge=false)
    # Generate DOT code for an edge
    dot = ""
    if edge.b != nothing
        b_id = findfirst(isequal(edge.b), edge.b.node.interfaces)
        b_label = "$b_id $(handle(edge.b))"
        b_object_id = objectid(edge.b.node)
    else
        b_object_id = "$(objectid(edge))2"
        dot *= "\t$b_object_id [shape=none, label=\"\", width=0.75, height=0.75]\n"
        b_label = ""
    end

    if edge.a != nothing
        a_id = findfirst(isequal(edge.a), edge.a.node.interfaces)
        a_label = "$a_id $(handle(edge.a))"
        a_object_id = objectid(edge.a.node)
    else
        a_object_id = "$(objectid(edge))1"
        dot *= "\t$a_object_id [shape=none, label=\"\", width=0.75, height=0.75]\n"
        a_label = ""
    end

    dot *= "\t$b_object_id -- $a_object_id "
    if is_external_edge
        dot *= "[taillabel=\"$(b_label)\", headlabel=\"$(a_label)\", style=\"dashed\" color=\"black\"]\n"
    else
        msg_label_a = haskey(msg_labels, edge.a) ? "<B><FONT COLOR=\"blue\">"*msg_labels[edge.a]*"</FONT></B>" : ""
        msg_label_b = haskey(msg_labels, edge.b) ? "<B><FONT COLOR=\"blue\">"*msg_labels[edge.b]*"</FONT></B>" : ""

        props = String[]
        if !isempty(b_label) || !isempty(msg_label_b)
            push!(props, "taillabel=<$(b_label) \n$(msg_label_b)>")
        end
        if !isempty(a_label) || !isempty(msg_label_a)
            push!(props, "headlabel=<$(a_label) \n$(msg_label_a)>")
        end
        push!(props, "label=<<FONT COLOR=\"red\">$(edge.variable.id)</FONT>>")
        if !isempty(props)
            dot *= "[$(join(props,", "))]"
        end
    end
end

function dot2gif(dot_graph::AbstractString)
    # Generate SVG image from DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    proc = open(`dot -Tgif`, "r+")
    write(proc.in, dot_graph)
    close(proc.in)
    return read(proc.out, String)
end

function dot2svg(dot_graph::AbstractString)
    # Generate SVG image from DOT graph
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    proc = open(`dot -Tsvg`, "r+")
    write(proc.in, dot_graph)
    close(proc.in)
    return read(proc.out, String)
end

function dot2png(dot_graph::AbstractString)
    # Generate PNG image from DOT graph
    validateGraphVizInstalled()
    proc = open(`dot -Tpng`, "r+")
    write(proc.in, dot_graph)
    close(proc.in)
    return read(proc.out, String)
end

function validateGraphVizInstalled()
    # Check if GraphViz is installed
    try
        (read(`dot -'?'`, String)[1:10] == "Usage: dot") || error()
    catch
        error("GraphViz is not installed correctly. Make sure GraphViz is installed. If you are on Windows, manually add the path to GraphViz to your path variable. You should be able to run 'dot' from the command line.")
    end
end

viewDotExternal(dot_graph::AbstractString) = (Sys.islinux() ? viewDotExternalInteractive(dot_graph::AbstractString) : viewDotExternalImage(dot_graph::AbstractString))

function viewDotExternalInteractive(dot_graph::AbstractString)
    # View a DOT graph in interactive viewer
    validateGraphVizInstalled() # Show an error if GraphViz is not installed correctly
    proc = open(`dot -Tx11`, "w")
    write(proc.in, dot_graph)
    close(proc.in)
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
