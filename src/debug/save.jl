export save, makeDump

function type_class(node::FactorNode)
    type = occursin("placeholder", string(node.id)) ? "data" :
           occursin("clamp", string(node.id)) ? "constant" : "factor"
    class = supertype(typeof(node)) == ForneyLab.DeltaFactor ? "deterministic" : "stochastic"
    (type, class)
end

function makeDump(fg::FactorGraph)
    dump = GraphDump()

    # get nodes
    n_id = collect(keys(fg.nodes))
    # get labels
    n_label = slug.(typeof.(collect(values(fg.nodes))))
    # get types and classes
    n_type_class = type_class.(collect(values(fg.nodes)))

    # get values for const node
    # change clamp and placeholder labels
    for (id, label, type_class) in zip(n_id, n_label, n_type_class)
        value = nothing
        type  = type_class[1]
        class = type_class[2]
        if type == "constant"
            idx   = findlast("clamp_", string(id))[end] + 1
            label = string(id)[idx:end]
            value = string(fg.nodes[id].value)
        elseif type == "data"
            idx   = findlast("placeholder_", string(id))[end] + 1
            label = string(id)[idx:end]
        end
        pushNode!(dump, id, label, type, class, value)
    end

    # get edges
    for fg_edge in fg.edges
        for edge in fg_edge.variable.edges
            id     = string(edge.a.node.id, "_", edge.b.node.id)
            label  = string(edge.variable.id)
            a      = handle(edge.a)
            b      = handle(edge.b)
            source = string(edge.a.node.id)
            target = string(edge.b.node.id)
            pushEdge!(dump, id, label, a, b, source, target)

        end
    end

    return dump
end

function save(dump::GraphDump; output = nothing)
    json_string = JSON.json(dump)

    open(output !== nothing ? output : string(@__DIR__,"/graph.json"), "w") do f
        JSON.print(f, json_string)
    end
end
