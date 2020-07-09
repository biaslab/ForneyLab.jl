export save!

function type_class(node::FactorNode)
    type = occursin("placeholder", string(node.id)) ? "data" :
           occursin("clamp", string(node.id)) ? "constant" : "factor"
    class = supertype(typeof(node)) == ForneyLab.DeltaFactor ? "deterministic" : "stochastic"
    (type, class)
end

function save!(gd::GraphDump, fg::FactorGraph)
    # get nodes
    n_id = collect(keys(fg.nodes))
    # get labels
    n_label = slug.(typeof.(collect(values(fg.nodes))))
    # get types and classes
    n_type_class = type_class.(collect(values(fg.nodes)))
    d_nodes = [Dict("id" => n_id[i], "label" => n_label[i],
                    "type" => n_type_class[i][1],
                    "class" => n_type_class[i][2]) for i in 1:length(n_id)]
    # get interfaces
    for (index, d_id) in enumerate(d_nodes)
        d_nodes[index]["interfaces"] = collect(keys(fg.nodes[d_id["id"]].i))
    end

    # get values for const node
    # change clamp and placeholder labels
    for d_node in d_nodes
        if d_node["type"] == "constant"
            d_node["value"] = string(fg.nodes[d_node["id"]].value)
            idx = findlast("clamp_", string(d_node["id"]))[end] + 1
            d_node["label"] = string(d_node["id"])[idx:end]
        elseif d_node["type"] == "data"
            idx = findlast("placeholder_", string(d_node["id"]))[end] + 1
            d_node["label"] = string(d_node["id"])[idx:end]
        end
    end

    # get edges
    d_edges = Vector{Dict{String, Any}}()
    for fg_edge in fg.edges
        for edge in fg_edge.variable.edges
            push!(d_edges, Dict("source" => string(edge.a.node.id), "target" => string(edge.b.node.id),
                                "a" => handle(edge.a), "b" => handle(edge.b),
                                "label" => string(edge.variable.id),
                                "id" => string(edge.a.node.id, "_", edge.b.node.id)))
        end
    end
    unique!(d_edges)

    data[:nodes] = d_nodes
    data[:edges] = d_edges

    json_string = JSON.json(data)

    open(output !== nothing ? output : string(@__DIR__,"/graph.json"), "w") do f
        JSON.print(f, json_string)
    end
end

#save() = save(currentGraph())
