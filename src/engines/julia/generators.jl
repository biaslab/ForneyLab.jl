"""
Generate complete algorithm code
"""
function algorithmString(algo::Dict)
    algo = "begin\n\n"
    for subgraph in algo[:subgraphs]
        algo *= subgraphString(subgraph)
        algo *= "\n\n"
    end
    algo *= "\nend # block"

    return algo
end

"""
Generate algorithm code for a single subgraph
"""
function subgraphString(subgraph::Dict)
    subgraph_str = ""
    if haskey(subgraph, :optimize) && subgraph[:optimize]
        subgraph_str *= optimizeString(subgraph)
        subgraph_str *= "\n\n"
    end

    if haskey(subgraph, :initialization)
        subgraph_str *= initializationString(subgraph)
        subgraph_str *= "\n\n"
    end

    n_messages = length(subgraph[:messages])
    subgraph_str *= "function step$(subgraph[:name])!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, $n_messages))\n\n"
    subgraph_str *= messagesString(subgraph[:messages])
    subgraph_str *= "\n"
    subgraph_str *= marginalsString(subgraph[:marginals])
    subgraph_str *= "return marginals\n\n"
    subgraph_str *= "end"

    return subgraph_str
end

"""
Generate template code for optimize block
"""
function optimizeString(subgraph::Dict)
    name_str = subgraph[:name]

    optim_str =  "# You have created an algorithm that requires updates for (a) clamped parameter(s).\n"
    optim_str *= "# This algorithm requires the definition of a custom `optimize!` function that updates the parameter value(s)\n"
    optim_str *= "# by altering the `data` dictionary in-place. The custom `optimize!` function may be based on the mockup below:\n\n"
    optim_str *= "# function optimize$(name_str)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init$(name_str)())\n"
    optim_str *= "# \t...\n"
    optim_str *= "# \treturn data\n"
    optim_str *= "# end"

    return optim_str
end

"""
Generate code for initialization block
"""
function initializationString(subgraph::Dict)
    init_str = "function init$(subgraph[:name])()\n\n"
    for message in subgraph[:initialization][:messages]
        init_str *= "messages[$(message[:schedule_index])] = Message($(vagueString(message)))\n"
    end
    init_str *= "\nreturn messages\n\n"
    init_str *= "end"

    return init_str
end

"""
Generate code for vague initializations
"""
function vagueString(message::Dict)
    family_str = typeString(message[:family])
    dims = message[:dimensionality]
    if dims == ()
        vague_str = "vague($family_str)"
    else
        vague_str = "vague($family_str, $dims)"
    end

    return vague_str
end

"""
Cast a Type to a String and remove any module prefixes
"""
typeString(type::Type) = split(string(type), '.')[end]

"""
Construct code for messages computation block
"""
function messagesString(messages::Dict)
    messages_str = ""
    for message in messages
        rule_str = typeString(message[:message_update_rule])
        inbounds_str = inboundsString(message[:inbounds])
        messages_str *= "messages[$(message[:schedule_index])] = rule$(rule_str)($inbounds_str)\n"
    end

    return messages_str
end

"""
Concatenate code for separate inbounds to a string
"""
function inboundsString(inbounds::Dict)
    inbounds_str = String[]
    for inbound in inbounds
        push!(inbounds_str, inboundString(inbound))
    end
    return join(inbounds_str, ", ")
end

"""
Generate code for a specific inbound
"""
function inboundString(inbound::Dict)
    if haskey(inbound, :schedule_index) # message inbound
        inbound_str = "messages[$(inbound[:schedule_index])]"
    elseif haskey(inbound, :marginal_id) # marginal inbound
        inbound_str = "marginals[$(inbound[:marginal_id])]"
    elseif haskey(inbound, :buffer_id) # placeholder inbound
        inbound_str = bufferString(inbound)
    elseif haskey(inbound, :value) # value inbound
        inbound_str = valueString(inbound)
    else
        inbound_str = "nothing"
    end

    return inbound_str
end

"""
Generate code for an inbound read from a buffer
"""
function bufferString(inbound::Dict)
    dist_or_msg_str = inbound[:dist_or_msg]
    variate_type_str = inbound[:variate_type]
    buffer_id_str = inbound[:buffer_id]
    inbound_str = "$dist_or_msg_str($variate_type_str, PointMass, m=data[:$buffer_id_str]"
    if haskey(inbound, buffer_index)
        inbound_str *= "[$(inbound[:buffer_index])]"
    end
    inbound_str *= ")"

    return inbound_str
end

"""
Generate code for a clamped inbound
"""
function valueString(inbound::Dict)
    dist_or_msg_str = inbound[:dist_or_msg]
    variate_type_str = inbound[:variate_type]
    buffer_id_str = inbound[:buffer_id]
    inbound_str = "$dist_or_msg_str($variate_type_str, PointMass, m="
    if size(node[:value]) == (1,1)
        inbound_str *= "mat($(node[:value][1]))"
    elseif isa(node[:value], Diagonal)
        inbound_str *= "Diagonal($(node[:value].diag))"
    elseif isa(node[:value], AbstractMatrix) && (size(node[:value])[2] == 1)
        inbound_str *= "[$(join(node[:value], " "))]'"
    else
        inbound_str *= "$(node[:value])"
    end
    inbound_str *= ")"

    return inbound_str
end

"""
Generate code for marginals computation block
"""
function marginalsString(marginals::Dict)
    marginals_str = ""
    for marginal in marginals
        marginals_str *= "marginals[:$(marginals[:marginal_id])] = "
        if marginal[:marginal_update_rule] == Nothing
            inbound = marginals[:inbounds][1]
            marginals_str *= "messages[$(inbound[:schedule_index])].dist"
        elseif marginal[:marginal_update_rule] == Product
            inbound1 = marginals[:inbounds][1]
            inbound2 = marginals[:inbounds][2]
            marginals_str *= "messages[$(inbound1[:schedule_index])].dist * messages[$(inbound2[:schedule_index])].dist"
        else
            rule_str = typeString(marginal[:marginal_update_rule])
            inbounds_str = inboundsString(marginal[:inbounds])
            code *= "rule$(rule_str)($inbounds_str)"
        end
        marginals_str *= "\n"
    end
    isempty(marginals_str) || (marginals_str *= "\n")

    return marginals_str
end

function freeEnergyString(algo::Dict)
    free_energy = algo[:free_energy]

    fe_str  = "function freeEnergy$(free_energy[:name])(data::Dict, marginals::Dict)\n\n"
    fe_str *= "F = 0.0\n\n"
    fe_str *= energiesString(free_energy[:average_energies])
    fe_str *= "\n"
    fe_str *= entropiesString(free_energy[:entropies])
    fe_str *= "\nreturn F\n\n" 
    fe_str *= "end"

    return fe_str
end

function energiesString(average_energies::Dict)
    energies_str = ""
    for energy in average_energies
        node_str = typeString(energy[:node])
        inbounds_str = inboundsString(energy[:inbounds])
        energies_str *= "F += averageEnergy($node_str, $inbounds_str))\n"
    end

    return energies_str
end

function entropiesString(entropies::Dict)
    entropies_str = ""
    for entropy in entropies
        if entropy[:conditional]
            inbounds_str = inboundsString(entropy[:inbounds])
            entropies_str *= "F -= conditionalDifferentialEntropy($inbounds_str)\n"
        else
            entropies_str *= "F -= differentialEntropy(marginals[:$(entropy[:marginal_id])])\n"
        end
    end

    return entropies_str
end