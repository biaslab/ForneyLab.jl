export algorithmSourceCode

"""
Generate Julia code for message passing
and optional free energy evaluation
"""
function algorithmSourceCode(algo::InferenceAlgorithm; free_energy=false, debug::Bool=false)
    algo_code = "begin\n\n"
    for (_, pf) in algo.posterior_factorization
        algo_code *= posteriorFactorSourceCode(pf; debug = debug)
        algo_code *= "\n\n"
    end

    if free_energy
        algo_code *= freeEnergySourceCode(algo; debug = debug)
        algo_code *= "\n\n"
    end

    algo_code *= "end # block"

    return algo_code
end

"""
Generate Julia code for free energy evaluation
"""
function freeEnergySourceCode(algo::InferenceAlgorithm; debug::Bool=false)
    algo.posterior_factorization.free_energy_flag || error("Required quantities for free energy evaluation are not computed by the algorithm. Make sure to flag free_energy=true upon algorithm construction to schedule computation of required quantities.")

    fe_code  = "function freeEnergy$(algo.id)(data::Dict, marginals::Dict; dump::Union{GraphDump, Nothing}=nothing)\n\n"
    fe_code *= "F = 0.0\n\n"
    fe_code *= energiesSourceCode(algo.average_energies, debug = debug)
    fe_code *= "\n"
    fe_code *= entropiesSourceCode(algo.entropies, debug = debug)
    fe_code *= "\nif dump !== nothing"
    fe_code *= "\n    push!(dump.steps, AlgorithmStep())"
    fe_code *= "\nend"
    fe_code *= "\nreturn F\n\n"
    fe_code *= "end"

    return fe_code
end

"""
Generate Julia code for message passing on a single posterior factor
"""
function posteriorFactorSourceCode(pf::PosteriorFactor; debug::Bool = false)
    pf_code = ""
    if pf.optimize
        pf_code *= optimizeSourceCode(pf)
        pf_code *= "\n\n"
    end

    if pf.initialize
        pf_code *= initializationSourceCode(pf)
        pf_code *= "\n\n"
    end

    n_entries = length(pf.schedule)
    pf_code *= "function step$(pf.algorithm_id)$(pf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, $n_entries); dump::Union{GraphDump, Nothing} = nothing)\n\n"
    if debug
        pf_code *= "\npush!(dump.steps[end].messages, Vector{MessageData}())\n"
        pf_code *= "\npush!(dump.steps[end].marginals, Vector{MarginalData}())\n"
    end
    pf_code *= scheduleSourceCode(pf.schedule; debug = debug)
    pf_code *= "\n"
    pf_code *= marginalTableSourceCode(pf.marginal_table; debug = debug)
    pf_code *= "return marginals\n\n"
    pf_code *= "end"

    return pf_code
end

"""
Generate template code for optimize block
"""
function optimizeSourceCode(pf::PosteriorFactor)
    optim_code =  "# You have created an algorithm that requires updates for (a) clamped parameter(s).\n"
    optim_code *= "# This algorithm requires the definition of a custom `optimize!` function that updates the parameter value(s)\n"
    optim_code *= "# by altering the `data` dictionary in-place. The custom `optimize!` function may be based on the mockup below:\n\n"
    optim_code *= "# function optimize$(pf.algorithm_id)$(pf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init$(pf.algorithm_id)$(pf.id)())\n"
    optim_code *= "# \t...\n"
    optim_code *= "# \treturn data\n"
    optim_code *= "# end"

    return optim_code
end

"""
Generate code for initialization block (if required)
"""
function initializationSourceCode(pf::PosteriorFactor)
    init_code = "function init$(pf.algorithm_id)$(pf.id)()\n\n"

    n_messages = length(pf.schedule)

    init_code *= "messages = Array{Message}(undef, $n_messages)\n\n"

    for entry in pf.schedule
        if entry.initialize
            init_code *= "messages[$(entry.schedule_index)] = Message($(vagueSourceCode(entry)))\n"
        end
    end
    init_code *= "\nreturn messages\n\n"
    init_code *= "end"

    return init_code
end

"""
Generate code for evaluating the average energy
"""
function energiesSourceCode(average_energies::Vector; debug::Bool=false)
    energies_code = ""
    for energy in average_energies
        count_code = countingNumberSourceCode(energy[:counting_number])
        node_code = removePrefix(energy[:node])
        inbounds_code = inboundsSourceCode(energy[:inbounds])

        if debug
            energies_code *= "v = $(count_code)averageEnergy($node_code, $inbounds_code)\n"
            energies_code *= "pushScore!(dump, \"$(energy[:node_id])\", v,  \"average_energy\")\n"
            energies_code *= "F += v\n"
        else
            energies_code *= "F += $(count_code)averageEnergy($node_code, $inbounds_code)\n"
        end
    end

    return energies_code
end

"""
Generate code for evaluating the entropy
"""
function entropiesSourceCode(entropies::Vector; debug::Bool=false)
    entropies_code = ""
    for entropy in entropies
        count_code = countingNumberSourceCode(entropy[:counting_number])
        inbound_code = inboundSourceCode(entropy[:inbound])

        if debug
            entropies_code *= "v = $(count_code)differentialEntropy($inbound_code)\n"
            for edge in entropy[:target].edges
                edge_id = string(edge.a.node.id, "_", edge.b.node.id)
                entropies_code *= "pushScore!(dump, \"$(edge_id)\", v, \"entropy\")\n"
            end
            entropies_code *= "F -= v\n"
        else
            entropies_code *= "F -= $(count_code)differentialEntropy($inbound_code)\n"
        end
    end

    return entropies_code
end

"""
Generate code for counting number
"""
function countingNumberSourceCode(cnt::Int64)
    if cnt == 1
        count_code = ""
    else
        count_code = "$(cnt)*"
    end

    return count_code
end

"""
Construct code for message updates
"""
function scheduleSourceCode(schedule::Schedule; debug::Bool = false)
    schedule_code = ""
    for entry in schedule
        rule_code = removePrefix(entry.message_update_rule)
        inbounds_code = inboundsSourceCode(entry.inbounds)
        schedule_code *= "messages[$(entry.schedule_index)] = rule$(rule_code)($inbounds_code)\n"

        if debug
            edge   = entry.interface.edge
            edgeid = string(edge.a.node.id, "_", edge.b.node.id)
            type   = entry.interface.node.id === edge.a.node.id ? "forward" : "backward"
            schedule_code *= "\npushMessage!(dump, \"$(edgeid)\", \"$(type)\", messages[$(entry.schedule_index)])\n"
        end
    end

    return schedule_code
end

"""
Generate code for marginal updates
"""
function marginalTableSourceCode(table::MarginalTable; debug::Bool = false)
    table_code = ""
    for entry in table
        table_code *= "marginals[:$(entry.marginal_id)] = "
        if entry.marginal_update_rule == Nothing
            inbound = entry.inbounds[1]
            table_code *= "messages[$(inbound.schedule_index)].dist"
        elseif entry.marginal_update_rule == Product
            inbound1 = entry.inbounds[1]
            inbound2 = entry.inbounds[2]
            table_code *= "messages[$(inbound1.schedule_index)].dist * messages[$(inbound2.schedule_index)].dist"
        else
            rule_code = removePrefix(entry.marginal_update_rule)
            inbounds_code = inboundsSourceCode(entry.inbounds)
            table_code *= "rule$(rule_code)($inbounds_code)"
        end

        if debug
            edgeIDs = map(e -> string(e.a.node.id, "_", e.b.node.id), entry.target.edges)
            table_code *= "\npushMarginal!(dump, \"$(entry.marginal_id)\", $(edgeIDs), marginals[:$(entry.marginal_id)])\n"
        end

        table_code *= "\n"
    end
    isempty(table_code) || (table_code *= "\n")

    return table_code
end

"""
Generate code for vague initializations
"""
function vagueSourceCode(entry::ScheduleEntry)
    family_code = removePrefix(entry.family)
    dims = entry.dimensionality
    if dims == ()
        vague_code = "vague($family_code)"
    else
        vague_code = "vague($family_code, $dims)"
    end

    return vague_code
end

"""
Generate code for message/marginal/energy/entropy computation inbounds
"""
function inboundsSourceCode(inbounds::Vector)
    inbounds_code = String[]
    for inbound in inbounds
        push!(inbounds_code, inboundSourceCode(inbound))
    end
    return join(inbounds_code, ", ")
end

"""
Generate code for a single inbound (overloaded for specific inbound type)
"""
inboundSourceCode(inbound::Nothing) = "nothing"
inboundSourceCode(inbound::ScheduleEntry) = "messages[$(inbound.schedule_index)]"
inboundSourceCode(inbound::MarginalEntry) = "marginals[:$(inbound.marginal_id)]"
function inboundSourceCode(inbound::Dict) # Custom inbound
    keyword_flag = true # Default includes keyword in custom argument
    if haskey(inbound, :keyword)
        keyword_flag = inbound[:keyword]
    end

    inbound_code = ""
    for (key, val) in inbound
        if key != :keyword
            if keyword_flag
                inbound_code = "$(removePrefix(key))=$(removePrefix(val))"
            else
                inbound_code = removePrefix(val)
            end
        end
    end

    return inbound_code
end
function inboundSourceCode(inbound::Clamp{V}) where V<:VariateType # Buffer or value inbound
    dist_or_msg_code = removePrefix(inbound.dist_or_msg)
    variate_type_code = removePrefix(V)

    inbound_code = "$dist_or_msg_code($variate_type_code, PointMass, m="
    if isdefined(inbound, :buffer_id)
        # Inbound is read from buffer
        inbound_code *= "data[:$(inbound.buffer_id)]"
        if isdefined(inbound, :buffer_index) && (inbound.buffer_index > 0)
            inbound_code *= "[$(inbound.buffer_index)]"
        end
    else
        # Inbound is read from clamp node value
        inbound_code *= valueSourceCode(inbound.value)
    end
    inbound_code *= ")"

    return inbound_code
end

inboundSourceCode(inbound::Number) = string(inbound)

"""
Convert a value to parseable Julia code
"""
valueSourceCode(val::Union{Vector, Number}) = string(val)
valueSourceCode(val::Diagonal) = "Diagonal($(val.diag))"
function valueSourceCode(val::AbstractMatrix)
    # If the dimensionality in at least one direction is 1, a matrix needs to be
    # constructed explicitly; otherwise the output Julia code will construct a vector.
    d = size(val)
    if d == (1,1)
        val_code = "mat($(val[1]))" # Shorthand notation for a (1,1) matrix reshape
    elseif 1 in d
        val_code = "reshape($(vec(val)), $d)" # Explicit reshape
    else
        val_code = string(val)
    end

    return val_code
end

"""
Remove module prefixes from types and functions
"""
removePrefix(arg::Any) = arg # Do not remove prefix in general
removePrefix(tup::Tuple) = string(tup)
removePrefix(num::Number) = string(num)
removePrefix(vect::Vector) = string(vect)
removePrefix(type::Type) = split(string(type), '.')[end]
removePrefix(func::Function) = split(string(func), '.')[end]
