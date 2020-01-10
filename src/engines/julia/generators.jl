export algorithmSourceCode

"""
Generate Julia code for message passing
and optional free energy evaluation
"""
function algorithmSourceCode(algo::Algorithm; free_energy=false)
    algo_str = "begin\n\n"
    for (id, rf) in algo.recognition_factors
        algo_str *= recognitionFactorSourceCode(rf)
        algo_str *= "\n\n"
    end

    if free_energy
        algo_str *= freeEnergySourceCode(algo)
        algo_str *= "\n\n"
    end

    algo_str *= "\nend # block"

    return algo_str
end

"""
Generate Julia code for free energy evaluation
"""
function freeEnergySourceCode(algo::Algorithm)
    fe_str  = "function freeEnergy(data::Dict, marginals::Dict)\n\n"
    fe_str *= "F = 0.0\n\n"
    fe_str *= energiesSourceCode(algo.average_energies)
    fe_str *= "\n"
    fe_str *= entropiesSourceCode(algo.entropies)
    fe_str *= "\nreturn F\n\n" 
    fe_str *= "end"

    return fe_str
end

"""
Generate Julia code for message passing on a single recognition factor
"""
function recognitionFactorSourceCode(rf::RecognitionFactor)
    rf_str = ""
    if isdefined(rf, :optimize) && rf.optimize
        rf_str *= optimizeString(rf)
        rf_str *= "\n\n"
    end

    if isdefined(rf, :initialize) && rf.initialize
        rf_str *= initializationSourceCode(rf)
        rf_str *= "\n\n"
    end

    n_entries = length(rf.schedule)
    rf_str *= "function step$(rf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, $n_entries))\n\n"
    rf_str *= scheduleSourceCode(rf.schedule)
    rf_str *= "\n"
    rf_str *= marginalTableSourceCode(rf.marginal_table)
    rf_str *= "return marginals\n\n"
    rf_str *= "end"

    return rf_str
end

"""
Generate template code for optimize block
"""
function optimizeString(rf::RecognitionFactor)
    optim_str =  "# You have created an algorithm that requires updates for (a) clamped parameter(s).\n"
    optim_str *= "# This algorithm requires the definition of a custom `optimize!` function that updates the parameter value(s)\n"
    optim_str *= "# by altering the `data` dictionary in-place. The custom `optimize!` function may be based on the mockup below:\n\n"
    optim_str *= "# function optimize$(rf.id)!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=init$(rf.id)())\n"
    optim_str *= "# \t...\n"
    optim_str *= "# \treturn data\n"
    optim_str *= "# end"

    return optim_str
end

"""
Generate code for initialization block (if required)
"""
function initializationSourceCode(rf::RecognitionFactor)
    init_str = "function init$(rf.id)()\n\n"
    for entry in rf.schedule
        if isdefined(entry, :initialize) && entry.initialize
            init_str *= "messages[$(entry.schedule_index)] = Message($(vagueSourceCode(entry)))\n"
        end
    end
    init_str *= "\nreturn messages\n\n"
    init_str *= "end"

    return init_str
end

"""
Generate code for evaluating the average energy
"""
function energiesSourceCode(average_energies::Vector)
    energies_str = ""
    for energy in average_energies
        node_str = removePrefix(energy[:node])
        inbounds_str = inboundsString(energy[:inbounds])
        energies_str *= "F += averageEnergy($node_str, $inbounds_str))\n"
    end

    return energies_str
end

"""
Generate code for evaluating the entropy
"""
function entropiesSourceCode(entropies::Vector)
    entropies_str = ""
    for entropy in entropies
        inbounds_str = inboundsString(entropy[:inbounds])
        if haskey(entropy, :conditional) && entropy[:conditional]
            entropies_str *= "F -= conditionalDifferentialEntropy($inbounds_str)\n"
        else
            entropies_str *= "F -= differentialEntropy($inbounds_str)\n"
        end
    end

    return entropies_str
end

"""
Construct code for message updates
"""
function scheduleSourceCode(schedule::Schedule)
    schedule_str = ""
    for entry in schedule
        rule_str = removePrefix(entry.message_update_rule)
        inbounds_str = inboundsString(entry.inbounds)
        schedule_str *= "messages[$(entry.schedule_index)] = rule$(rule_str)($inbounds_str)\n"
    end

    return schedule_str
end

"""
Generate code for marginal updates
"""
function marginalTableSourceCode(table::MarginalTable)
    table_str = ""
    for entry in table
        table_str *= "marginals[:$(entry.marginal_id)] = "
        if entry.marginal_update_rule == Nothing
            inbound = entry.inbounds[1]
            table_str *= "messages[$(inbound.schedule_index)].dist"
        elseif entry.marginal_update_rule == Product
            inbound1 = entry.inbounds[1]
            inbound2 = entry.inbounds[2]
            table_str *= "messages[$(inbound1.schedule_index)].dist * messages[$(inbound2.schedule_index)].dist"
        else
            rule_str = removePrefix(entry.marginal_update_rule)
            inbounds_str = inboundsString(entry.inbounds)
            table_str *= "rule$(rule_str)($inbounds_str)"
        end
        table_str *= "\n"
    end
    isempty(table_str) || (table_str *= "\n")

    return table_str
end

"""
Generate code for vague initializations
"""
function vagueSourceCode(entry::ScheduleEntry)
    family_str = removePrefix(entry.family)
    dims = entry.dimensionality
    if dims == ()
        vague_str = "vague($family_str)"
    else
        vague_str = "vague($family_str, $dims)"
    end

    return vague_str
end

"""
Generate code for message/marginal/energy/entropy computation inbounds
"""
function inboundsString(inbounds::Vector)
    inbounds_str = String[]
    for inbound in inbounds
        push!(inbounds_str, inboundSourceCode(inbound))
    end
    return join(inbounds_str, ", ")
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

    inbound_str = ""
    for (key, val) in inbound
        if key != :keyword
            if keyword_flag
                inbound_str = "$(removePrefix(key))=$(removePrefix(val))"
            else
                inbound_str = removePrefix(val)
            end
        end
    end

    return inbound_str
end
function inboundSourceCode(inbound::Clamp{V}) where V<:VariateType # Buffer or value inbound
    dist_or_msg_str = removePrefix(inbound.dist_or_msg)
    variate_type_str = removePrefix(V)

    inbound_str = "$dist_or_msg_str($variate_type_str, PointMass, m="
    if isdefined(inbound, :buffer_id)
        # Inbound is read from buffer
        inbound_str *= "data[:$(inbound.buffer_id)]"
        if isdefined(inbound, :buffer_index) && (inbound.buffer_index > 0)
            inbound_str *= "[$(inbound.buffer_index)]"
        end
    else
        # Inbound is read from clamp node value
        inbound_str *= valueSourceCode(inbound.value)
    end
    inbound_str *= ")"

    return inbound_str
end

"""
Convert a value to parseable Julia code
"""
valueSourceCode(val::Union{Vector, Number}) = string(val)
valueSourceCode(val::Diagonal) = "Diagonal($(val.diag))"
function valueSourceCode(val::AbstractMatrix)
    if size(val) == (1,1)
        val_str = "mat($(val[1]))"
    else
        val_str = string(val)
    end

    return val_str
end

"""
Remove module prefixes from types and functions
"""
removePrefix(arg::Any) = arg # Do not remove prefix in general (might otherwise e.g. split floats)
removePrefix(type::Type) = split(string(type), '.')[end]
removePrefix(func::Function) = split(string(func), '.')[end]
