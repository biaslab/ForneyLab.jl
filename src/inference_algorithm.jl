export InferenceAlgorithm, GenericAlgorithm

abstract InferenceAlgorithm

type GenericAlgorithm <: InferenceAlgorithm
    execute::Function
end

Base.deepcopy(::InferenceAlgorithm) = error("deepcopy(::InferenceAlgorithm) is not possible. You should construct a new InferenceAlgorithm.")

function compile!(schedule::Schedule, algorithm::InferenceAlgorithm)
    # Generate execute functions for all schedule entries
    for schedule_entry in schedule
        compile!(schedule_entry, algorithm)
    end

    return schedule
end

compile!(schedule_entry::ScheduleEntry, algorithm::InferenceAlgorithm) = compile!(schedule_entry, Val{symbol(schedule_entry.message_calculation_rule)}, algorithm)