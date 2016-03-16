Base.deepcopy(::InferenceAlgorithm) = error("deepcopy(::InferenceAlgorithm) is not possible. You should construct a new InferenceAlgorithm.")

function compile!(schedule::Schedule, algorithm::InferenceAlgorithm)
    # Generate entry.execute for every entry in schedule
    for schedule_entry in schedule
        compile!(schedule_entry, algorithm)
    end

    return schedule
end

function compile!(schedule_entry::ScheduleEntry, algorithm::InferenceAlgorithm)
    return compile!(schedule_entry, Val{symbol(schedule_entry.rule)}, algorithm)
end