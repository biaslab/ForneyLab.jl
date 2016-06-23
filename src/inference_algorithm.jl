Base.deepcopy(::InferenceAlgorithm) = error("deepcopy(::InferenceAlgorithm) is not possible. You should construct a new InferenceAlgorithm.")

function compile!{T<:ScheduleEntry}(schedule::Vector{T}, algorithm::InferenceAlgorithm)
    # Generate entry.execute for every entry in schedule
    for schedule_entry in schedule
        compile!(schedule_entry, algorithm)
    end

    return schedule
end