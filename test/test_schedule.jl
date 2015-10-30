facts("General Schedule and ScheduleEntry tests") do
    FactorGraph()

    # ScheduleEntry
    schedule_entry1 = ScheduleEntry(MockNode(id=:mock1).i[:out])
    schedule_entry2 = ScheduleEntry(MockNode(id=:mock2).i[:out])
    @fact is(schedule_entry1.interface, n(:mock1).i[:out]) --> true
    @fact_throws deepcopy(schedule_entry1)

    # copy(::ScheduleEntry)
    schedule_entry1_copy = copy(schedule_entry1)
    @fact is(schedule_entry1, schedule_entry1_copy) --> false
    @fact is(schedule_entry1.interface, schedule_entry1_copy.interface) --> true
    @fact is(schedule_entry1.message_calculation_rule, schedule_entry1_copy.message_calculation_rule) --> true
    @fact isdefined(schedule_entry1_copy, :post_processing) --> isdefined(schedule_entry1, :post_processing)
    if isdefined(schedule_entry1_copy, :post_processing)
        @fact is(schedule_entry1_copy.post_processing, schedule_entry1.post_processing) --> true
    end

    # Schedule
    schedule = [schedule_entry1, schedule_entry2]
    @fact typeof(schedule) --> Schedule
    schedule_deepcopy = deepcopy(schedule)
    for idx=1:length(schedule)
        @fact is(schedule[idx], schedule_deepcopy[idx]) --> false
        @fact is(schedule[idx].interface, schedule_deepcopy[idx].interface) --> true
    end
end
