function generateSchedule(variable::Variable;
                    graph::FactorGraph=currentGraph())

    # Generate schedule
    dg = summaryDependencyGraph(graph)
    iface_list = children([edge.head, edge.tail], dg)

    # Infer message types
    schedule = inferTypes(iface_list, SumProductRule)

....

    return schedule

end