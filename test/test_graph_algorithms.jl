facts("Graph algorithms") do
    context("children(...) should find all children of a vertex and return a topologically sorted array.") do
        fg = FactorGraph()
        Edge(MockNode(1; id=:n1).interfaces[1], MockNode(2; id=:n2).interfaces[1])
        Edge(n(:n2).interfaces[2], MockNode(3; id=:n3).interfaces[1])
        Edge(MockNode(1; id=:n4).interfaces[1], n(:n3).interfaces[2])
        Edge(n(:n3).interfaces[3], MockNode(1; id=:n5).interfaces[1])
        dg = summaryDependencyGraph(fg)

        # No special filters
        result = children(n(:n3).interfaces[2], dg)
        @fact Set{Interface}(result) --> Set{Interface}([   n(:n3).interfaces[2];
                                                            n(:n5).interfaces[1];
                                                            n(:n2).interfaces[2];
                                                            n(:n1).interfaces[1]])
        @fact findfirst(result, n(:n1).interfaces[1]) --> less_than(findfirst(result, n(:n2).interfaces[2]))
        @fact findfirst(result, n(:n3).interfaces[2]) --> length(result)

        # With breakers
        result = children(n(:n3).interfaces[2], dg, breakers=Set([n(:n1).interfaces[1]]))
        @fact Set{Interface}(result) --> Set{Interface}([   n(:n3).interfaces[2];
                                                            n(:n5).interfaces[1];
                                                            n(:n2).interfaces[2]])

        # With restricted vertex set
        result = children([n(:n3).interfaces[2]], dg, restrict_to=Set([   n(:n3).interfaces[2];
                                                                        n(:n5).interfaces[1];
                                                                        n(:n2).interfaces[2];
                                                                        n(:n4).interfaces[1]]))
        @fact Set{Interface}(result) --> Set{Interface}([   n(:n3).interfaces[2];
                                                            n(:n5).interfaces[1];
                                                            n(:n2).interfaces[2]])

        # Loopy dependency graph
        initializeLoopyGraph()
        fg = currentGraph()
        dg = summaryDependencyGraph(fg)
        @fact_throws ArgumentError children(n(:inhibitor).i[:in], dg)
        result = children(n(:inhibitor).i[:in], dg, allow_cycles=true)
        @fact Set{Interface}(result) --> Set{Interface}([   n(:noise).i[:out];
                                                            n(:add).i[:in1];
                                                            n(:driver).i[:in];
                                                            n(:inhibitor).i[:in]])
    end
end
