#####################
# Unit tests
#####################

facts("FactorGraph unit tests") do
    context("FactorGraph() should initialize a factor graph") do
        fg = FactorGraph()
        @fact typeof(fg) => FactorGraph
        @fact length(fg.inference_schemes) => 1
        @fact typeof(fg.active_scheme) => InferenceScheme
        @fact fg.active_scheme => fg.inference_schemes[1]
        @fact fg.nodes => Set{Node}()
        @fact fg.edges => Set{Edge}()
        @fact current_graph => fg # Global should be set
    end

    context("currentGraph() should return a current graph object") do
        FactorGraph() # Reset
        @fact currentGraph() => current_graph
        my_graph = currentGraph()  # Make local pointer to global variable
        @fact typeof(my_graph) => FactorGraph
    end

    context("setCurrentGraph() should set a new current graph object") do
        my_first_graph = FactorGraph() # Reset
        my_second_graph = FactorGraph()
        @fact my_first_graph == current_graph => false
        setCurrentGraph(my_first_graph)
        @fact my_first_graph == current_graph => true
    end
end