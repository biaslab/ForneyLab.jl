@doc "Factor node. Subtypes: $(format(subtypes(Node)))" Node
@doc "Probability distribution. Subtypes: $(format(subtypes(ProbabilityDistribution)))" ProbabilityDistribution
@doc "Univariate probability distribution. Subtypes: $(format(subtypes(UnivariateProbabilityDistribution)))" UnivariateProbabilityDistribution
@doc "Multivariate probability distribution. Subtypes: $(format(subtypes(MultivariateProbabilityDistribution)))" MultivariateProbabilityDistribution
@doc "Inference algorithm. Subtypes: $(format(subtypes(InferenceAlgorithm)))" InferenceAlgorithm
@doc "Sum-product algorithm. Subtypes: $(format(subtypes(AbstractSumProduct)))" AbstractSumProduct

function ruleMethodDocstring(rule::Function, sig)
    # Return docstring for rule with signature sig
    for m in Docs.modules
        if haskey(Base.Docs.meta(m), rule)
            docs = Base.Docs.meta(m)[rule]
            if isa(docs, Base.Docs.FuncDoc)
                if sig in keys(docs.meta)
                    return docs.meta[sig]
                end
            end
        end
    end

    return Base.Markdown.parse("**Undocumented**")
end


"""
rules(node_type, [rule], [outbound=3])

Display all message calculation rules for node_type.
Optional argument rule can specify a rule, like sumProductRule!.
Optionally, you can filter on the outbound interface id.

    rules(GaussianNode)
    rules(AdditionNode, sumProductRule!)
    rules(AdditionNode, outbound=2)
"""
function rules(node_type::DataType, rule::Function; outbound::Int=0)
    (node_type <: Node) || error("node_type should be a subtype of Node")

    results = []
    for m in methods(rule)
        if m.sig.parameters[1] <: node_type
            if outbound==0 || Type{Val{outbound}}==m.sig.parameters[2]
                docstring = ruleMethodDocstring(rule, m.sig)
                push!(results, (m.sig, docstring))
            end
        end
    end

    if !isempty(results)
        sort!(results, lt = (a, b) -> Base.Docs.type_morespecific(first(a), first(b)))
        for (msig, docstring) in results
            println("\n", replace(string(rule), "ForneyLab.", ""))
            println(replace(string(msig), "ForneyLab.", ""))
            display(docstring)
        end
    end

end

function rules(node_type::DataType; outbound::Int=0)
    (node_type <: Node) || error("node_type should be a subtype of Node")

    for rule in [sumProductRule!; variationalRule!; expectationRule!]
        rules(node_type, rule, outbound=outbound)
    end
end