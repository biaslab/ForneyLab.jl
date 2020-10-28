# Getting started
This page provides the necessary information you need to get started with ForneyLab. We will show the general approach to solving inference problems with ForneyLab by means of a running example: inferring the bias of a coin.

## Installation
Install ForneyLab through the Julia package manager:
```julia
] add ForneyLab
```
!!! note
    If you want to use the graph visualization functions, you need to have [GraphViz](http://www.graphviz.org/) installed.

## Example: Inferring the bias of a coin
The ForneyLab approach to solving inference problems consists of three phases:

1. [Model specification](@ref getting-started-model-specification): ForneyLab offers a domain-specific language to specify your probabilistic model.
2. [Message-passing algorithm generation](@ref): Given a model, ForneyLab contructs an algorithm that infers your quantities of interest.
3. [Message-passing inference execution](@ref): ForneyLab compiles your algorithm to executable (Julia) code to which you may feed data and prior statistics. Executing the algorithm returns (approximate) posterior beliefs for your quantities of inferest.

### Coin flip simulation
Let's start by gathering some data. One approach could be flipping a coin N times and recording each outcome. Here, however, we will simulate this process by sampling some values from a Bernoulli distribution. Each sample can be thought of as the outcome of single flip which is either heads or tails (1 or 0). We will assume that our virtual coin is biased, and lands heads up on 75% of the trials (on average).

```@example 2
N = 25          # number of coin tosses
p = 0.75        # p parameter of the Bernoulli distribution
sbernoulli(n, p) = [(rand() < p) ? 1 : 0 for _ = 1:n] # define Bernoulli sampler
dataset = sbernoulli(N, p); # run N Bernoulli trials
print("dataset = ") ; show(dataset)
```

### [Model specification](@id getting-started-model-specification)
In a Bayesian setting, the next step is to specify our probabilistic model. This amounts to specifying the joint probability of the random variables of the system.

#### Likelihood
We will assume that the outcome of each coin flip is governed by the Bernoulli distribution, i.e.

``y_i \sim Bernoulli(\theta)``,

where ``y_i=1`` represents "heads", ``y_i=0`` represents "tails", and ``θ \in [0,1]`` is the underlying probability of the coin landing heads up for a single coin flip.

#### Prior
We will choose the conjugate prior of the Bernoulli likelihood function defined above, namely the beta distribution, i.e.

``\theta \sim Beta(a, b)``,

where ``a`` and ``b`` are the hyperparameters that encode our prior beliefs about the possible values of ``\theta``. We will assign values to the hyperparameters in a later step.   

#### Joint probability
The joint probability is given by the multiplication of the likelihood and the prior, i.e.

``P(\{y_i\}, θ) = \prod_{i=1}^N P(\{y_i\} | θ) P(θ)``.

Now let's see how to specify this model using ForneyLab's syntax.

```@setup 2
import ForneyLab: draw
using ForneyLab: dot2svg, genDot, FactorGraph, PosteriorFactor

struct SVG
    code :: String
end
Base.show(io::IO, ::MIME"image/svg+xml", b::SVG) = write(io, b.code)

draw(f::FactorGraph) = SVG(dot2svg(genDot(nodes(f), edges(f))))
function draw(rf::PosteriorFactor)
    subgraph_nodes = nodes(rf.internal_edges)
    external_edges = setdiff(edges(subgraph_nodes), rf.internal_edges)
    SVG(dot2svg(genDot(subgraph_nodes, rf.internal_edges, external_edges=external_edges)))
end
```

```@example 2
using ForneyLab
g = FactorGraph()       # create a factor graph
a = placeholder(:a)     # define hyperparameter a as placeholder
b = placeholder(:b)     # define hyperparameter b as placeholder
@RV θ ~ Beta(a, b)      # prior
@RV y ~ Bernoulli(θ)    # likelihood
placeholder(y, :y)      # define y as a placeholder for data
draw(g)                 # draw the factor graph
```
As you can see, ForneyLab offers a model specification syntax that resembles closely to the mathematical equations defined above. Placeholders are used to indicate variables that take specific values at a later date. For example, the way we feed observations into the model is by iteratively assigning each of the observations in our dataset to the random variable `y`. Perhaps less obvious is the fact that the hyperparameters `a` and `b` are also defined as placeholders. The reason is that we will use them to input our current belief about `θ` for every observation that is processed. In section [Message-passing inference execution](@ref) we will see how this is done.

### Message-passing algorithm generation
Once we have defined our model, the next step is to instruct ForneyLab to generate a message-passing algorithm that solves our given inference problem. To do this, we need to specify which type of algorithm we want to use. In this case we will use *belief propagation*, also known as the *sum-product algorithm*. Once we execute the following code, we see that a function called `step!(...)` becomes available in the current scope. This function contains the sum-product message-passing algorithm.
```@example 2
# Generate a message passging sum-product algorithm that infers theta
algo = messagePassingAlgorithm(θ) # derive a sum-product algorithm to infer θ
algo_code = algorithmSourceCode(algo) # convert the algorithm to Julia code
algo_expr = Meta.parse(algo_code) # parse the algorithm into a Julia expression
eval(algo_expr) # evaluate the functions contained in the Julia expression
nothing # hide
```

```julia
:(function step!(data::Dict, marginals::Dict=Dict(), messages::Vector{Message}=Array{Message}(undef, 2))
      #= none:3 =#
      messages[1] = ruleSPBetaOutNPP(nothing, Message(Univariate, PointMass, m=data[:a]), Message(Univariate, PointMass, m=data[:b]))
      #= none:4 =#
      messages[2] = ruleSPBernoulliIn1PN(Message(Univariate, PointMass, m=data[:y]), nothing)
      #= none:6 =#
      marginals[:θ] = (messages[1]).dist * (messages[2]).dist
      #= none:8 =#
      return marginals
  end)
```

### Message-passing inference execution
The last step is to execute the message-passing algorithm. In order to do this, we first need to assign values to the hyperparameters ``a`` and ``b`` which characterize our prior beliefs ``p(\theta)`` about the bias of the coin. Then, we need to feed the observations, one at a time, to the algorithm together with our current belief (prior) ``p(\theta)`` about the bias of the coin. The important thing to note here is that the posterior distribution after processing one observation ``p(\theta|y_{i-1})`` becomes the prior for the processing of the next observation.

```@example 2
# Create a marginals dictionary, and initialize hyperparameters
a = 2.0
b = 7.0
marginals = Dict(:θ => ProbabilityDistribution(Beta, a=a, b=b))

for i in 1:N
    # Feed in datapoints 1 at a time
    data = Dict(:y => dataset[i],
                :a => marginals[:θ].params[:a],
                :b => marginals[:θ].params[:b])

    step!(data, marginals)
end
```

### Results
The plot below shows the result of the inference procedure. We see how the
posterior is a “compromise” between the prior and likelihood, as mandated by Bayesian inference.
```@example 2
using Plots, LaTeXStrings, SpecialFunctions; theme(:default)
plot(fillalpha=0.3, fillrange = 0, leg=false, xlabel=L"\theta", yticks=nothing)
BetaPDF(α, β) = x ->  x^(α-1)*(1-x)^(β-1)/beta(α, β) # beta distribution definition
BernoulliPDF(z, N) = θ -> θ^z*(1-θ)^(N-z) # Bernoulli distribution definition

rθ = range(0, 1, length=100)
p1 = plot(rθ, BetaPDF(a, b), title="Prior", fillalpha=0.3, fillrange = 0, ylabel=L"P(\theta)", c=1,)
p2 = plot(rθ, BernoulliPDF(sum(dataset), N), title="Likelihood", fillalpha=0.3, fillrange = 0, ylabel=L"P(D|\theta)", c=2)
p3 = plot(rθ, BetaPDF(marginals[:θ].params[:a], marginals[:θ].params[:b]), title="Posterior", fillalpha=0.3, fillrange = 0, ylabel=L"P(\theta|D)", c=3)
plot(p1, p2, p3, layout=@layout([a; b; c]))
```

## Where to go next?
There are a set of [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo) available in ForneyLab's repository that demonstrate the more advanced features of ForneyLab. Alternatively, you can head to the [User guide](@ref) which provides more detailed information of how to use ForneyLab to solve inference problems.
