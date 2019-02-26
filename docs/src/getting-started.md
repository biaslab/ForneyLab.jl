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

1. [Model specification](@ref getting-started-model-specification): ForneyLab provides a simple meta-language to specify probabilistic models.
2. [Message-passing algorithm generation](@ref): Given a model, it is ForneyLab's job to provide you with an algorithm that runs inference on your quantities of interest.
3. [Message-passing inference execution](@ref): Feed observations and a prior belief of your quantities of interest to the algorithm and you will get an updated posterior in return.

### Coin flip simulation
Let's start by gathering some data. One approach could be flipping a coin N times and recording each outcome. Here, however, we will simulate this process by sampling a random variable with a Bernoulli distribution. Each sample can be thought of as the outcome of single flip of a coin which is either heads or tails (1 or 0). We will assume that our virtual coin has an underlying probability of 75% of landing heads up.

```@example 1
N = 25          # number of coin tosses
p = 0.75        # p parameter of the Bernoulli distribution
sbernoulli(n, p) = [(rand() < p) ? 1 : 0 for _ = 1:n] # define Bernoulli sampler
dataset = sbernoulli(N, p); # run N Bernoulli trials
print("dataset = ") ; show(dataset)
```

### [Model specification](@id getting-started-model-specification)
The next step is to specify our probabilistic model. In a Bayesian setting, this amounts to specifying the joint probability of the random variables of the system.

#### Likelihood

We will suppose that the outcome of each coin flip is modelled by the Bernoulli distribution, i.e.

``y_i \sim Bernoulli(\theta)``,

where ``y_i=1`` represents "heads", ``y_i=0`` represents "tails", and ``θ \in [0,1]`` is the underlying probability of the coin landing heads up for a single coin flip.

#### Prior
In order to simplify computation we will choose a prior distribution that is conjugate to the Bernoulli likelihood function defined above, which turns out to be the beta distribution, i.e.

``\theta \sim Beta(a, b)``,

where ``a`` and ``b`` are the hyper-parameters that encode our prior beliefs about the possible values of ``\theta``. We will assign values to the hyper-parameters in a later step.   

#### Joint probability
The joint probability is given by the multiplication of the likelihood and the prior, i.e.

``P(\{y_i\}, θ) = P(\{y_i\} | θ) P(θ)``.

Now let's see how to specify this model using ForneyLab's syntax.

```@example 1
using ForneyLab
g = FactorGraph()       # create a factor graph
a = placeholder(:a)     # define hyperparameter a as placeholder
b = placeholder(:b)     # define hyperparameter b as placeholder
@RV θ ~ Beta(a, b)      # prior
@RV y ~ Bernoulli(θ)    # likelihood
placeholder(y, :y)      # define y as a placeholder for data
draw(g)                 # draw the factor graph
```
As you can see, ForneyLab offers a model specification syntax that resembles closely to the mathematical equations defined above.

### Message-passing algorithm generation
Once we have defined our model, the next step is to instruct ForneyLab to generate a message-passing algorithm that solves our given inference problem. To do this, we need to specify which type of algorithm we want to use. In this case we will use *belief propagation*, also known as the *sum-product algorithm*. Once we execute the following code, we see that a function called `step!(...)` becomes available in the current scope. This function contains the sum-product message-passing algorithm.
```@example 1
# Generate a message passging sum-product algorithm that infers theta
algo_str = sumProductAlgorithm(θ) # ForneyLab returns the algorithm as a string
algorithm = Meta.parse(algo_str) # parse the algorithm into a Julia expression
eval(algorithm); # evaluate the functions contained in the Julia expression
print("algorithm = ") ; show(algorithm)
```

### Message-passing inference execution
The last step is to execute the message-passing algorithm obtained in the previous section. In order to do this, we first need to assign values to the hyper-parameters ``a`` and ``b`` which characterize our prior beliefs ``p(\theta)`` about the bias of the coin. Then, we need to feed the observations, one at a time, to the algorithm together with our current belief  about ``\theta`` acting as prior.

```@example 1
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
```@example 1
using Plots, LaTeXStrings, SpecialFunctions;
pyplot(fillalpha=0.3, leg=false, xlabel=L"\theta", yticks=nothing)
BetaPDF(α, β) = x ->  x^(α-1)*(1-x)^(β-1)/beta(α, β) # beta distribution definition
BernoulliPDF(z, N) = θ -> θ^z*(1-θ)^(N-z) # Bernoulli distribution definition

rθ = range(0, 1, length=100)
p1 = plot(rθ, BetaPDF(a, b), title="Prior", fill=true, ylabel=L"P(\theta)", c=1,)
p2 = plot(rθ, BernoulliPDF(sum(dataset), N), title="Likelihood", fill=true, ylabel=L"P(D|\theta)", c=2)
p3 = plot(rθ, BetaPDF(marginals[:θ].params[:a], marginals[:θ].params[:b]), title="Posterior", fill=true, ylabel=L"P(\theta|D)", c=3)
plot(p1, p2, p3, layout=@layout([a; b; c]))
```

## Where to go next?

There are a set of [demos](https://github.com/biaslab/ForneyLab.jl/tree/master/demo) available in ForneyLab's repository that demonstrate the more advanced features of ForneyLab. Alternatively, you can head to the [User guide](@ref) which provides more detailed information of how to use ForneyLab to solve inference problems.
