# Internals

This page summarizes some of ForneyLab's internal structures. It is mainly tailored for software developers interested in a high-level overview of the inner workings of the package. Coding style conventions can be found in `STYLEGUIDE.md`.

## Directory structure

ForneyLab's directories and files are structured as follows:

- `/demo/`: demos in Jupyter (iJulia) notebook format (`.ipynb`)
- `/docs/`: documentation source and build location
- `/src/`: ForneyLab source files
    - `algorithms/`: inference algorithm implementations
    - `engines/`: rendering of message passing schedules to executable code
        - `julia/`: Julia engine and update rule implementations
    - `factor_nodes/`: all node-specific files
    - `update_rules/`: message passing update rules
- `/test/`: test files with directory structure similar to `/src/`.


## InferenceAlgorithm data structure

A ForneyLab `InferenceAlgorithm` is structured as follows:

- `algorithm::InferenceAlgorithm`: specifies everything required for algorithm generation
    - `posterior_factorization::PosteriorFactorization`: specifies factorization of the posterior
        - `posterior_factors::Vector{PosteriorFactor}` (per item):
            - `id::Symbol`: posterior factor id
            - `optimize::Bool`: require optimization block
            - `initialize::Bool`: require initialization block
            - `schedule::Schedule` (per `ScheduleEntry` item):
                - `schedule_index::Int`: position of entry in schedule
                - `message_update_rule::Type`: update rule type for message computation
                - `initialize::Bool`: require message initialization
                - `family::FactorFunction`: family of message distribution (for initialization)
                - `dimensionality::Tuple`: dimensionality of message distribution (for initialization)
                - `inbounds::Vector` (per item):
                    - `inbound::Union`: inbound, see below
            - `marginal_table::MarginalTable` (per `MarginalEntry` item):
                - `marginal_id::Symbol`: identifier for the marginal
                - `marginal_update_rule::Union{Nothing, Product, Type}`: update rule type for marginal computation
                - `inbounds::Vector` (per item):
                    - `inbound::Union`: inbound, see below
    - `average_energies::Dict` (per item):
        - `node::Type`: node type for average energy computation
        - `counting_number::Int64`: counting number for average energy
        - `inbounds::Vector` (per item):
            - `inbound::Union`: inbound, see below
    - `entropies::Dict` (per item):
        - `counting_number::Int64`: counting number for (joint) entropy
        - `inbound::Union`: inbound, see below

Inbounds are of type `Union{Nothing, ScheduleEntry, MarginalEntry, Dict, Clamp}`.


## Update rules naming convention

The name of an update rule is composed of several parts:
1. The word `rule`
2. Type of algorithm
    - `SP`: sum-product
    - `VB`: variational Bayes
    - `SVB`: structured variational Bayes
    - `M`: marginal (joint)
3. Type of factor node
4. Interface of the outgoing message
5. Types of incoming messages (absent for `VB` rules)
    - `N`: Nothing
    - `P`: point mass
    - `D`: distribution
    - `[I]`: first letter of the message's probability distribution


###### Example 1: `ruleSPGaussianMeanPrecisionMPNP`
1. `rule` : update rule
2. `SP` : sum-product algorithm
3. `GaussianMeanPrecision`: Gaussian mean precision factor node
4. `M`: outgoing message through 'Mean' interface
5. `PNP`: incoming message types are: point mass, Nothing and point mass



###### Example 2: `ruleVBBernoulliOut`
1. `rule`: update rule
2. `VB`: variational Bayes algorithm
3. `Bernoulli`: Bernoulli factor node
4. `Out`: outgoing message through 'Out' interface
5. `-`


###### Example 3: `ruleEPProbitIn1GB`
1. `rule`: update rule
2. `EP`: expectation propagation algorithm
3. `Probit`: probit factor node
4. `In1`: outgoing message through 'in1' interface
5. `GB`: incoming message types are: Gaussian and Bernoulli
* *Note that EP update rules do not have `N` (nothing) in the set of incoming messages given that in EP there is an incoming message through the interface of the outgoing message that is being calculated.*
