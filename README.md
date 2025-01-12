# `mixedpower`
>Power for linear mixed effects models in Python

Calculating power in linear mixed effects designs is non-trivial. [G*Power](https://www.psychologie.hhu.de/arbeitsgruppen/allgemeine-psychologie-und-arbeitspsychologie/gpower) doesn't even do it. Yet many experiments have crossed random factors: for instance, participants and stimuli (targets). This "library" provides approximate power calculations for such designs, specifically focusing on a “CCC” layout:

```
C = Condition is within Participant
C = Condition is within Target
C = Condition is within both Participant and Target (i.e. repeated across both)
```

>[!NOTE]
>The code generated in this database is a Python conversion of [Jake Westfall's R code](https://github.com/jake-westfall/two_factor_power/tree/master). [Here](https://psycnet.apa.org/doiLanding?doi=10.1037%2Fxge0000014) is the accompanying publication. See [here](https://pmc.ncbi.nlm.nih.gov/articles/PMC6646942/) for additional resources.

# Installation 
To install 

```
pip install git+https://github.com/w-decker/mixedpower.git
```

# Usage

The core functionality of `mixedpower` is comprised of two methods: calculating power directly and solving for experimental variables based on desired power. 

## Calculating power directly

```python
import mixedpower as mp

args = dict(
    design='CCC', # totally crossed design (currently the only design supported)
    cohens_d=0.5, # effect size
    resid=0.3, # residual variabce
    target_intercept=0.2,
    participant_intercept=0.2,
    participant_x_target=0.1, # interaction
    target_slope=0.1,
    participant_slope=0.1,
    n_participants=40,
    n_targets=40,
    alpha=0.05
)

power, _ = mp.power(**args)
print(f'Empirical power: {power:.5f}')
```

## Solving for sample size

```python
import mixedpower as mp

args = dict(
    p=0.8, # desired power
    cohens_d = 0.5,
    resid=0.3,
    target_intercept=0.2,
    participant_intercept=0.2,
    participant_x_target = 0.1,
    target_slope=0.1,
    participant_slope=0.1,
    n_targets = 30,
    code=1, 
    alpha=0.05,
)

n_participants, _ = mp.solve(variable='n_participants', **args)
print(f'Number of participants: {n_participants}')
```

# Details

## Variance components

The variance components defined in the program are used to approximate the standard error of the condition effect.

## Degrees of freedom

This packages uses the [Welch-Satterthwaite](https://en.wikipedia.org/wiki/Welch%E2%80%93Satterthwaite_equation) approximation to calculate degrees of freedom.
