# Posterior summaries of STAGE transition points (m50)

Returns posterior summaries of \\m\_{50}\\, the value of \\x\\ where
\\P(y=1 \| x) = 0.5\\. Under equal class priors, \\m\_{50}\\ is the
direct sampling parameter of the STAGE model.

## Usage

``` r
# S3 method for class 'stage_fit'
transition_point(object, population = NULL, ...)
```

## Arguments

- object:

  A `stage_fit` object as returned by
  [`fit_stage()`](https://anhsmith.github.io/stage/reference/fit_stage.md).

- population:

  Optional population index (integer) or label (character). If `NULL`,
  return global and all populations.

- ...:

  Currently ignored. Included for method compatibility.

## Value

A named list. Always contains `global` (posterior summary of the global
`m50`). For J \> 1 models, also contains `pop1`, `pop2`, etc. for each
population-specific transition point. Each summary is a named numeric
vector with `mean`, `median`, `q2.5`, `q97.5`.
