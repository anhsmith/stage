# Posterior summaries of STAGE transition points (m50)

Posterior summaries of STAGE transition points (m50)

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

## Value

A named list or numeric vector with posterior summary stats.
