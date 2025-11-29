# Predict class probabilities from a STAGE model

Predict class probabilities from a STAGE model

## Usage

``` r
# S3 method for class 'stage_fit'
predict(object, newdata, group = NULL, type = c("prob", "class"), ...)
```

## Arguments

- object:

  A `stage_fit` object as returned by
  [`fit_stage()`](https://anhsmith.github.io/stage/reference/fit_stage.md).

- newdata:

  Numeric vector of new x values.

- group:

  Optional group for each new x (length 1 or `length(newdata)`). If
  `NULL` and there are multiple populations in the model, population 1
  is used for all predictions.

- type:

  `"prob"` for P(y = 1 \| x), `"class"` for 0/1 classification.

- ...:

  Currently ignored.

## Value

Numeric vector of probabilities or classes.
