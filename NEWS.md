# svycdiff 0.2.0

* **Added**: support for binary and count outcomes via generalized estimating equations

* **Added**: Augmented inverse probability weighted (AIPW) estimator to main function via `id_form = "DR"`

* **Improved**: Expanded simulation function to allow for binary and count outcomes, stratified clustered sampling, and increased parameters to specify the true propensity, selection, and outcome models.

* **Updated**: Main function to be faster for larger sample sizes. Separated pieces of estimation workflow into helper functions (helpers.R) for better compartmentalization.

# svycdiff 0.1.0

* Initial CRAN submission.
