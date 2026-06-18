###############################################################
## Run circular benchmark analyses
###############################################################

source("R/circular_density_models.R")
source("R/circular_benchmark_models.R")

## Example usage:
## theta <- your_circular_sample_in_radians
## res <- compare_circular_functional_models(theta, dataset_name = "dataset")
## write.csv(res$table, "outputs/circular_dataset_results.csv", row.names = FALSE)

message("Circular benchmark functions loaded. Provide datasets and call compare_circular_functional_models().")
