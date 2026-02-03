#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix compare_geno_cpp(NumericMatrix data) {
  int num_samples = data.nrow();
  int num_markers = data.ncol();

  // Initialize the result matrix
  NumericMatrix matching_matrix(num_samples, num_samples);

  // Loop through all pairs of samples
  for (int i = 0; i < num_samples; ++i) {
    for (int j = i; j < num_samples; ++j) {
      double matching_markers = 0.0;
      double valid_markers_count = 0.0;

      // Compare markers for the current pair
      for (int k = 0; k < num_markers; ++k) {
        double val1 = data(i, k);
        double val2 = data(j, k);

        // Check for NA values (R's NA is a specific double value)
        if (!R_IsNA(val1) && !R_IsNA(val2)) {
          valid_markers_count++;
          if (val1 == val2) {
            matching_markers++;
          }
        }
      }

      double proportion = 0.0;
      if (valid_markers_count > 0) {
        proportion = matching_markers / valid_markers_count;
      } else {
        proportion = NA_REAL; // Set to NA if no valid markers
      }

      // Store the result
      matching_matrix(i, j) = proportion;
      matching_matrix(j, i) = proportion;
    }
  }

  return matching_matrix;
}
