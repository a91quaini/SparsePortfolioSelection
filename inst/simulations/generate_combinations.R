# Function to compute all combinations of integers from 0 to n-1 in k places
# and write them to a CSV file
compute_combinations <- function(n, k, file_name = "combinations.csv") {
  # Check that k is not larger than n
  if (k > n) {
    stop("k must be less than or equal to n")
  }

  # Create a vector of integers from 0 to n-1
  numbers <- 0:(n - 1)

  # Compute the combinations using combn().
  # The combn() function returns a matrix with 'k' rows (each column is one combination).
  # Transpose it so that each row becomes one combination with k columns.
  comb_matrix <- t(combn(numbers, k))

  # Construct the file path (assuming the 'src' directory exists in your package)
  file_path <- file.path("src", file_name)

  # Write the matrix to a CSV file. 'row.names = FALSE' avoids adding row numbers.
  write.csv(comb_matrix, file = file_path, row.names = FALSE)

  # Optionally, print a message indicating that the file was saved successfully.
  message("Combinations saved to: ", file_path)

  # Return the matrix (in case you need to use it within R)
  return(comb_matrix)
}

# Create all combinations of n choose k with n=20 and k=5,10,15
result <- compute_combinations(20, 5, "combinations_20choose5")
result <- compute_combinations(20, 10, "combinations_20choose10")
result <- compute_combinations(20, 15, "combinations_20choose15")
