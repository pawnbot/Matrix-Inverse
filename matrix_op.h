#ifndef MATRIX_OP_H
#define MATRIX_OP_H

/*
 * Generate (n x n) matrix A
 */
void generate_matrix(
  int function,
  int size,
  int block_size,
  int node_index,
  int node_count,               
  double *data
  );
  
/*
 * Read matrix A from file filename
 * Return negative value if error
 */
int read_matrix_file(
  const char* filename,
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );

/*
 * Print matrix A to file filename
 * Return negative value if error
 */
int print_matrix_file(
  const char* filename,
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );

/*
 * Print part of matrix A on display
 */
void print_matrix(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );

// Convert down triangular matrix to symmetric
void downtr_to_symm(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );
    
// Calculates inf matrix norm
double inf_norm_matrix(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );
  
// Calculates matrix res = A * B^t
// if second overlaps with first returns -1
int tmatrix_mul(
  int size,
  int block_size,
  int node_index,
  int node_count,
  const double *first,
  double *second,
  double *res
  );
  
/*
 * Calculates A^{-1} for SPD matrix A, data will be filled with A^{-1}
 * Return negative value if error
 */
int spd_inverse(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  );

#endif // MATRIX_OP_H
