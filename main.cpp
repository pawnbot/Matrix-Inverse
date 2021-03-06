#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "common_defs.h"
#include "timer.h"
#include "inline_func.h"
#include "matrix_op.h"

#define MATRIX_INDEX 1
// #define PRINT_RESULT_TO_FILE
const char* output = "result.txt";

char buf[1024];

int main(int argc, char *argv[]) {
  int node_index, node_count;
  int size;
  int block_size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &node_index);
  MPI_Comm_size(MPI_COMM_WORLD, &node_count);

  if (node_index == 0) {
    sprintf(buf, "result%d.txt", node_count);
    freopen(buf, "w", stdout);
  }

  // Timer started
  MPI_Barrier(MPI_COMM_WORLD);
  if (node_index == 0)
    timer_start();
  

  if (argc != 3 && argc != 4) {
    if (node_index == 0)
      printf("Usage: %s matrix_size block_size [input_file]\n", argv[0]);
    MPI_Finalize();
    return -1;
  } else {
    size = atoi(argv[1]);
    block_size = atoi(argv[2]);
    if ((size < 2) || (block_size < 1) || (block_size > size)) {
      if (node_index == 0) {
        printf("Wrong format matrix_size > 1, matrix_size >= block_size > 0\n");
        printf("Usage: %s matrix_size block_size [input_file]\n", argv[0]);
      }
      MPI_Finalize();
      return -2;
    }
  }
  int memlen = size + node_count * block_size - 1;
  memlen /= node_count * block_size;
  memlen *= block_size * size;
  
  double *data = new double [memlen];
  memset(data, 0, memlen * sizeof(double));
  if (argc == 3) {
    generate_matrix(MATRIX_INDEX, size, block_size, node_index, node_count, data);
  } else {
    if (read_matrix_file(argv[3], size, block_size, node_index, node_count, data)) {
      if (node_index == 0)
        printf("Cannot open(read from) file: %s\n", argv[3]);
      delete[] data;
      MPI_Finalize();
      return -3;
    }
  }
  print_matrix(size, block_size, node_index, node_count, data);
  // Initialization done!
  MPI_Barrier(MPI_COMM_WORLD);
  if (node_index == 0)
    print_full_time("on init");

  if (spd_inverse(size, block_size, node_index, node_count, data)) {
    if (node_index == 0)
      printf("Method cannot be applied!\n");
    delete[] data;
    MPI_Finalize();
    return -4;
  }
  
  // Algorithm done!
  MPI_Barrier(MPI_COMM_WORLD);
  if (node_index == 0)
    print_full_time("on algorithm");
  
  // Printing result on console (and in file)
  downtr_to_symm(size, block_size, node_index, node_count, data);
  print_matrix(size, block_size, node_index, node_count, data);
  
  double *workspace = new double[memlen];
  memset(workspace, 0, memlen * sizeof(double));
  // For this two matrices inverse is known
  if (MATRIX_INDEX == 1 || MATRIX_INDEX == 2) {
    if (MATRIX_INDEX == 1)
      generate_matrix(2, size, block_size, node_index, node_count, workspace);
    else
      generate_matrix(1, size, block_size, node_index, node_count, workspace);
    sub_array(memlen, data, workspace);
    double err_norm = inf_norm_matrix(size, block_size, 
                                      node_index, node_count, workspace);
    err_norm = 0;
    // if (node_index == 0)
    //   printf("Error = %11.5le\n", err_norm);
  }
  // Restore input matrix to calculate residual  
  if (argc == 3) {
    generate_matrix(MATRIX_INDEX, size, block_size, node_index, node_count, workspace);
  } else {
    if (read_matrix_file(argv[3], size, block_size, node_index, node_count, workspace)) {
      if (node_index == 0)
        printf("Cannot open(read from) file: %s\n", argv[3]);
      delete[] data;
      delete[] workspace;
      MPI_Finalize();
      return -5;
    }
  }

#ifdef PRINT_RESULT_TO_FILE  
  if (print_matrix_file(output, size, block_size, node_index, node_count, data)) {
    if (node_index == 0)
      printf("Cannot open(print to) file: %s\n", output);
    delete[] data;
    delete[] workspace;
    MPI_Finalize();
    return -6;
  }
#endif

  double *residual = new double[memlen];
  tmatrix_mul(size, block_size, node_index, node_count, workspace, data, residual);
  generate_matrix(-2, size, block_size, node_index, node_count, residual);
  generate_matrix(0, size, block_size, node_index, node_count, workspace);
  sub_array(memlen, workspace, residual);
  downtr_to_symm(size, block_size, node_index, node_count, residual);

  double res_norm = inf_norm_matrix(size, block_size, node_index, node_count, residual);
  if (node_index == 0)
    printf("Residual = %11.5le\n", res_norm);

  delete[] data;
  delete[] workspace;
  delete[] residual;
  MPI_Finalize();
  return 0;
}

