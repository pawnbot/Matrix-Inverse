#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "matrix_op.h"
#include "common_defs.h"
#include "inline_func.h"

/*
 * Decompose part A, using Cholesky Decomposition A = res * res^t
 * Answer ll be stored in matrix array
 * Return negative value if error
 */
static int decompose_part_cholesky (
  int n,                    // order of system
  int size,                 // line size of original matrix
  double *matrix            // (n x n) symmetric part
  ) 
{ 
  // min on diagonal
  const double min_diagonal = FLOATING_POINT_MINIMAL_FOR_DIVISION;
  
  double *pa = matrix;
  for (int i = 0; i < n; ++i, pa += size) {
    for (int k = 0; k < i; ++k)
      pa[i] -= pa[k] * pa[k];   
   
    if (pa[i] < min_diagonal)
      return -1;
    
    pa[i] = sqrt(pa[i]);

    double mult = 1 / pa[i];
    for (int j = i + 1; j < n; ++j) {
      // r_{j i} = r_{i i}^{-1} * (a_{j i} - sum_{k=1}^{i-1}{r_{j k} * r_{i k})
      matrix[j * size + i] -= scalar_product(i, matrix + j * size, matrix + i * size);
      matrix[j * size + i] *= mult;
    }
  } 
  return 0;
}

/*
 * Inverse downtr matrix part A
 * Answer ll be stored in matrix array
 * Return negative value if error
 */
static int inverse_downtr_part (
  int n,                    // order of system
  int size,                 // line size of matrix
  double *matrix,           // (n x n) downtr part
  double *result            // (n x n) downtr buf
  ) 
{ 
  // min on diagonal
  const double min_diagonal = FLOATING_POINT_MINIMAL_FOR_DIVISION;
  
  //result = I
  double *pa = result;
  for (int row = 0; row < n; ++row, pa += size) {
    for (int col = 0; col < n; ++col) {
      if (row == col)
        pa[col] = 1;
      else
        pa[col] = 0;
    }
  }
  
  pa = matrix;
  for (int k = 0; k < n; ++k, pa += size) {
    for (int j = 0; j < k; ++j) {
      sub_scaled_array(j + 1, result + j * size, pa[j], result + k * size);
    }
    if (fabs(pa[k]) < min_diagonal)
      return -1;
    scale_array(k + 1, result + k * size, 1 / pa[k]);
    // copy result to matrix
    memcpy(pa, result + k * size, n * sizeof(double));
  }

  return 0;
}

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
  )
{
  // generate_matrix(2, size, block_size, node_index, node_count, data);
  int proc_id = 0;
  int pos = 0;
  double *buf = new double [block_size * size];
  double *block = new double [block_size * block_size];
  
  for (int cs = 0; cs < size; cs += block_size) {
    int ce = (cs + block_size < size ? cs + block_size : size);
    if (proc_id == node_index) {
      for (int ks = 0; ks < cs; ks += block_size) {
        // R_{i i} R_{i i}^t = A_{i i} - sum_{k} R{i k} * R{i k}^t
        int ke = (ks + block_size < size ? ks + block_size : size);
        sub_tpart_mul(ce - cs, ce - cs, ke - ks, size, size, size,
                      data + pos * size + ks,
                      data + pos * size + ks,
                      data + pos * size + cs);
      }
      // Decompose R_{i i} R{i i}^t
      if (decompose_part_cholesky(ce - cs, size, data + pos * size + cs)) {
        data[pos * size + cs] = -100500;
      } else {
        // Inverse down triangular matrix R_{i i}
        if (inverse_downtr_part(ce - cs, size, data + pos * size + cs, buf)) {
          data[pos * size + cs] = -100500;
        }
      }
      // copy line to buf
      for (int row = 0; row < ce - cs; ++row) {
        memcpy(buf + row * ce, data + pos * size, ce * sizeof(double));
        ++pos;
      }
    }

    MPI_Bcast(buf, (ce - cs) * ce, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);

    if (buf[cs] < -1) {
      delete[] buf;
      delete[] block;
      return -1;
    }
    int line_pos = 0;
    for (int rs  = node_index * block_size; rs < size;
             rs += node_count * block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      if (rs < ce) {
        line_pos += re - rs;
        continue;
      }
      for (int ks = 0; ks < cs; ks += block_size) {
        // R_{j i} R_{i i}^t = A_{j i} - sum_{k} R{j k} * R{i k}^t
        int ke = (ks + block_size < size ? ks + block_size : size);
        sub_tpart_mul(re - rs, ce - cs, ke - ks, size, ce, size,
                      data + line_pos * size + ks,
                      buf + ks,
                      data + line_pos * size + cs);
      }
      // R{j i} = (R_{j i} R_{i i}^t) * R{i i}^-t
      memset(block, 0, block_size * block_size * sizeof(double));
      add_tpart_mul(re - rs, ce - cs, ce - cs, size, ce, block_size,
                    data + line_pos * size + cs,
                    buf + cs,
                    block);
      for (int row = 0; row < re - rs; ++row) {
        memcpy(data + line_pos * size + cs, block + row * block_size, 
               (ce - cs) * sizeof(double));
        ++line_pos;
      }
    }
    proc_id = (proc_id + 1) % node_count; 
  }
  // decomposition done!
  int memlen = size + node_count * block_size - 1;
  memlen /= node_count * block_size;
  memlen *= block_size * size;
  double *workspace = new double [memlen];
  generate_matrix(0, size, block_size, node_index, node_count, workspace);
  proc_id = 0;
  pos = 0;
  for (int rs  = 0; rs < size; rs += block_size) {
    int re = (rs + block_size < size ? rs + block_size : size);
    if (node_index == proc_id) {
      for (int ks = 0; ks < re; ks += block_size) {
        int ke = (ks + block_size < re ? ks + block_size : re);
        // block = E_{r k}
        copy_to_block(re - rs, ke - ks, size, workspace + pos * size + ks, block);
        // E_{r k} = R_{r r}^{-1} E_{r k}
        part_mul(re - rs, re - rs, ke - ks, size, ke - ks, size,
                 data + pos * size + rs,
                 block,
                 workspace + pos * size + ks);
      }
      for (int row = 0; row < re - rs; ++row) {
        memcpy(buf + row * re, workspace + pos * size, re * sizeof(double));
        ++pos;
      }
    }

    MPI_Bcast(buf, (re - rs) * re, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);

    int line_pos = 0;
    for (int ks  = node_index * block_size; ks < size; 
             ks += node_count * block_size, line_pos += block_size) {
      int ke = (ks + block_size < size ? ks + block_size : size);
      if (ks < re)
        continue;

      for (int cs = 0; cs < re; cs += block_size) {
        int ce = (cs + block_size < re ? cs + block_size : re);
        // E_{k c} -= R{k r} * E{r c}
        sub_part_mul(ke - ks, re - rs, ce - cs, size, re, size,
                     data + line_pos * size + rs,
                     buf + cs,
                     workspace + line_pos * size + cs);
      }
    }
    proc_id = (proc_id + 1) % node_count;
  }
  // backtrace done!
  memset(data, 0, memlen * sizeof(double));
  downtr_to_symm(size, block_size, node_index, node_count, workspace);
  generate_matrix(-1, size, block_size, node_index, node_count, workspace);
  proc_id = 0;
  pos = 0;
  int line_size = 0;
  for (int rs  = 0; rs < size; rs += block_size) {
    int re = (rs + block_size < size ? rs + block_size : size);
    line_size = size - rs;
    if (proc_id == node_index) {
      for (int row = 0; row < re - rs; ++row) {
        memcpy(buf + row * line_size, workspace + pos * size + rs, 
               line_size * sizeof(double));
        ++pos;
      }
    }

    MPI_Bcast(buf, (re - rs) * line_size, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);

    int line_pos = 0;
    for (int ks  = node_index * block_size; ks < size; 
             ks += node_count * block_size, line_pos += block_size) {
      int ke = (ks + block_size < size ? ks + block_size : size);
      if (ks < rs)
        continue;
      for (int cs = ks; cs < size; cs += block_size) {
        int ce = (cs + block_size < size ? cs + block_size : size);
        // data_{k r} += workspace_{k c} * buf{r c}
        add_tpart_mul(ke - ks, re - rs, ce - cs, size, line_size , size,
                      workspace + line_pos * size + cs,
                      buf + cs - rs,
                      data + line_pos * size + rs);
      }
    }
    proc_id = (proc_id + 1) % node_count;
  }
  
  delete[] buf;
  delete[] block;
  delete[] workspace;
  return 0;
}
