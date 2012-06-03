#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <new>

#include "matrix_op.h"
#include "common_defs.h"
#include "inline_func.h"

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
  )
{
  double *pmatrix = data;

  for (int rs = node_index * block_size; rs < size; rs += node_count * block_size) {
    int re = (rs + block_size < size ? rs + block_size : size);
    for (int row = rs; row < re; ++row) {
      for (int col = 0; col < size; ++col) {
        switch (function) {
          case 1:
            *(pmatrix++) = size - max(row, col);
            break;
          case 2:
            if (row == col) {
              if (row == 0)
                *(pmatrix++) = 1;
              else
                *(pmatrix++) = 2;
            } else {
              if (abs(row - col) == 1) {
                *(pmatrix++) = -1;
              } else {
                *(pmatrix++) = 0;
              }
            }
            break;
          case 3:
            *(pmatrix++) = 1. / (row + col + 1);
            break;
          case -1:
            if (row > col) {
              *(pmatrix++) = 0;
            } else {
              pmatrix++;
            }
            break;
          case -2:
            if (row < col) {
              *(pmatrix++) = 0;
            } else {
              pmatrix++;
            }
            break;
          case 0:
          default:
            if (row == col)
              *(pmatrix++) = 1;
            else
              *(pmatrix++) = 0;
        }
      }
    }
  }
}

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
  )
{
  int status = 0;
  int pos = 0;
  int proc_id = 0;
  if (node_index == 0) {
    FILE *input = NULL;
    status = ((input = fopen(filename, "r")) == NULL);
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status) {
      if (input)
        fclose(input);
      return -1;
    }
    double *buf = new double [block_size * size];
    for (int rs = 0; rs < size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      for (int row = rs; row < re; ++row) {
        for (int col = 0; col < size; ++col) {
          status = (fscanf(input, "%lf", &buf[(row - rs) * size + col]) != 1);
          if (status) {
            row = re;
            col = size;
          }
        }
      }
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status) {
        delete[] buf;
        fclose(input);
        return -2;
      }
      if (proc_id != node_index) {
        MPI_Send(buf, (re - rs) * size, MPI_DOUBLE, proc_id,
                 0, MPI_COMM_WORLD);
      } else {
        memcpy(data + pos * size, buf, (re - rs) * size * sizeof(double));
        pos += re - rs;
      }
      proc_id = (proc_id + 1) % node_count;
    }
    delete[] buf;
    fclose(input);
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status)
      return -1;
    MPI_Status recv_status;
    for (int rs = 0; rs < size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status)
        return -2;
      if (proc_id == node_index) {
        MPI_Recv(data + pos * size, (re - rs) * size, MPI_DOUBLE, 0,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &recv_status);
        pos += re - rs;
      }
      proc_id = (proc_id + 1) % node_count;
    }
  }
  return 0;
}

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
  )
{
  int status = 0;
  int pos = 0;
  int proc_id = 0;
  if (node_index == 0) {
    FILE *output = NULL;
    status = ((output = fopen (filename, "w")) == NULL);
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status) {
      if (output)
        fclose(output);
      return -1;
    }
    double *buf = new double [block_size * size];
    MPI_Status recv_status;
    for (int rs = 0; rs < size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      if (proc_id != node_index) {
        MPI_Recv(buf, (re - rs) * size, MPI_DOUBLE, proc_id,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &recv_status);
      } else {
        memcpy(buf, data + pos * size, (re - rs) * size * sizeof(double));
        pos += re - rs;
      }
      for (int row = rs; row < re; ++row) {
        for (int col = 0; col < size; ++col) {
          status = (fprintf(output, "%lf ", buf[(row - rs) * size + col]) < 0);
          if (status)
            col = size;
        }
        if (status == 0)
          status = (fprintf(output, "\n") < 0);
        if (status)
          row = re;
      }
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status) {
        delete[] buf;
        fclose(output);
        return -2;
      }
      proc_id = (proc_id + 1) % node_count;
    }
    delete[] buf;
    fclose(output);
  } else {
    MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (status)
      return -1;
    for (int rs = 0; rs < size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      if (proc_id == node_index) {
        MPI_Send(data + pos * size, (re - rs) * size, MPI_DOUBLE, 0,
                 0, MPI_COMM_WORLD);
        pos += re - rs;
      }
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status)
        return -2;
      proc_id = (proc_id + 1) % node_count;
    }
  }
  return 0;
}

/*
 * Print part of matrix A on display
 */
void print_matrix(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  )
{
  int status = 0;
  int pos = 0;
  int proc_id = 0;
  int print_size = min(PRINT_LEN, size);

  if (node_index == 0) {
    double *buf = new double [block_size * size];
    MPI_Status recv_status;
    for (int rs = 0; rs < print_size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      if (proc_id != node_index) {
        MPI_Recv(buf, (re - rs) * size, MPI_DOUBLE, proc_id,
                 MPI_ANY_TAG, MPI_COMM_WORLD, &recv_status);
      } else {
        memcpy(buf, data + pos * size, (re - rs) * size * sizeof(double));
        pos += re - rs;
      }
      for (int row = rs; row < min(re, print_size); ++row) {
        for (int col = 0; col < print_size; ++col) {
          status = (printf("%7.3lf ", buf[(row - rs) * size + col]) < 0);
          if (status)
            col = size;
        }
        if (status == 0)
          status = (printf("\n") < 0);
        if (status)
          row = re;
      }
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status) {
        delete[] buf;
        return;
      }
      proc_id = (proc_id + 1) % node_count;
    }
    delete[] buf;
  } else {
    for (int rs = 0; rs < print_size; rs += block_size) {
      int re = (rs + block_size < size ? rs + block_size : size);
      if (proc_id == node_index) {
        MPI_Send(data + pos * size, (re - rs) * size, MPI_DOUBLE, 0,
                 0, MPI_COMM_WORLD);
        pos += re - rs;
      }
      MPI_Bcast(&status, 1, MPI_INT, 0, MPI_COMM_WORLD);
      if (status)
        return;
      proc_id = (proc_id + 1) % node_count;
    }
  }
}

// Convert down triangular matrix to symmetric
void downtr_to_symm(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  )
{
  int proc_id = 0;
  double *buf = new double [block_size * size];
  int pos_memcpy = 0;
  for (int rs = 0; rs < size; rs += block_size) {
    int re = (rs + block_size < size ? rs + block_size : size);
    if (node_index == proc_id) {
      for (int row = rs; row < re; ++row) {
        memcpy(buf + (row - rs) * re, data + pos_memcpy * size, re * sizeof(double));
        ++pos_memcpy;
      }
    }
    MPI_Bcast(buf, (re - rs) * re, MPI_DOUBLE, proc_id, MPI_COMM_WORLD);
    double *pdata = data;
    double *pbuf = buf;
    for (int cs  = node_index * block_size; cs < re;
             cs += node_count * block_size) {
      int ce = (cs + block_size < re ? cs + block_size : re);
      for (int col = cs; col < ce; ++col) {
        pbuf = buf + max(col + 1 - rs, 0) * re;
        for (int row = max(col + 1, rs); row < re; ++row, pbuf += re) {
          pdata[row] = pbuf[col];
        }
        pdata += size;
      }
    }
    proc_id = (proc_id + 1) % node_count;
  }
  delete[] buf;
}

// Calculates inf matrix norm
double inf_norm_matrix(
  int size,
  int block_size,
  int node_index,
  int node_count,
  double *data
  )
{
  int pos = 0;
  double buf = 0;
  for (int rs = node_index * block_size; rs < size; rs += node_count * block_size) {
    int re = (rs + block_size < size ? rs + block_size : size);
    for (int row = rs; row < re; ++row) {
      buf = max(buf, sum_abs_array(size, data + pos * size));
      ++pos;
    }
  }
  double result = 0;
  MPI_Allreduce(&buf, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  return result;
}

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
  )
{
  int memlen = size + node_count * block_size - 1;
  memlen /= node_count * block_size;
  memlen *= block_size * size;
  memset(res, 0, memlen * sizeof(double));

  if (second <= first && first < second + memlen) {
    return -1;
  }
  int proc_id = 0;
  MPI_Status status;

  int left  = (node_index - 1 + node_count) % node_count;
  int right = (node_index + 1) % node_count;
  for (int index = 0; index < node_count; ++index) {
    proc_id = (node_index + index) % node_count;
    int pos_first = 0;
    for (int is  = node_index * block_size; is < size;
         is += node_count * block_size) {
      int ie = (is + block_size < size ? is + block_size : size);
      int pos_second = 0;
      for (int js  = proc_id * block_size;
               js < size;
               js += node_count * block_size) {
        int je = (js + block_size < size ? js + block_size : size);
        for (int ks = 0; ks < size; ks += block_size) {
          int ke = (ks + block_size < size ? ks + block_size : size);
          add_tpart_mul(ie - is, je - js, ke - ks, size, size, size,
                        first + pos_first * size + ks,
                        second + pos_second * size + ks,
                        res + pos_first * size + js);
        }
        pos_second += je - js;
      }
      pos_first += ie - is;
    }
    if (node_count > 1) {
      MPI_Sendrecv_replace(second, memlen, MPI_DOUBLE, left, 0, right, MPI_ANY_TAG,
                           MPI_COMM_WORLD, &status);
    }
  }
  return 0;
}

