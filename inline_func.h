#ifndef INLINE_FUNC_H
#define INLINE_FUNC_H

// swap two elements
template<typename T>
static inline void 
swap(T& a, T& b)
{
    T t = a;
    a = b;
    b = t;
}

// minimum of two elements
template<typename T>
static inline T 
min(const T a, const T b) 
{
    return (a < b ? a : b);
}

// maximum of two elements
template<typename T>
static inline T 
max(const T a, const T b) 
{
    return (a > b ? a : b);
}

// swap elements of arrays: src <-> dst
template<typename T>
static inline void
swap_arrays(
  int n,         // array size
  T *src,        // first array
  T *dst         // second array
  )
{
  int i, m;
  T temp;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++) {
    temp = dst[i]; dst[i] = src[i]; src[i] = temp;
  }

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    temp = dst[i]; dst[i] = src[i]; src[i] = temp;
    temp = dst[i+1]; dst[i+1] = src[i+1]; src[i+1] = temp;
    temp = dst[i+2]; dst[i+2] = src[i+2]; src[i+2] = temp;
    temp = dst[i+3]; dst[i+3] = src[i+3]; src[i+3] = temp;
    temp = dst[i+4]; dst[i+4] = src[i+4]; src[i+4] = temp;
    temp = dst[i+5]; dst[i+5] = src[i+5]; src[i+5] = temp;
    temp = dst[i+6]; dst[i+6] = src[i+6]; src[i+6] = temp;
    temp = dst[i+7]; dst[i+7] = src[i+7]; src[i+7] = temp;
  }
}

// Subtract arrays: dst = src - dst
template<typename T>
static inline void
sub_arrays(
  int n,         // array size
  const T *src,  // input array
  T *dst         // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    dst[i] = src[i] - dst[i];

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    dst[i]     = src[i] - dst[i];
    dst[i + 1] = src[i + 1] - dst[i + 1];
    dst[i + 2] = src[i + 2] - dst[i + 2];
    dst[i + 3] = src[i + 3] - dst[i + 3];
    dst[i + 4] = src[i + 4] - dst[i + 4];
    dst[i + 5] = src[i + 5] - dst[i + 5];
    dst[i + 6] = src[i + 6] - dst[i + 6];
    dst[i + 7] = src[i + 7] - dst[i + 7];
  }
}

// Subtract arrays: dst -= src
template<typename T>
static inline void
sub_array( 
  int n,              // array size
  const T *src,       // input array
  T *dst              // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; ++i)
    dst[i] -= src[i];

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    dst[i]     -= src[i];
    dst[i + 1] -= src[i + 1];
    dst[i + 2] -= src[i + 2];
    dst[i + 3] -= src[i + 3];
    dst[i + 4] -= src[i + 4];
    dst[i + 5] -= src[i + 5];
    dst[i + 6] -= src[i + 6];
    dst[i + 7] -= src[i + 7];
  }
}

// Add array src to dst: dst += src
template<typename T>
static inline void
add_array( 
  int n,              // array size
  const T *src,       // input array
  T *dst              // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; ++i)
    dst[i] += src[i];

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    dst[i]     += src[i];
    dst[i + 1] += src[i + 1];
    dst[i + 2] += src[i + 2];
    dst[i + 3] += src[i + 3];
    dst[i + 4] += src[i + 4];
    dst[i + 5] += src[i + 5];
    dst[i + 6] += src[i + 6];
    dst[i + 7] += src[i + 7];
  }
}

// Subtract arrays: dst = dst - src * mult
template<typename T>
static inline void
sub_scaled_array( 
  int n,              // array size
  const T *src,       // input array
  T mult,             // multiplier
  T *dst              // output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    dst[i] -= mult * src[i];

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    dst[i]     -= mult * src[i];
    dst[i + 1] -= mult * src[i + 1];
    dst[i + 2] -= mult * src[i + 2];
    dst[i + 3] -= mult * src[i + 3];
    dst[i + 4] -= mult * src[i + 4];
    dst[i + 5] -= mult * src[i + 5];
    dst[i + 6] -= mult * src[i + 6];
    dst[i + 7] -= mult * src[i + 7];
  }
}

// Scale array a by mult
template<typename T>
static inline void
scale_array( 
  int n,         // array size
  T *a,          // input/output array
  T mult         // multiplicator
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  for (i = 0; i < m; i++)
    a[i] *= mult;

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    a[i]     *= mult;
    a[i + 1] *= mult;
    a[i + 2] *= mult;
    a[i + 3] *= mult;
    a[i + 4] *= mult;
    a[i + 5] *= mult;
    a[i + 6] *= mult;
    a[i + 7] *= mult;
  }
}

// Scale array res by array mult
template<typename T>
static inline void
scale_array_by_array (
  int n,           // array size
  const T *mult,   // mult array
  T *res           // input/output array
  )
{
  int i, m;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  for (i = 0; i < m; i++)
    res[i] *= mult[i];

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    res[i]     *= mult[i];
    res[i + 1] *= mult[i + 1];
    res[i + 2] *= mult[i + 2];
    res[i + 3] *= mult[i + 3];
    res[i + 4] *= mult[i + 4];
    res[i + 5] *= mult[i + 5];
    res[i + 6] *= mult[i + 6];
    res[i + 7] *= mult[i + 7];
  }
}

// Compute euclidean scalar product of arrays
template<typename T>
static inline T
scalar_product (
  int n,            // array size
  const T *a,       // (n) array
  const T *b        // (n) array
  )
{
  int i, m;
  T res = 0;

  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  // Rest of loop
  for (i = 0; i < m; i++)
    res += a[i] * b[i];

  // Unrolled loop
  for (i = m; i < n; i += 8)
    res += a[i + 0] * b[i + 0] + a[i + 1] * b[i + 1]
         + a[i + 2] * b[i + 2] + a[i + 3] * b[i + 3]
         + a[i + 4] * b[i + 4] + a[i + 5] * b[i + 5]
         + a[i + 6] * b[i + 6] + a[i + 7] * b[i + 7];

  return res;
}

// Sum of absolute values in array
template<typename T>
static inline T
sum_abs_array(
  int n,               // array size
  const T *a           // input array
  )
{
  int i, m;
  T res = 0;
  
  // Clean-up loop so remaining vector length is a multiple of 8.
  m = n & 7;
  for (i = 0; i < m; i++)
    res += fabs (a[i]);

  // Unrolled loop
  for (i = m; i < n; i += 8) {
    res += fabs (a[i]);
    res += fabs (a[i + 1]);
    res += fabs (a[i + 2]);
    res += fabs (a[i + 3]);
    res += fabs (a[i + 4]);
    res += fabs (a[i + 5]);
    res += fabs (a[i + 6]);
    res += fabs (a[i + 7]);
  }
  
  return res;
}

// Copy part of matrix to block (n x m)
template<typename T>
static inline void
copy_to_block(
  int n,
  int m, 
  int line_size, 
  const T *src, 
  T *dst  // (n x m) block
  )
{
  for (int row = 0; row < n; ++row) {
    memcpy(dst + row * m, src + row * line_size, m * sizeof(T));
  }
}

// Copy block (n x m) to part of matrix
template<typename T>
static inline void
copy_from_block(
  int n,
  int m, 
  int line_size, 
  const T *src,  // (n x m) block
  T *dst
  )
{
  for (int row = 0; row < n; ++row) {
    memcpy(dst + row * line_size, src + row * m, m * sizeof(T));
  }
}

// res += A * B^t
template<typename T>
static inline void
add_tpart_mul(
  int n,
  int m,
  int k,
  int size_A,
  int size_B,
  int size_res,
  const T *A,  // (n x k) matrix part
  const T *B,  // (m x k) matrix part
  T *res       // (n x m) matrix part
  )
{
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < m; ++col)
      res[row * size_res + col] += scalar_product(k, A + row * size_A, B + col * size_B);
  }
}

// res -= A * B^t
template<typename T>
static inline void
sub_tpart_mul(
  int n,
  int m,
  int k,
  int size_A,
  int size_B,
  int size_res,
  const T *A,  // (n x k) matrix part
  const T *B,  // (m x k) matrix part
  T *res       // (n x m) matrix part
  )
{
  for (int row = 0; row < n; ++row) {
    for (int col = 0; col < m; ++col)
      res[row * size_res + col] -= scalar_product(k, A + row * size_A, B + col * size_B);
  }
}

// res = A * B
template<typename T>
static inline void
part_mul(
  int n,
  int m,
  int k,
  int size_A,
  int size_B,
  int size_res,
  const T *A,  // (n x m) matrix part
  const T *B,  // (m x k) matrix part
  T *res       // (n x k) matrix part
  )
{
  T *pres = res;
  const T *pa = A;
  for (int row = 0; row < n; ++row) {
    memset(pres, 0, k * sizeof(T));
    for (int i = 0; i < m; ++i) {
      sub_scaled_array(k, B + i * size_B, -pa[i], pres);
    }
    pres += size_res;
    pa += size_A;
  }
}

// res -= A * B
template<typename T>
static inline void
sub_part_mul(
  int n,
  int m,
  int k,
  int size_A,
  int size_B,
  int size_res,
  const T *A,  // (n x m) matrix part
  const T *B,  // (m x k) matrix part
  T *res       // (n x k) matrix part
  )
{
  T *pres = res;
  const T *pa = A;
  for (int row = 0; row < n; ++row) {
    for (int i = 0; i < m; ++i) {
      sub_scaled_array(k, B + i * size_B, pa[i], pres);
    }
    pres += size_res;
    pa += size_A;
  }
}

#endif // INLINE_FUNC_Ha
