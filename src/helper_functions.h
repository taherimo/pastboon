#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include "uthash.h"
#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdbool.h>

/*
 * This block of code was taken from the BoolNet package.
 * Original file: BoolNet/src/common.h
*/

// The start of the code taken from BoolNet

#define BITS_PER_BLOCK_32 (sizeof(unsigned int) * 8)

#define GET_BIT(x, i) (((x) & ((unsigned long long)1 << (i))) != 0)

#define SET_BIT(x, i) ((x) | ((unsigned long long)1 << (i)))

#define CLEAR_BIT(x, i) ((x) & (~((unsigned long long)1 << (i))))

#define SET_BIT_TO_VAL(x, i, v)                                                \
  (((x) & (~((unsigned long long)1 << (i)))) | ((v) << (i)))

#define GET_BIT_ARRAY(x, i)                                                    \
  (((*(&(x) + i / BITS_PER_BLOCK_32)) &                                        \
    ((unsigned int)1 << (i % BITS_PER_BLOCK_32))) != 0)

#define SET_BIT_ARRAY(x, i)                                                    \
  (*(&(x) + i / BITS_PER_BLOCK_32) |=                                          \
   ((unsigned int)1 << (i % BITS_PER_BLOCK_32)))

#define CLEAR_BIT_ARRAY(x, i)                                                  \
  (*(&(x) + i / BITS_PER_BLOCK_32) &=                                          \
   (~((unsigned int)1 << (i % BITS_PER_BLOCK_32))))

typedef struct {

  void *ptr;
  UT_hash_handle hh;

} AllocatedMemory;

extern AllocatedMemory *memoryMap;

static inline void *CALLOC(size_t n, size_t sz) {
  void *ptr = calloc(n, sz);

  if (ptr == NULL)
    error("Out of memory!");

  AllocatedMemory *m = calloc(1, sizeof(AllocatedMemory));
  m->ptr = ptr;
  HASH_ADD_PTR(memoryMap, ptr, m);
  return ptr;
}

static inline void FREE(void *ptr) {
  AllocatedMemory *m;
  HASH_FIND_PTR(memoryMap, &ptr, m);
  HASH_DEL(memoryMap, m);
  free(m);
  free(ptr);
}

extern void bin2decC(int *dec, int *bin, int *numBits);

extern void dec2binC(int *bin, int *dec, int *numBits);

// Returns a random double in [0,1)
static inline double doublerand_1(void) { return unif_rand(); }

// Returns a random integer value in [0,maxVal-1]
static inline unsigned int intrand(unsigned int maxVal) {
  return (unsigned int)(unif_rand() * maxVal);
}

// The end of the code taken from BoolNet

extern int areArraysEqual(unsigned int arr1[], unsigned int arr2[],
                          unsigned int size);

// Returns a random integer value in [0,UINT_MAX]
static inline unsigned int intrand_fullrange(void) {
  return (unsigned int)(unif_rand() * (UINT_MAX + 1));
}

#endif
