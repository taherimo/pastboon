#include "helper_functions.h"
#include <stdio.h>
#include <string.h>

/*
 * This block of code was taken from the BoolNet package.
 * Original file: BoolNet/src/common.c
 */

// The start of the code taken from BoolNet

AllocatedMemory *memoryMap = NULL;

void freeAllMemory(void) {
  AllocatedMemory *m, *tmp;
  HASH_ITER(hh, memoryMap, m, tmp) {
    HASH_DEL(memoryMap, m);
    free(m->ptr);
    free(m);
  }
  // Rprintf("Freed all memory\n");
}

void bin2decC(int *dec, int *bin, int *numBits) {
  // clear output first
  unsigned int numElts;
  if (*numBits % BITS_PER_BLOCK_32 == 0)
    numElts = *numBits / BITS_PER_BLOCK_32;
  else
    numElts = *numBits / BITS_PER_BLOCK_32 + 1;

  memset(dec, 0, numElts * sizeof(int));

  // decode input and write binary integers
  unsigned int *unsigned_dec = (unsigned int *)dec;
  unsigned int i;

  for (i = 0; i < *numBits; ++i) {
    unsigned_dec[i / BITS_PER_BLOCK_32] |=
        ((unsigned int)bin[i] << (i % BITS_PER_BLOCK_32));
  }
}

void dec2binC(int *bin, int *dec, int *numBits) {
  unsigned int i;
  unsigned int *unsigned_dec = (unsigned int *)dec;

  for (i = 0; i < *numBits; ++i)
    if ((unsigned_dec[i / BITS_PER_BLOCK_32] &
         ((unsigned int)1 << (i % BITS_PER_BLOCK_32))) != 0)
      bin[i] = 1;
    else
      bin[i] = 0;
}

// The end of the code taken from BoolNet

int areArraysEqual(unsigned int arr1[], unsigned int arr2[],
                   unsigned int size) {
  for (int i = 0; i < size; i++) {
    if (arr1[i] != arr2[i]) {
      return 0;
    }
  }
  return 1;
}
