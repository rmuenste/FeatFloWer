#include <cuda_runtime_api.h>

#include <stddef.h>
#include <stdio.h>

void *ff_hypre_malloc_managed(size_t bytes)
{
  void *pointer = NULL;
  cudaError_t error;

  if (bytes == 0)
  {
    return NULL;
  }

  error = cudaMallocManaged(&pointer, bytes, cudaMemAttachGlobal);
  if (error != cudaSuccess)
  {
    fprintf(stderr, "cudaMallocManaged failed: %s\n", cudaGetErrorString(error));
    return NULL;
  }
  return pointer;
}

int ff_hypre_free_managed(void *pointer)
{
  cudaError_t error;

  if (pointer == NULL)
  {
    return 0;
  }

  error = cudaFree(pointer);
  if (error != cudaSuccess)
  {
    fprintf(stderr, "cudaFree failed: %s\n", cudaGetErrorString(error));
    return (int) error;
  }
  return 0;
}
