
#include <stdlib.h>
#include <stdio.h>
#include <float.h>

#include <cuda.h>
#include <cuda_runtime.h>

#define BLOCK_SIZE 512

// Vector sum of C = A + B

__global__ void vecadd(double * A, double * B, double *C, int N) {
   int i = blockIdx.x*blockDim.x + threadIdx.x;
   if(i<N) { C[i] = A[i] + B[i]; }
}

int main(int argc, char** argv) {
   double * hostA;
   double * hostB;
   double * hostC;
   double * hostCompare;
   double * deviceA;
   double * deviceB;
   double * deviceC;
   int N = 10000;
   int i;
   double delta;
   int num_bad_values;

   // Allocate host matrices.
   hostA = (double *) malloc(N * sizeof(double));
   hostB = (double *) malloc(N * sizeof(double));
   hostC = (double *) malloc(N * sizeof(double));
   hostCompare = (double *) malloc(N * sizeof(double));

   // Fill-in values for hostA and hostB.
   for(i = 0; i < N; i++) {
     hostA[i] = 1.0 + i;
     hostB[i] = 1.0/(1.0 + i);
   }
   // Calculate on host for comparison.
   for(i = 0; i < N; i++) {
     hostCompare[i] = hostA[i] + hostB[i];
   }
   
   // Allocate GPU memory.
   cudaMalloc(&deviceA, N * sizeof(double));
   cudaMalloc(&deviceB, N * sizeof(double));
   cudaMalloc(&deviceC, N * sizeof(double));


   // Copy memory to the GPU.
   cudaMemcpy(deviceA, hostA,
              N * sizeof(double),
              cudaMemcpyHostToDevice);
   cudaMemcpy(deviceB, hostB,
              N * sizeof(double),
              cudaMemcpyHostToDevice);

    dim3 dimBlock(BLOCK_SIZE, 1, 1);
    dim3 dimGrid(((N - 1)/BLOCK_SIZE) + 1, 1, 1);

   // Launch the kernel.
   vecadd<<<dimGrid, dimBlock>>>(deviceA, deviceB, deviceC, N);

   cudaThreadSynchronize();

   // Copy the GPU memory back to the CPU .
   cudaMemcpy(hostC, deviceC,
              N * sizeof(double),
              cudaMemcpyDeviceToHost);

   delta = DBL_EPSILON;
   num_bad_values = 0;
   for(i = 0; i < N; i++) {
      if(fabs(hostC[i] - hostCompare[i]) > delta) {
	 num_bad_values++;
     }
   }

   fprintf(stdout," Number of bad values in C = %lu\n", num_bad_values);

   // Free the GPU memory.
   cudaFree(deviceA);
   cudaFree(deviceB);
   cudaFree(deviceC);
   // Free the CPU memory.
   free(hostA);
   free(hostB);
   free(hostC);
   free(hostCompare);

   return 0;
}

