
#include <stdlib.h>
#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>

// Compute C = A*B

__global__ void matrixMultiply(float * A, float * B, float *C,
                    int numAColumns, int numBColumns,
                    int numCRows, int numCColumns) {
   int Row = blockIdx.y * blockDim.y + threadIdx.y;
   int Col = blockIdx.x * blockDim.x + threadIdx.x;
   if( (Row < numCRows) && (Col < numCColumns) ) {
      float Cvalue = 0.0;
      for(int i = 0; i < numAColumns; ++i) {
         Cvalue += A[Row*numAColumns + i]*B[i*numBColumns + Col];
      }
      C[Row*numCColumns + Col] = Cvalue;
   } // if( (Row < numCRows) && (Col < numCColumns) )
}

int main(int argc, char** argv) {
   float * hostA;  // The A matrix
   float * hostB;  // The B matrix
   float * hostC;  // The output C matrix
   float * hostCompare;
   float * deviceA;
   float * deviceB;
   float * deviceC;
   int numARows; // number of rows in the matrix A
   int numAColumns; // number of columns in the matrix A
   int numBRows; // number of rows in the matrix B
   int numBColumns; // number of columns in the matrix B
   int numCRows; // number of rows in the matrix C
   int numCColumns; // number of columns in the matrix C

   int i, j, k, offsetA, offsetB, offsetC;
   size_t num_bad_values;
   float sum, delta;

   numARows = 100;
   numAColumns = 50;
   numBRows = 50;
   numBColumns = 200;
   
   // Set numCRows and numCColumns
   numCRows = numARows;
   numCColumns = numBColumns;

   // Allocate host matrices.
   hostA = (float *) malloc(numARows * numAColumns * sizeof(float));
   hostB = (float *) malloc(numBRows * numBColumns * sizeof(float));
   hostC = (float *) malloc(numCRows * numCColumns * sizeof(float));
   hostCompare = (float *) malloc(numCRows * numCColumns * sizeof(float));

   // Fill-in values for hostA and hostB.
   for(i = 0; i < numARows; i++) {
     for(j = 0; j < numAColumns; j++) {
       offsetA = j + i*numAColumns;
       hostA[offsetA] = 2.0*i + 3.0*j;
     }
   }
   for(i = 0; i < numBRows; i++) {
     for(j = 0; j < numBColumns; j++) {
       offsetB = j + i*numBColumns;
       hostB[offsetB] = 5.0*i + 7.0*j;
     }
   }
   // Calculate matrix multiplication on CPU.
   for(i = 0; i < numCColumns; i++) {
     for(j = 0; j < numCRows; j++) {
       offsetC = i + j*numCColumns;
       sum = 0.0;
       for(k = 0; k < numAColumns; k++) {
	 offsetA = k + j*numAColumns;
	 offsetB = i + k*numBColumns;
	 sum += hostA[offsetA]*hostB[offsetB];
       }
       hostCompare[offsetC] = sum;
     }
   }
   
   // Allocate GPU memory.
   cudaMalloc(&deviceA, numARows * numAColumns * sizeof(float));
   cudaMalloc(&deviceB, numBRows * numBColumns * sizeof(float));
   cudaMalloc(&deviceC, numCRows * numCColumns * sizeof(float));


   // Copy host memory to the GPU.
   cudaMemcpy(deviceA, hostA,
              numARows * numAColumns * sizeof(float),
              cudaMemcpyHostToDevice);
   cudaMemcpy(deviceB, hostB,
              numBRows * numBColumns * sizeof(float),
              cudaMemcpyHostToDevice);

   // Initialize the grid and block dimensions.
   // 16x16 is a typical value, it could be changed.
   dim3 threadsPerBlock(16, 16);
   dim3 numBlocks((numCColumns + threadsPerBlock.x - 1)/threadsPerBlock.x,
                  (numCRows + threadsPerBlock.y - 1)/threadsPerBlock.y);

   // Launch the kernel.
   matrixMultiply<<<numBlocks, threadsPerBlock>>>(deviceA, deviceB, deviceC,
                    numAColumns, numBColumns,
                    numCRows, numCColumns);

   cudaThreadSynchronize();

   // Copy the GPU memory back to the CPU .
   cudaMemcpy(hostC, deviceC,
              numCRows * numCColumns * sizeof(float),
              cudaMemcpyDeviceToHost);

   delta = 0.0000001;
   num_bad_values = 0;
   for(i = 0; i < numCColumns; i++) {
     for(j = 0; j < numCRows; j++) {
       offsetC = i + j*numCColumns;
       if(fabs(hostC[offsetC] - hostCompare[offsetC]) > delta) {
	 num_bad_values++;
       }
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

