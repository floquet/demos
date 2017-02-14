
module load cuda
module swap compiler/intel/15.0.3 compiler/gcc/4.8.5
nvcc -o matmult.exe \
 -gencode arch=compute_35,code=sm_35 \
  matmult.cu

