
module load cuda
module swap compiler/intel/15.0.3 compiler/gcc/4.8.5
nvcc -c vecadd_ocl.cpp
g++ -o vecadd_ocl.exe \
  -L/p/home/apps/cuda/7.5/lib64 -lOpenCL \
   vecadd_ocl.o

