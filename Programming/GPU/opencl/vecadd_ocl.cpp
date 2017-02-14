
#include <CL/cl.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>

const char * get_error_string(cl_int err) {
         switch(err) {
             case 0: return "CL_SUCCESS";
             case -1: return "CL_DEVICE_NOT_FOUND";
             case -2: return "CL_DEVICE_NOT_AVAILABLE";
             case -3: return "CL_COMPILER_NOT_AVAILABLE";
             case -4: return "CL_MEM_OBJECT_ALLOCATION_FAILURE";
             case -5: return "CL_OUT_OF_RESOURCES";
             case -6: return "CL_OUT_OF_HOST_MEMORY";
             case -7: return "CL_PROFILING_INFO_NOT_AVAILABLE";
             case -8: return "CL_MEM_COPY_OVERLAP";
             case -9: return "CL_IMAGE_FORMAT_MISMATCH";
             case -10: return "CL_IMAGE_FORMAT_NOT_SUPPORTED";
             case -11: return "CL_BUILD_PROGRAM_FAILURE";
             case -12: return "CL_MAP_FAILURE";

             case -30: return "CL_INVALID_VALUE";
             case -31: return "CL_INVALID_DEVICE_TYPE";
             case -32: return "CL_INVALID_PLATFORM";
             case -33: return "CL_INVALID_DEVICE";
             case -34: return "CL_INVALID_CONTEXT";
             case -35: return "CL_INVALID_QUEUE_PROPERTIES";
             case -36: return "CL_INVALID_COMMAND_QUEUE";
             case -37: return "CL_INVALID_HOST_PTR";
             case -38: return "CL_INVALID_MEM_OBJECT";
             case -39: return "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
             case -40: return "CL_INVALID_IMAGE_SIZE";
             case -41: return "CL_INVALID_SAMPLER";
             case -42: return "CL_INVALID_BINARY";
             case -43: return "CL_INVALID_BUILD_OPTIONS";
             case -44: return "CL_INVALID_PROGRAM";
             case -45: return "CL_INVALID_PROGRAM_EXECUTABLE";
             case -46: return "CL_INVALID_KERNEL_NAME";
             case -47: return "CL_INVALID_KERNEL_DEFINITION";
             case -48: return "CL_INVALID_KERNEL";
             case -49: return "CL_INVALID_ARG_INDEX";
             case -50: return "CL_INVALID_ARG_VALUE";
             case -51: return "CL_INVALID_ARG_SIZE";
             case -52: return "CL_INVALID_KERNEL_ARGS";
             case -53: return "CL_INVALID_WORK_DIMENSION";
             case -54: return "CL_INVALID_WORK_GROUP_SIZE";
             case -55: return "CL_INVALID_WORK_ITEM_SIZE";
             case -56: return "CL_INVALID_GLOBAL_OFFSET";
             case -57: return "CL_INVALID_EVENT_WAIT_LIST";
             case -58: return "CL_INVALID_EVENT";
             case -59: return "CL_INVALID_OPERATION";
             case -60: return "CL_INVALID_GL_OBJECT";
             case -61: return "CL_INVALID_BUFFER_SIZE";
             case -62: return "CL_INVALID_MIP_LEVEL";
             case -63: return "CL_INVALID_GLOBAL_WORK_SIZE";
             default: return "Unknown OpenCL error";
         }
}

void get_build_error(cl_program prg, void * devId) {   
	char *ch;
        cl_device_id _devId;
	size_t length ;
	
	memcpy((cl_device_id *)&_devId,
               (const void *)devId, sizeof(cl_device_id));
    
	clGetProgramBuildInfo(prg, _devId, CL_PROGRAM_BUILD_LOG, 0, NULL, &length);
	ch = (char *)malloc(length * sizeof(char));
	clGetProgramBuildInfo(prg, _devId, CL_PROGRAM_BUILD_LOG, length, ch, &length);
	
	if (length > 2)
		printf("Got OpenCL kernel build error...\n%s", ch);
	
	free(ch);
}

const char * KernelSource =                          "\n" \
"__kernel void vecAdd(  __global float *a,            \n" \
"                       __global float *b,            \n" \
"                       __global float *c,            \n" \
"                       int n)                        \n" \
"{                                                    \n" \
"    int id = get_global_id(0);                       \n" \
"                                                     \n" \
"    if (id < n) {                                    \n" \
"        c[id] = a[id] + b[id];                       \n" \
"    }                                                \n" \
"}                                                    \n" \
"\n";

int main(int argc, char ** argv) {
   int inputLength;
   float * hostInput1;
   float * hostInput2;
   float * hostOutput;
   float * hostCompare;

   cl_int cl_err = CL_SUCCESS;
   cl_context clcntx;
   cl_device_id clDevice;
   cl_command_queue clcmdq;
   cl_mem cmDevIn1;
   cl_mem cmDevIn2;
   cl_mem cmDevOut;
   cl_program clpgm;
   cl_kernel clkrn;
   size_t parmsz;
   int i;
   float epsilon;

   inputLength = 100;
   hostInput1 = (float *) malloc(inputLength * sizeof(float));
   hostInput2 = (float *) malloc(inputLength * sizeof(float));
   hostOutput = (float *) malloc(inputLength * sizeof(float));
   hostCompare = (float *) malloc(inputLength * sizeof(float));

   for(i = 0; i < inputLength; i++) {
      hostInput1[i] = i + 1.0;
      hostInput2[i] = inputLength*2.0 - i;
      hostCompare[i] = hostInput1[i] + hostInput2[i];
   }

   cl_platform_id platform;
   cl_uint num_platforms;
   cl_err = clGetPlatformIDs(1, &platform, &num_platforms);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, from clGetPlatformIDs line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   cl_err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &clDevice, NULL);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, from clGetDeviceIDs line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   clcntx = clCreateContext(0, 1, &clDevice, NULL, NULL, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, from clCreateContext line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   clcmdq = clCreateCommandQueue(clcntx, clDevice, 0, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   cmDevIn1 = clCreateBuffer(clcntx, CL_MEM_READ_ONLY,
                             sizeof(cl_float)*inputLength,
                             NULL, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cmDevIn2 = clCreateBuffer(clcntx, CL_MEM_READ_ONLY,
                             sizeof(cl_float)*inputLength,
                             NULL, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cmDevOut = clCreateBuffer(clcntx, CL_MEM_WRITE_ONLY,
                             sizeof(cl_float)*inputLength,
                             NULL, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   clpgm = clCreateProgramWithSource(clcntx, 1,
                                     (const char **) &KernelSource,
                                     NULL, &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
 
   char cl_compiler_flags[4096];
   sprintf(cl_compiler_flags, "-cl-mad-enable");

   cl_err = clBuildProgram(clpgm, 0, NULL, cl_compiler_flags, NULL, NULL);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error clBuildProgram, line %u , file %s\n", __LINE__, __FILE__);
     get_build_error(clpgm, &clDevice);
     return 1;
   }

   clkrn = clCreateKernel(clpgm, "vecAdd", &cl_err);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     fprintf(stderr,"   %s\n", get_error_string(cl_err));
     fprintf(stderr,"KernelSource = %s", KernelSource);
     return 1;
   }

   cl_err = clSetKernelArg(clkrn, 0, sizeof(cl_mem), (void*)&cmDevIn1);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cl_err = clSetKernelArg(clkrn, 1, sizeof(cl_mem), (void*)&cmDevIn2);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cl_err = clSetKernelArg(clkrn, 2, sizeof(cl_mem), (void*)&cmDevOut);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cl_err = clSetKernelArg(clkrn, 3, sizeof(cl_int), (void*)&inputLength);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

   cl_err = clEnqueueWriteBuffer(clcmdq, cmDevIn1, CL_FALSE, 0,
                                 sizeof(cl_float)*inputLength,
                                 (const void *) hostInput1,
                                 0, NULL, NULL);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }
   cl_err = clEnqueueWriteBuffer(clcmdq, cmDevIn2, CL_FALSE, 0,
                                 sizeof(cl_float)*inputLength,
                                 (const void *) hostInput2,
                                 0, NULL, NULL);
   if(cl_err != CL_SUCCESS) {
     fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     return 1;
   }

    int threadsPerBlock = 256;
    int blocksPerGrid =
           (inputLength + threadsPerBlock - 1) / threadsPerBlock;

  cl_event event = NULL;
  size_t global_size = threadsPerBlock*blocksPerGrid;
  size_t local_size = threadsPerBlock;
  cl_err = clEnqueueNDRangeKernel(clcmdq, clkrn, 1, NULL,
                                   &global_size, &local_size,
                                   0, NULL, NULL);
  if(cl_err != CL_SUCCESS) {
    fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
    return 1;
  }

  cl_err = clFinish(clcmdq);
  if(cl_err != CL_SUCCESS) {
    fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
     fprintf(stderr,"   %s\n", get_error_string(cl_err));
    return 1;
  }

  cl_err = clEnqueueReadBuffer(clcmdq, cmDevOut, CL_TRUE, 0,
                                sizeof(cl_float)*inputLength,
                                (void *) hostOutput,
                                0, NULL, NULL);
  if(cl_err != CL_SUCCESS) {
    fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
    return 1;
  }

  epsilon = 0.0000001;
  float difference;
  int there_is_problem = 0;
   for(i = 0; i < inputLength; i++) {
     difference = fabs(hostCompare[i] - hostOutput[i]);
     if(difference > epsilon) { there_is_problem = 1; }
   }
   if(there_is_problem != 0) {
    fprintf(stderr,"Error, line %u , file %s\n", __LINE__, __FILE__);
    fprintf(stderr,"   comparison does not match\n");
   }
   else {
     fprintf(stdout,"Comparison Successful.\n");
   }

   free(hostInput1);
   free(hostInput2);
   free(hostOutput);

   clReleaseMemObject(cmDevIn1);
   clReleaseMemObject(cmDevIn2);
   clReleaseMemObject(cmDevOut);

   clReleaseKernel(clkrn);
   clReleaseProgram(clpgm);
   clReleaseCommandQueue(clcmdq);
   clReleaseContext(clcntx);

   return 0;
}


