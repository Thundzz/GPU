
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#ifdef __APPLE__
#  include <OpenCL/opencl.h>
#else
#  include <CL/opencl.h>
#endif
#include "util.h"

#define MAX_PLATFORMS 3
#define MAX_DEVICES   5

#define KERNEL_NAME  "heat"

#define ITERATIONS    100

unsigned SIZE = 4 * 1024 * 1024;
unsigned TILE = 16;

float *input_data, *output_data;

cl_mem input_buffer;  // device memory used for input data
cl_mem output_buffer;  // device memory used for output data

struct timeval tv1,tv2;

static void alloc_buffers_and_user_data(cl_context context)
{
  // CPU side
  input_data = calloc((SIZE+2), sizeof(float));
  output_data = calloc((SIZE+2), sizeof(float));

  for(int i = 0; i < SIZE; i++)
    input_data[i] = (float)i;

  // Allocate buffers inside device memory
  //
  input_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(float) * (SIZE+2), NULL, NULL);
  if (!input_buffer)
    error("Failed to allocate input buffer");

  output_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(float) * (SIZE+2), NULL, NULL);
  if (!output_buffer)
    error("Failed to allocate output buffer");
}



/* Version CPU pour comparer le resultat */
void heat(float* next, const float* previous) {
  for(unsigned int i=0;i<SIZE;i++)
    next[i+1] = (previous[i] + previous[i+1] * 2 + previous[i+2]) / 4;
}

static void check_output_data(void)
{
  /* Version cpu pour comparaison */
  float* reference = (float*) calloc((SIZE+2), sizeof(float));

  memset(reference, 0, (SIZE+2) * sizeof(float)) ;

  gettimeofday(&tv1,NULL);
 
 
  for(int i=0;i<ITERATIONS;++i) {
    heat(reference, input_data);
    heat(input_data, reference);
  }
  gettimeofday(&tv2,NULL);
  float time=((float)TIME_DIFF(tv1,tv2)) / 1000;
  
  printf("Seq :\t%f\n", time);
  
  /* Comparaison des resultats */
  unsigned int errors=0;
  for(unsigned int i=0;i<SIZE;i++){
    if((reference[i]-output_data[i])/reference[i] > 1e-6)
      errors++;
  }
  if(errors)
    fprintf(stderr,"%d erreurs !\n", errors);
  else
    fprintf(stderr,"pas d'erreurs, cool !\n"); 
}

static void free_buffers_and_user_data(void)
{
  free(input_data);
  free(output_data);
  
  clReleaseMemObject(input_buffer);
  clReleaseMemObject(output_buffer);
}

static void send_input(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueWriteBuffer(queue, input_buffer, CL_TRUE, 0,
			     sizeof(float) * (SIZE+2), input_data, 0, NULL, NULL);
  check(err, "Failed to write to input array");

  err = clEnqueueWriteBuffer(queue, output_buffer, CL_TRUE, 0,
			     sizeof(float) * (SIZE+2), output_data, 0, NULL, NULL);
  check(err, "Failed to write to output array");
}

static void retrieve_output(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0,
			    sizeof(float) * (SIZE+2), output_data, 0, NULL, NULL );  
  check(err, "Failed to read output array");
}

int main(int argc, char** argv)
{
  cl_platform_id pf[MAX_PLATFORMS];
  cl_uint nb_platforms = 0;
  cl_int err;                            // error code retuh_dataned from api calls
  cl_device_type device_type = CL_DEVICE_TYPE_ALL;

  // Filter args
  //
  argv++;
  while (argc > 1) {
    if(!strcmp(*argv, "-g") || !strcmp(*argv, "--gpu-only")) {
      if(device_type != CL_DEVICE_TYPE_ALL)
	error("--gpu-only and --cpu-only can not be specified at the same time\n");
      device_type = CL_DEVICE_TYPE_GPU;
    } else if(!strcmp(*argv, "-c") || !strcmp(*argv, "--cpu-only")) {
      if(device_type != CL_DEVICE_TYPE_ALL)
	error("--gpu-only and --cpu-only can not be specified at the same time\n");
      device_type = CL_DEVICE_TYPE_CPU;
    } else if(!strcmp(*argv, "-s") || !strcmp(*argv, "--size")) {
      unsigned i;
      int r;
      char c;

      r = sscanf(argv[1], "%u%[mMkK]", &SIZE, &c);

      if (r == 2) {
	if (c == 'k' || c == 'K')
	  SIZE *= 1024;
	else if (c == 'm' || c == 'M')
	  SIZE *= 1024 * 1024;
      }

      argc--; argv++;
    } else
      break;
    argc--; argv++;
  }

  if(argc > 1)
    TILE = atoi(*argv);

  // Get list of OpenCL platforms detected
  //
  err = clGetPlatformIDs(3, pf, &nb_platforms);
  check(err, "Failed to get platform IDs");

  printf("%d OpenCL platforms detected\n", nb_platforms);

  // For each platform do
  //
  for (cl_int p = 0; p < nb_platforms; p++) {
    cl_uint num;
    int platform_valid = 1;
    char name[1024], vendor[1024];
    cl_device_id devices[MAX_DEVICES];
    cl_uint nb_devices = 0;
    cl_context context;                 // compute context
    cl_program program;                 // compute program
    cl_kernel kernel;

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_NAME, 1024, name, NULL);
    check(err, "Failed to get Platform Info");

    err = clGetPlatformInfo(pf[p], CL_PLATFORM_VENDOR, 1024, vendor, NULL);
    check(err, "Failed to get Platform Info");

    printf("Platform %d: %s - %s\n", p, name, vendor);

    // Get list of devices
    //
    err = clGetDeviceIDs(pf[p], device_type, MAX_DEVICES, devices, &nb_devices);
    printf("nb devices = %d\n", nb_devices);

    if(nb_devices == 0)
      continue;

    // Create compute context with "device_type" devices
    //
    context = clCreateContext (0, nb_devices, devices, NULL, NULL, &err);
    check(err, "Failed to create compute context");

    // Load program source into memory
    //
    const char	*opencl_prog;
    opencl_prog = file_load(KERNEL_FILE);

    // Attach program source to context
    //
    program = clCreateProgramWithSource(context, 1, &opencl_prog, NULL, &err);
    check(err, "Failed to create program");

    // Compile program
    //
    {
      char flags[1024];

      sprintf (flags,
	       "-cl-mad-enable -cl-fast-relaxed-math -DSIZE=%d -DTILE=%d -DTYPE=%s",
	       SIZE, TILE, "float");

      err = clBuildProgram (program, 0, NULL, flags, NULL, NULL);
      if(err != CL_SUCCESS) {
	size_t len;

	// Display compiler log
	//
	clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, 0, NULL, &len);
	{
	  char buffer[len+1];

	  fprintf(stderr, "--- Compiler log ---\n");
	  clGetProgramBuildInfo(program, devices[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, NULL);
	  fprintf(stderr, "%s\n", buffer);
	  fprintf(stderr, "--------------------\n");
	}
	if(err != CL_SUCCESS)
	  error("Failed to build program!\n");
      }
    }

    // Create the compute kernel in the program we wish to run
    //
    kernel = clCreateKernel(program, KERNEL_NAME, &err);
    check(err, "Failed to create compute kernel");

    // Allocate and initialize input data
    //
    alloc_buffers_and_user_data(context);

    // Iterate over devices
    //
    for(cl_int dev = 0; dev < nb_devices; dev++) {
      cl_command_queue queue;

      char name[1024];
      cl_device_type dtype;

      err = clGetDeviceInfo(devices[dev], CL_DEVICE_NAME, 1024, name, NULL);
      check(err, "Cannot get type of device");
      err = clGetDeviceInfo(devices[dev], CL_DEVICE_TYPE, sizeof(cl_device_type), &dtype, NULL);
      check(err, "Cannot get type of device");

      printf("\tDevice %d : %s [%s]\n", dev, (dtype == CL_DEVICE_TYPE_GPU) ? "GPU" : "CPU", name);

      // Create a command queue
      //
      queue = clCreateCommandQueue(context, devices[dev], CL_QUEUE_PROFILING_ENABLE, &err);
      check(err,"Failed to create command queue");

      // Write our data set into device buffer
      //
      send_input(queue);

      // Execute kernel
      //
      {
	cl_event prof_event;
	cl_ulong start, end;
	struct timeval t1,t2;
	double timeInMicroseconds;
	size_t global[1] = { SIZE };                      // global domain size for our calculation
	size_t local[1]  = { TILE };                       // local domain size for our calculation

	printf("\t%zu Threads in workgroups of %zu\n", global[0], local[0]);


	gettimeofday (&t1, NULL);


  for(int i = 0; i<ITERATIONS; i++) // Iterations are done inside the kernel
      {
	// Set the arguments to our compute kernel
	//
	err = 0;
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buffer);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_buffer);
	check(err, "Failed to set kernel arguments!\n");

	err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global, local, 0, NULL, NULL);
	check(err, "Failed to execute kernel!\n");

	// Now the other way...
	//
	err = 0;
	err |= clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_buffer);
	err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buffer);
	check(err, "Failed to set kernel arguments!\n");

	err = clEnqueueNDRangeKernel(queue, kernel, 1, NULL, global, local, 0, NULL,&prof_event);
	check(err, "Failed to execute kernel!\n");
      }
	// Wait for the command commands to get serviced before reading back results
	//
	clFinish(queue);

	gettimeofday (&t2,NULL);

	// Check performance
	//
	timeInMicroseconds = (double)TIME_DIFF(t1, t2) / ITERATIONS;

	printf("\tComputation performed in %lf µs over device #%d\n",
	       timeInMicroseconds,
	       dev);

	clReleaseEvent(prof_event);
      }

      // Read back the results from the device to verify the output
      //
      retrieve_output(queue);

      // Validate computation
      //
      check_output_data();

      clReleaseCommandQueue(queue);
    }

    // Cleanup
    //
    free_buffers_and_user_data();

    clReleaseKernel(kernel);
    clReleaseProgram(program);
    clReleaseContext(context);
  }


  return 0;
}



