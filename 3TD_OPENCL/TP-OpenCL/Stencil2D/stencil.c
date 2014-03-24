
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

#define ITERATIONS    100

// Application-specific data
//
#define KERNEL_NAME  "stencil"

unsigned SIZE = 1024;
unsigned TILE = 16;

#define C0           0.75
#define C1           0.25
#define BORDER       1
#define PADDING      (64/sizeof(float) - 2 * BORDER)

int LINESIZE;
int OFFSET;
int TOTALSIZE;

float *input_data, *output_data, *tmp;

cl_mem input_buffer;  // device memory used for input data
cl_mem output_buffer;  // device memory used for output data


static void alloc_buffers_and_user_data(cl_context context)
{
  // CPU side
  input_data  = malloc(TOTALSIZE * sizeof(float));
  output_data = malloc(TOTALSIZE * sizeof(float));
  tmp         = malloc(TOTALSIZE * sizeof(float));

  srand(1234);
  for(int y = 0; y < SIZE; y++)
    for(int x = 0; x < SIZE; x++) {
      input_data[OFFSET + LINESIZE * y + x] = rand();
    }

  // Allocate buffers inside device memory
  //
  input_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(float) * TOTALSIZE, NULL, NULL);
  if (!input_buffer)
    error("Failed to allocate input buffer");

  output_buffer = clCreateBuffer(context,  CL_MEM_READ_WRITE,  sizeof(float) * TOTALSIZE, NULL, NULL);
  if (!output_buffer)
    error("Failed to allocate output buffer");
}


static void stencil(float *A, float* B)
{
  A += OFFSET;
  B += OFFSET;

  for(int y = 0; y < SIZE; y++)
    for(int x = 0; x < SIZE; x++) {
      B[ y * LINESIZE + x] = 
	C0 * A[ y * LINESIZE + x] + 
	C1 * ( A[ y * LINESIZE + x - 1] +
	       A[ y * LINESIZE + x + 1] +
	       A[ (y - 1) * LINESIZE + x] +
	       A[ (y + 1) * LINESIZE + x]);
    }
}

static void check_output_data(void)
{
  for(int i = 0; i < ITERATIONS/2; i++) {
    stencil(input_data, tmp);
    stencil(tmp, input_data);
  }

  unsigned correct = 0;

  for(int y = 0; y < SIZE; y++)
    for(int x = 0; x < SIZE; x++) {
      float expected = tmp[OFFSET + LINESIZE * y + x];
      float observed = output_data[OFFSET + LINESIZE * y + x];

      if(abs((observed - expected)/expected) < 1e-6)
	correct++;
    }
  printf("\tComputed '%d/%d' correct values!\n", correct, SIZE * SIZE);
}

static void free_buffers_and_user_data(void)
{
  free(input_data);
  free(output_data);
  free(tmp);

  clReleaseMemObject(input_buffer);
  clReleaseMemObject(output_buffer);
}

static void send_input(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueWriteBuffer(queue, input_buffer, CL_TRUE, 0,
			     sizeof(float) * TOTALSIZE, input_data, 0, NULL, NULL);
  check(err, "Failed to write to input array");
}

static void retrieve_output(cl_command_queue queue)
{
  cl_int err;

  err = clEnqueueReadBuffer(queue, output_buffer, CL_TRUE, 0,
			    sizeof(float) * TOTALSIZE, output_data, 0, NULL, NULL );  
  check(err, "Failed to read output array");
}

int main(int argc, char** argv)
{
  cl_platform_id pf[MAX_PLATFORMS];
  cl_uint nb_platforms = 0;
  cl_int err;                            // error code returned from api calls
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

  LINESIZE =  SIZE + PADDING + 2 * BORDER;
  OFFSET =    LINESIZE + 64/sizeof(float);
  TOTALSIZE = LINESIZE * (SIZE + 2 * BORDER);

  printf("OFFSET=%d\nLINESIZE=%d\nTOTALSIZE=%d\n", OFFSET, LINESIZE, TOTALSIZE);

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
	       "-cl-mad-enable -cl-fast-relaxed-math -DSIZE=%d -DTILE=%d -DTYPE=%s "
	       "-DOFFSET=%d -DLINESIZE=%d -DC0=%f -DC1=%f",
	       SIZE, TILE, "float", OFFSET, LINESIZE, C0, C1);

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
	cl_ulong start, end;
	struct timeval t1,t2;
	double timeInMicroseconds;
	size_t global[2] = { SIZE, SIZE };  // global domain size for our calculation
	size_t local[2]  = { TILE, TILE };   // local domain size for our calculation

	printf("\t%dx%d Threads in workgroups of %dx%d\n", global[0], global[1], local[0], local[1]);


	gettimeofday (&t1, NULL);

	for (unsigned iter = 0; iter < ITERATIONS/2; iter++) {

	  // Set kernel arguments
	  //
	  err = 0;
	  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &input_buffer);
	  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &output_buffer);
	  check(err, "Failed to set kernel arguments");
	  err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global, local,
				       0, NULL, NULL);
	  check(err, "Failed to execute kernel");

	  // Now the other way
	  //
	  err = 0;
	  err  = clSetKernelArg(kernel, 0, sizeof(cl_mem), &output_buffer);
	  err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &input_buffer);
	  check(err, "Failed to set kernel arguments");
	  err = clEnqueueNDRangeKernel(queue, kernel, 2, NULL, global, local,
				       0, NULL, NULL);
	  check(err, "Failed to execute kernel");

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
