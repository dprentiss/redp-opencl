/**********************************************************************
Copyright ©2013 Advanced Micro Devices, Inc. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

•	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
•	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
********************************************************************/

#include <sstream>
#include <float.h>
#include "redp.hpp"
#include <ctime>

/*
* Allocate and initialize memory 
* on the host. Print input array. 
*/
int initializeHost(void)
{
  // width               = (int) pow(2, NUM_OPTIONS);
  width               = (int) pow(2, 6);
  // width               = 128;
  output              = NULL;

  /////////////////////////////////////////////////////////////////
  // Allocate and initialize memory used by host 
  /////////////////////////////////////////////////////////////////
  std::cout << "\nAllocate and initialize memory used by host... \n"; 

  cl_uint optionFloatBytes = width * sizeof(cl_float);
  /*
  cl_uint nodeIntBytes = NUM_NODES * sizeof(cl_uint);
  cl_uint nodeFloatBytes = NUM_NODES * sizeof(cl_float);
  cl_uint arcIntBytes = NUM_ARCS * sizeof(cl_uint);
  cl_uint arcFloatBytes = NUM_ARCS * sizeof(cl_float);
  cl_uint arcBoolBytes = NUM_ARCS * sizeof(cl_bool);
  */

  output = (cl_float *) malloc(optionFloatBytes);
  if(!output)
  {
    std::cout << "Error: Failed to allocate input memory on host\n";
    return REDP_FAILURE;
  }

  for(cl_uint i = 0; i < width; i++)
    output[i] = 0;
  
  return REDP_SUCCESS;
}

/*
* Converts the contents of a file into a string
*/
std::string convertToString(const char *filename)
{
  size_t size;
  char*  str;
  std::string s;

  std::fstream f(filename, (std::fstream::in | std::fstream::binary));

  if(f.is_open())
  {
    size_t fileSize;
    f.seekg(0, std::fstream::end);
    size = fileSize = (size_t)f.tellg();
    f.seekg(0, std::fstream::beg);

    str = new char[size+1];
    if(!str)
    {
      f.close();
      std::cout << "Memory allocation failed";
      return NULL;
    }

    f.read(str, fileSize);
    f.close();
    str[size] = '\0';

    s = str;
    delete[] str;
    return s;
  }
  else
  {
    std::cout << "\nFile containg the kernel code(\".cl\") not found. Please copy the required file in the folder containg the executable.\n";
    exit(1);
  }
  return NULL;
}

/*
* \brief OpenCL related initialization 
*        Create Context, Device list, Command Queue
*        Create OpenCL memory buffer objects
*        Load CL file, compile, link CL source 
*		  Build program and kernel objects
*/
int initializeCL(void)
{
  cl_int status = 0;
  size_t deviceListSize;

  //////////////////////////////////////////////////////////////////// 
  // STEP 1 Getting Platform.
  //////////////////////////////////////////////////////////////////// 

  /*
  * Have a look at the available platforms and pick either
  * the AMD one if available or a reasonable default.
  */

  std::cout << "Getting platforms... \n";

  cl_uint numPlatforms;
  cl_platform_id platform = NULL;
  status = clGetPlatformIDs(0, NULL, &numPlatforms);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: Getting Platforms. (clGetPlatformsIDs)\n";
    return REDP_FAILURE;
  }

  if(numPlatforms > 0)
  {
    cl_platform_id* platforms = new cl_platform_id[numPlatforms];
    status = clGetPlatformIDs(numPlatforms, platforms, NULL);
    if(status != CL_SUCCESS)
    {
      std::cout << "Error: Getting Platform Ids. (clGetPlatformsIDs)\n";
      return REDP_FAILURE;
    }
    for(unsigned int i=0; i < numPlatforms; ++i)
    {
      char pbuff[100];
      status = clGetPlatformInfo(
        platforms[i],
        CL_PLATFORM_VENDOR,
        sizeof(pbuff),
        pbuff,
        NULL);
      if(status != CL_SUCCESS)
      {
        std::cout << "Error: Getting Platform Info.(clGetPlatformInfo)\n";
        return REDP_FAILURE;
      }
      platform = platforms[i];
      if(!strcmp(pbuff, "Advanced Micro Devices, Inc."))
      {
        break;
      }
    }
    delete platforms;
  }

  if(NULL == platform)
  {
    std::cout << "NULL platform found so Exiting Application." << std::endl;
    return REDP_FAILURE;
  }


  //////////////////////////////////////////////////////////////////// 
  // STEP 2 Creating context using the platform selected
  //      Context created from type includes all available
  //      devices of the specified type from the selected platform 
  //////////////////////////////////////////////////////////////////// 


  /*
  * If we could find our platform, use it. Otherwise use just available platform.
  */

  std::cout << "Creating context... \n";
  cl_context_properties cps[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };

  context = clCreateContextFromType(cps, 
    CL_DEVICE_TYPE_GPU, 
    NULL, 
    NULL, 
    &status);
  if(status != CL_SUCCESS) 
  {  
    std::cout << "Error: Creating Context. (clCreateContextFromType)\n";
    return REDP_FAILURE; 
  }


  //////////////////////////////////////////////////////////////////// 
  // STEP 3
  //      3.1 Query context for the device list size,
  //      3.2 Allocate that much memory using malloc or new
  //      3.3 Again query context info to get the array of device
  //              available in the created context
  //////////////////////////////////////////////////////////////////// 

  // First, get the size of device list data
  status = clGetContextInfo(context, 
    CL_CONTEXT_DEVICES, 
    0, 
    NULL, 
    &deviceListSize);
  if(status != CL_SUCCESS) 
  {  
    std::cout <<
      "Error: Getting Context Info \
      (device list size, clGetContextInfo)\n";
    return REDP_FAILURE;
  }

  devices = (cl_device_id *)malloc(deviceListSize);
  if(devices == 0)
  {
    std::cout << "Error: No devices found.\n";
    return REDP_FAILURE;
  }

  // Now, get the device list data
  status = clGetContextInfo(
    context, 
    CL_CONTEXT_DEVICES, 
    deviceListSize, 
    devices, 
    NULL);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Getting Context Info \
      (device list, clGetContextInfo)\n";
    return REDP_FAILURE;
  }

  //////////////////////////////////////////////////////////////////// 
  // STEP 4 Creating command queue for a single device
  //      Each device in the context can have a 
  //      dedicated commandqueue object for itself
  //////////////////////////////////////////////////////////////////// 

  std::cout << "Creating command queue... \n";
  commandQueue = clCreateCommandQueue(
    context, 
    devices[0], 
    0, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Creating Command Queue. (clCreateCommandQueue)\n";
    return REDP_FAILURE;
  }

  /////////////////////////////////////////////////////////////////
  // STEP 5 Creating cl_buffer objects from host buffer
  //          These buffer objects can be passed to the kernel
  //          as kernel arguments
  /////////////////////////////////////////////////////////////////
  std::cout << "Creating buffer objects... \n";

  pointerBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_NODES,
    pointer, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (pointerBuffer)\n";
    return REDP_FAILURE;
  }

  demandBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_NODES,
    demand, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (demandBuffer)\n";
    return REDP_FAILURE;
  }

  riskABuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * NUM_NODES,
    riskA, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (riskABuffer)\n";
    return REDP_FAILURE;
  }

  riskBBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * NUM_NODES,
    riskB, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (riskBBuffer)\n";
    return REDP_FAILURE;
  }

  tailBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_ARCS,
    tail, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (tailBuffer)\n";
    return REDP_FAILURE;
  }

  headBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_ARCS,
    head, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (headBuffer)\n";
    return REDP_FAILURE;
  }

  ffttBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * NUM_ARCS,
    fftt, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (ffttBuffer)\n";
    return REDP_FAILURE;
  }

  capBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * NUM_ARCS,
    cap, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (capBuffer)\n";
    return REDP_FAILURE;
  }

  optionPointerBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * (NUM_OPTIONS + 1),
    optionPointer, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (optionPointerBuffer)\n";
    return REDP_FAILURE;
  }

  optionsBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * optionArcNum,
    options, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (optionsBuffer)\n";
    return REDP_FAILURE;
  }

  exclusionsBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_OPTIONS * NUM_OPTIONS,
    exclusions, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (exclusionsBuffer)\n";
    return REDP_FAILURE;
  }

  maxCapsBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * optionArcNum,
    maxCaps, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (maxCapsBuffer)\n";
    return REDP_FAILURE;
  }

  costsBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint) * NUM_OPTIONS,
    costs, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (costsBuffer)\n";
    return REDP_FAILURE;
  }

  scenarioProbsBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * NUM_SCENARIOS,
    scenarioProbs, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (scenarioProbsBuffer)\n";
    return REDP_FAILURE;
  }

  /*
  budgetBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint),
    &budget, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (budgetBuffer)\n";
    return REDP_FAILURE;
  }

  sinkNodeBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_uint),
    &sinkNode, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (sinkNodeBuffer)\n";
    return REDP_FAILURE;
  }
  */

  outputBuffer = clCreateBuffer(
    context, 
    CL_MEM_READ_WRITE | CL_MEM_USE_HOST_PTR,
    sizeof(cl_float) * width,
    output, 
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: clCreateBuffer (outputBuffer)\n";
    return REDP_FAILURE;
  }


  /////////////////////////////////////////////////////////////////
  // STEP 6. Building Kernel
  //      6.1 Load CL file, using basic file i/o
  //      6.2 Build CL program object
  //      6.3 Create CL kernel object
  /////////////////////////////////////////////////////////////////
  std::cout <<"Building kernel... \n";
  const char * filename  = "kernelB.cl";
  std::cout << "    Converting " << filename << " to string... \n";
  std::string  sourceStr = convertToString(filename);
  const char * source    = sourceStr.c_str();
  size_t sourceSize[]    = { strlen(source) };

  std::cout <<"    Creating program... \n";
  program = clCreateProgramWithSource(
    context, 
    1, 
    &source,
    sourceSize,
    &status);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Loading Binary into cl_program \
      (clCreateProgramWithBinary)\n";
    return REDP_FAILURE;
  }

  // create a cl program executable for all the devices specified
  std::stringstream buildSS;
  std::string buildStr;
  const char* buildOptions;

  buildSS << "-D NUM_NODES=" << NUM_NODES << " -D NUM_ARCS=" << NUM_ARCS
    << " -D NUM_OPTIONS=" << NUM_OPTIONS << " -D NUM_SCENARIOS=" << NUM_SCENARIOS
    << " -D FLOAT_MAX=INFINITY";
  buildStr = buildSS.str();
  buildOptions = buildStr.c_str();

  std::cout << "    Building program... \n";
  status = clBuildProgram(program, 1, devices, buildOptions, NULL, NULL);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Building Program (clBuildProgram)\n";
    return REDP_FAILURE; 
  }

  // get a kernel object handle for a kernel with the given name
  std::cout << "    Creating kernel from program... \n";
  kernel = clCreateKernel(program, "kernelB", &status);
  if(status != CL_SUCCESS) 
  {  
    std::cout << "Error: Creating Kernel from program. (clCreateKernel)\n";
    return REDP_FAILURE;
  }

  return REDP_SUCCESS;
}


/*
* \brief Run OpenCL program 
*
*        Bind host variables to kernel arguments 
*        Run the CL kernel
*/
int runCLKernels(void)
{
  cl_int status;
  cl_uint maxDims;
  cl_event events[2];
  size_t globalThreads[2];
  size_t localThreads[2];
  size_t maxWorkGroupSize;
  size_t maxWorkItemSizes[3];

  //////////////////////////////////////////////////////////////////// 
  // STEP 7 Analyzing proper workgroup size for the kernel
  //          by querying device information
  //    7.1 Device Info CL_DEVICE_MAX_WORK_GROUP_SIZE
  //    7.2 Device Info CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
  //    7.3 Device Info CL_DEVICE_MAX_WORK_ITEM_SIZES
  //////////////////////////////////////////////////////////////////// 

  std::cout << "Getting workgroup size... \n";

  /**
  * Query device capabilities. Maximum 
  * work item dimensions and the maximmum
  * work item sizes
  */ 
  status = clGetDeviceInfo(
    devices[0], 
    CL_DEVICE_ADDRESS_BITS, 
    sizeof(size_t), 
    (void*)&maxWorkGroupSize, 
    NULL);

  status = clGetDeviceInfo(
    devices[0], 
    CL_DEVICE_MAX_WORK_GROUP_SIZE, 
    sizeof(size_t), 
    (void*)&maxWorkGroupSize, 
    NULL);
  if(status != CL_SUCCESS) 
  {  
    std::cout << "Error: Getting Device Info. (clGetDeviceInfo)\n";
    return REDP_FAILURE;
  }

  status = clGetDeviceInfo(
    devices[0], 
    CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, 
    sizeof(cl_uint), 
    (void*)&maxDims, 
    NULL);
  if(status != CL_SUCCESS) 
  {  
    std::cout << "Error: Getting Device Info. (clGetDeviceInfo)\n";
    return REDP_FAILURE;
  }

  status = clGetDeviceInfo(
    devices[0], 
    CL_DEVICE_MAX_WORK_ITEM_SIZES, 
    sizeof(size_t)*maxDims,
    (void*)maxWorkItemSizes,
    NULL);
  if(status != CL_SUCCESS) 
  {  
    std::cout << "Error: Getting Device Info. (clGetDeviceInfo)\n";
    return REDP_FAILURE;
  }

  globalThreads[0] = width;
  globalThreads[1] = NUM_SCENARIOS;
  localThreads[0]  = 1;
  localThreads[1]  = NUM_SCENARIOS;

  if(localThreads[0] > maxWorkGroupSize ||
    localThreads[0] > maxWorkItemSizes[0] )
  {
    std::cout << "Unsupported: Device does not support requested number of work items.";
    return REDP_FAILURE;
  }

  //////////////////////////////////////////////////////////////////// 
  // STEP 8 Set appropriate arguments to the kernel
  //////////////////////////////////////////////////////////////////// 

  std::cout << "Setting kernel arguments... \n";
  int paramCounter = 0;

  // the output array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&outputBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (output)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the pointer array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&pointerBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (pointer)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the demand array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&demandBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (demand)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the tail array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&tailBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (tail)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the head array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&headBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (head)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the fftt array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&ffttBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (fftt)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // the cap array to the kernel
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&capBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (cap)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // optionPointer parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&optionPointerBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (optionPointer)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // options parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&optionsBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (options)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // exclusions parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&exclusionsBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (exclusions)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // maxCaps parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&maxCapsBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (maxCaps)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // costs parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&costsBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (costs)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // scenarioProbs parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&scenarioProbsBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (scenarioProbs)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // budget parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_uint), 
    &budget);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (budget)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // sinkNode parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_uint), 
    &sinkNode);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (sinkNode)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // riskA parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&riskABuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (riskA)\n";
    return REDP_FAILURE;
  }
  paramCounter++;

  // riskB parameter
  status = clSetKernelArg(
    kernel, 
    paramCounter, 
    sizeof(cl_mem), 
    (void *)&riskBBuffer);
  if(status != CL_SUCCESS) 
  { 
    std::cout << "Error: Setting kernel argument. (risknA)\n";
    return REDP_FAILURE;
  }
  paramCounter++;


  //////////////////////////////////////////////////////////////////// 
  // STEP 9 Enqueue a kernel run call.
  //          Wait till the event completes and release the event
  //////////////////////////////////////////////////////////////////// 

  std::cout << "Enqueuing kernel... \n";
  status = clEnqueueNDRangeKernel(
    commandQueue,
    kernel,
    2,
    NULL,
    globalThreads,
    localThreads,
    0,
    NULL,
    &events[0]);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Enqueuing kernel onto command queue. \
      (clEnqueueNDRangeKernel)\n" << status;
    return REDP_FAILURE;
  }


  /*
  std::cout << "Waiting for kernel... \n";
  // wait for the kernel call to finish execution
  status = clWaitForEvents(1, &events[0]);
  if(status != CL_SUCCESS) 
  { 
  std::cout <<
  "Error: Waiting for kernel run to finish. \
  (clWaitForEvents)\n";
  return REDP_FAILURE;
  }
  */

  status = clReleaseEvent(events[0]);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Release event object. \
      (clReleaseEvent)\n";
    return REDP_FAILURE;
  }

  //////////////////////////////////////////////////////////////////// 
  // STEP 10  Enqueue readBuffer to read the output back
  //  Wait for the event and release the event
  //////////////////////////////////////////////////////////////////// 
  status = clEnqueueReadBuffer(
    commandQueue,
    outputBuffer,
    CL_TRUE,
    0,
    width * sizeof(cl_uint),
    output,
    0,
    NULL,
    &events[1]);

  if(status != CL_SUCCESS) 
  { 
    std::cout << 
      "Error: clEnqueueReadBuffer failed. \
      (clEnqueueReadBuffer)\n";
    return REDP_FAILURE;
  }

  // Wait for the read buffer to finish execution
  status = clWaitForEvents(1, &events[1]);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Waiting for read buffer call to finish. \
      (clWaitForEvents)\n";
    return REDP_FAILURE;
  }

  status = clReleaseEvent(events[1]);
  if(status != CL_SUCCESS) 
  { 
    std::cout <<
      "Error: Release event object. \
      (clReleaseEvent)\n";
    return REDP_FAILURE;
  }
  return REDP_SUCCESS;
}


/*
* \brief Release OpenCL resources (Context, Memory etc.) 
*/
int cleanupCL(void)
{
  cl_int status;

  //////////////////////////////////////////////////////////////////// 
  // STEP 11  CLean up the opencl resources used 
  //////////////////////////////////////////////////////////////////// 

  status = clReleaseKernel(kernel);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseKernel \n";
    return REDP_FAILURE; 
  }
  status = clReleaseProgram(program);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseProgram\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(pointerBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (pointerBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(demandBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (demandBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(riskABuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (riskABuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(riskBBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (riskBBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(tailBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (tailBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(headBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (headBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(ffttBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (ffttBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(capBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (capBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(optionPointerBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (optionPointerBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(optionsBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (optionsBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseMemObject(outputBuffer);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseMemObject (outputBuffer)\n";
    return REDP_FAILURE; 
  }
  status = clReleaseCommandQueue(commandQueue);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseCommandQueue\n";
    return REDP_FAILURE;
  }
  status = clReleaseContext(context);
  if(status != CL_SUCCESS)
  {
    std::cout << "Error: In clReleaseContext\n";
    return REDP_FAILURE;
  }
  return REDP_SUCCESS;
}


/* 
* \brief Releases program's resources 
*/
void cleanupHost(void)
{
  if(pointer != NULL)
  {
    free(pointer);
    pointer = NULL;
  }
  if(demand != NULL)
  {
    free(demand);
    pointer = NULL;
  }
  if(riskA != NULL)
  {
    free(riskA);
    riskA = NULL;
  }
  if(riskB != NULL)
  {
    free(riskB);
    riskB = NULL;
  }
  if(tail != NULL)
  {
    free(tail);
    tail = NULL;
  }
  if(head != NULL)
  {
    free(head);
    head = NULL;
  }
  if(fftt != NULL)
  {
    free(fftt);
    fftt = NULL;
  }
  if(cap != NULL)
  {
    free(cap);
    cap = NULL;
  }
  if(optionPointer != NULL)
  {
    free(optionPointer);
    optionPointer = NULL;
  }
  if(options != NULL)
  {
    free(options);
    options = NULL;
  }
  if(exclusions != NULL)
  {
    free(exclusions);
    exclusions = NULL;
  }
  if(maxCaps != NULL)
  {
    free(maxCaps);
    maxCaps = NULL;
  }
  if(costs != NULL)
  {
    free(costs);
    costs = NULL;
  }
  if(scenarioProbs != NULL)
  {
    free(scenarioProbs);
    scenarioProbs = NULL;
  }
  if(output != NULL)
  {
    free(output);
    output = NULL;
  }
  if(devices != NULL)
  {
    free(devices);
    devices = NULL;
  }
}


void print1DArray(
  const std::string arrayName, 
  const float * arrayData, 
  const unsigned int length)
{
  cl_uint i;
  cl_uint numElementsToPrint = (128 < length) ? 128 : length;

  std::cout << std::endl;
  std::cout << arrayName << ":" << std::endl;
  for(i = 0; i < numElementsToPrint; ++i)
  {
    if(i > 0 && i % width == 0){
      std::cout << "\n\n";
    }
    std::cout << i << " " << arrayData[i] << std::endl;
  }
  std::cout << std::endl;

}

/*
void verify()
{

bool passed = true;
for(unsigned int i = 0; i < width; ++i)
if(input[i] * sinkNode != output[i])
passed = false;

if(passed == true)
std::cout << "Passed!\n" << std::endl;
else
std::cout << "Failed!\n" << std::endl;
}
*/

int main(int argc, char * argv[])
{
  std::ifstream networkFile;
  std::string line;
  std::string networkName;
  std::stringstream lineStream(std::ios_base::ate | std::ios_base::in | std::ios_base::out);
  cl_uint firstNode;
  cl_uint firstArc;
  cl_uint firstOption;
  cl_uint firstExclusion;
  cl_uint firstCap;
  cl_uint firstCost;
  cl_uint firstScenarioProb;
  cl_uint firstScenarioRiskA;
  cl_uint firstScenarioRiskB;
  cl_uint budgetLine;
  cl_uint sinkNodeLine;
  cl_uint i;
  cl_uint lineCount = 1;
  cl_uint *tmpOptions;

  if ( argc == 2 ) {
    std::cout << "\nReading "<< argv[1] << "\n";
    networkFile.open(argv[1], std::ios::in);
    if(networkFile.is_open()) {

      // get header and line numbers
      getline(networkFile, line);
      ++lineCount;
      lineStream.clear();
      lineStream.str(line);
      lineStream >> networkName >> firstNode >> firstArc >> firstOption >> firstExclusion
        >> firstCap >> firstCost >> firstScenarioProb >> firstScenarioRiskA >> firstScenarioRiskB 
        >> budgetLine >> sinkNodeLine;

      std::cout << "Parsing " << networkName << std::endl;

      NUM_NODES = firstArc - firstNode;
      NUM_ARCS = firstOption - firstArc;
      NUM_OPTIONS = firstExclusion - firstOption;
      NUM_SCENARIOS = firstScenarioRiskA - firstScenarioProb;

      // get and ignore comments
      for(i = 0; i < firstNode - 2; ++i) {
        getline(networkFile, line);
        ++lineCount;
      }

      std::cout << networkName << " has" << std::endl
        << "    " << NUM_NODES << " nodes" << std::endl
        << "    " << NUM_ARCS << " links" << std::endl
        << "    " << NUM_OPTIONS << " options" << std::endl
        << "    " << NUM_SCENARIOS << " scenarios" << std::endl;

      // get nodes and populate node arrays
      if(firstNode - lineCount != 0) {
        std::cout << "\nError reading network file header";
        return REDP_FAILURE;
      }
      cl_uint nodeIntBytes = NUM_NODES * sizeof(cl_uint);
      cl_uint nodeFloatBytes = NUM_NODES * sizeof(cl_float);
      pointer = (cl_uint *) malloc(nodeIntBytes);
      if(!pointer)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      demand = (cl_uint *) malloc(nodeIntBytes);
      if(!demand)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      riskA = (cl_float *) malloc(nodeFloatBytes);
      if(!riskA)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      riskB = (cl_float *) malloc(nodeFloatBytes);
      if(!riskB)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineCount < firstArc) {
        getline(networkFile, line);
        ++lineCount;
        lineStream.clear();
        lineStream.str(line);
        while(lineStream >> pointer[i] >> demand[i]) {
          std::cout << i + 1 << " " << pointer[i] << " " << demand[i] << "\n";
        }
        ++i;
      }

      // get arcs and populate arc arrays
      lineStream.clear();
      while(lineCount < firstOption) {
        getline(networkFile, line);
        ++lineCount;
        lineStream << line << " ";
      }
      cl_uint arcIntBytes = NUM_ARCS * sizeof(cl_uint);
      cl_uint arcFloatBytes = NUM_ARCS * sizeof(cl_float);
      tail = (cl_uint *) malloc(arcIntBytes);
      if(!tail)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      head = (cl_uint *) malloc(arcIntBytes);
      if(!head)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      fftt = (cl_float *) malloc(arcFloatBytes);
      if(!fftt)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      cap = (cl_float *) malloc(arcFloatBytes);
      if(!cap)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineStream >> tail[i] >> head[i] >> fftt[i] >> cap[i]) {
        std::cout << i + 1 << " " << tail[i] << " " << head[i] << " " << fftt[i] << " " << cap[i] << "\n";
        i++;
      }

      // get options
      cl_uint optionIntBytes = NUM_OPTIONS * sizeof(cl_uint);
      optionPointer = (cl_uint *) malloc(optionIntBytes + sizeof(cl_uint));
      if(!optionPointer) {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      tmpOptions = (cl_uint *) malloc(arcIntBytes);
      if(!tmpOptions) {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineCount < firstExclusion) {
        getline(networkFile, line);
        ++lineCount;
        lineStream.clear();
        lineStream.str(line);
        optionPointer[lineCount - firstOption - 1] = i;
        std::cout << optionPointer[lineCount - firstOption - 1] << " ";
        while(lineStream >> tmpOptions[i]) {
          ++i;
        }
      }
      optionPointer[NUM_OPTIONS] = i;
            std::cout << optionPointer[NUM_OPTIONS] << std::endl;

      optionArcNum = i;
      options = (cl_uint *) malloc(optionArcNum * sizeof(cl_uint));
      for(i = 0; i < optionArcNum; ++i) {
        options[i] = tmpOptions[i];
        std::cout << options[i] << " ";
      }
      std::cout << std::endl;
      /*
      for(i = 0; i < NUM_NODES; ++i) {
        std::cout << optionPointer[i % NUM_OPTIONS] << " ";
      }
      std::cout << std::endl;
      */

      // get exclutions
      cl_uint nextExclusion;
      exclusions = (cl_bool *) malloc(NUM_OPTIONS * NUM_OPTIONS * sizeof(cl_bool));
      if(!exclusions) {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      for(i = 0; i < NUM_OPTIONS * NUM_OPTIONS; ++i) {
        exclusions[i] = 1;
      }
      i = 0;
      while(lineCount < firstCap) {
        getline(networkFile, line);
        ++lineCount;
        lineStream.clear();
        lineStream.str(line);
        while(lineStream >> nextExclusion && nextExclusion > 0) {
          exclusions[(lineCount - firstExclusion - 1) * NUM_OPTIONS + nextExclusion - 1] = 0;
        }
      }
      for(i = 0; i < NUM_OPTIONS * NUM_OPTIONS; ++i) {
        if(i % NUM_OPTIONS == 0 && i > 0) {std::cout << std::endl;}
        std::cout << exclusions[i] << " ";
      }
      std::cout << std::endl;

      // get shelter caps
      maxCaps = (cl_uint *) malloc(optionIntBytes);
      if(!maxCaps)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineCount < firstCost) {
        lineStream.clear();
        getline(networkFile, line);
        ++lineCount;
        lineStream << line;
        lineStream >> maxCaps[i];
        std::cout << maxCaps[i] << std::endl;
        ++i;
      }

      // get option costs
      costs = (cl_uint *) malloc(optionIntBytes);
      if(!costs)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineCount < firstScenarioProb) {
        lineStream.clear();
        getline(networkFile, line);
        ++lineCount;
        lineStream << line;
        lineStream >> costs[i];
        std::cout << costs[i] << std::endl;
        ++i;
      }

      // get scenario probalities
      cl_uint scenarioFloatBytes = NUM_SCENARIOS * sizeof(cl_float);
      scenarioProbs = (cl_float *) malloc(scenarioFloatBytes);
      if(!costs)
      {
        std::cout << "Error: Failed to allocate input memory on host\n";
        return REDP_FAILURE;
      }
      i = 0;
      while(lineCount < firstScenarioRiskA) {
        lineStream.clear();
        getline(networkFile, line);
        ++lineCount;
        lineStream << line;
        lineStream >> scenarioProbs[i];
        std::cout << scenarioProbs[i] << std::endl;
        ++i;
      }

      // get risk A
      riskA = (cl_float *) malloc(NUM_SCENARIOS * nodeFloatBytes);
      i = 0;
      while(lineCount < firstScenarioRiskB) {
        getline(networkFile, line);
        ++lineCount;
        lineStream.clear();
        lineStream.str(line);
        while(lineStream >> riskA[i]) {
          std::cout << riskA[i] << " ";
          ++i;
        }
        std::cout << std::endl;
      }

      // get risk B
      riskB = (cl_float *) malloc(NUM_SCENARIOS * nodeFloatBytes);
      i = 0;
      while(lineCount < budgetLine) {
        getline(networkFile, line);
        ++lineCount;
        lineStream.clear();
        lineStream.str(line);
        while(lineStream >> riskB[i]) {
          std::cout << riskB[i] << " ";
          ++i;
        }
        std::cout << std::endl;
      }

      // get budget cap
      getline(networkFile, line);
      ++lineCount;
      lineStream.clear();
      lineStream.str(line);
      lineStream >> budget;
      std::cout << budget << std::endl;

      // get sink node
      getline(networkFile, line);
      ++lineCount;
      lineStream.clear();
      lineStream.str(line);
      lineStream >> sinkNode;
      std::cout << sinkNode << std::endl;

      networkFile.close();
    }
  }

  // Initialize Host application 
  if(initializeHost() != REDP_SUCCESS)
    return REDP_FAILURE;

  // Initialize OpenCL resources
  if(initializeCL() != REDP_SUCCESS)
    return REDP_FAILURE;

  // Run the CL program
  if(runCLKernels() != REDP_SUCCESS)
    return REDP_FAILURE;

  // Print output array
  print1DArray(std::string("Output"), output, width);

  // Releases OpenCL resources 
  if(cleanupCL()!= REDP_SUCCESS)
    return REDP_FAILURE;

  // Release host resources
  cleanupHost();

  return REDP_SUCCESS;
}