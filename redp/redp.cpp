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
#include "kernelB.h"
#include "time.h"


/*
* Allocate and initialize memory 
* on the host. Print input array. 
*/
int initializeHost(void)
{
  width               = (int) pow(2, NUM_OPTIONS);
  //width               = (int) pow(2, 15);
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

void print1DArray(
  const std::string arrayName, 
  const unsigned int * arrayData, 
  const unsigned int length)
{
  cl_uint i;
  cl_uint numElementsToPrint = (10000 < length) ? 10000 : length;

  std::cout << std::endl;
  std::cout << arrayName << ":" << std::endl;
  for(i = 0; i < numElementsToPrint; ++i)
  {
    if(i % NUM_NODES == 0){
      std::cout << "\n";
    }
    std::cout << arrayData[i] << " ";
  }
  std::cout << std::endl;

}

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
  cl_uint j;
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

      /*
      NUM_NODES = firstArc - firstNode;
      NUM_ARCS = firstOption - firstArc;
      NUM_OPTIONS = firstExclusion - firstOption;
      NUM_SCENARIOS = firstScenarioRiskA - firstScenarioProb;
      */

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
          std::cout << i << " " << pointer[i] << " " << demand[i] << "\n";
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
        fftt[i] = fftt[i] == 999999 ? CL_INFINITY : fftt[i];
        std::cout << i << " " << tail[i] << " " << head[i] << " " << fftt[i] << " " << cap[i] << "\n";
        i++;
      }

      // get options
      cl_uint optionIntBytes = NUM_OPTIONS * sizeof(cl_uint);
      optionPointer = (cl_uint *) malloc(optionIntBytes);
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
      optionPointer[NUM_OPTIONS] = i; // end of option array
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
        if(i % 6 == 0 && i > 0) {std::cout << std::endl;}
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

  initializeHost();

  clock_t startTime = clock();

  int decisionId;
  int scenarioId;
  float expectedRisk;

  for(decisionId = 0; decisionId < width; ++decisionId) {
    expectedRisk = 0;
    for(scenarioId = 0; scenarioId < NUM_SCENARIOS; ++scenarioId) {

      // set optionSwitch for this decision
      if(scenarioId == 0) { 
        std::cout << decisionId << ", ";
        for(i = 0; i < NUM_OPTIONS; ++i) {
          optionSwitch[NUM_OPTIONS - 1 - i] = (decisionId >> i) & 1;
          std::cout << optionSwitch[NUM_OPTIONS - 1 - i];
        }
      }

      // check exclusions for this decision
      for(i = 0; i < NUM_OPTIONS; ++i) {
        for(j = i; j < NUM_OPTIONS; ++j) {
          if(optionSwitch[i] && optionSwitch[j] && !exclusions[i * NUM_OPTIONS + j]) {
            if(scenarioId == 0) {
              output[decisionId] = -1;
              std::cout << ", exclusion"; 
            }
          }
        }
      }

      // check budget for this decision
      int totalOptionCost = 0;
      for(i = 0; i < NUM_OPTIONS; ++i) {
        totalOptionCost += optionSwitch[i] * costs[i];
      }
      if(totalOptionCost > budget) {
        if(scenarioId == 0) {
          output[decisionId] = -2;
          std::cout << ", budget";
        }
      }

      // set risk for this decision and scenario
      for(i = 0; i < NUM_ARCS; ++i) {
        int headRiskA = riskA[head[i]-1];
        int tailRiskA = riskA[tail[i]-1];
        int headRiskB = riskB[head[i]-1];
        int tailRiskB = riskB[tail[i]-1];
        arcRiskA[i] = headRiskA > tailRiskB ? headRiskA : tailRiskA;
        arcRiskB[i] = headRiskB > tailRiskB ? headRiskB : tailRiskB;
        if(!getMaxCap(i, options, optionPointer)) {
          arcRiskB[i] = 0.01;
        }
      }

      // set network parameters for decision
      for(i = 0; i < NUM_OPTIONS; ++i) {
        for(j = optionPointer[i]; j < optionPointer[NUM_OPTIONS] && j < optionPointer[i+1]; ++j) {
          if(optionSwitch[i] == 0) {
            // turn off arc options[j]
            arcRiskB[options[j]] = CL_INFINITY;
          }
        }
      }

      initNetwork(sinkNode, topoOrder, tail, head, pointer, demand, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow);
      // printStatus();
      // if(decisionId == 63) {printStatus();}
      equilibrateNetwork(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
      std::cout << " " << verifyEq(minCost,maxCost, demand, sinkNode);
      // printStatus();
      // if(decisionId == 28) {printStatus();}

      // calculate and return maximum risk for this decision
      maxRisk[scenarioId] = 0;
      for(i = 0; i < NUM_NODES; ++i) {
        if(minCost[i] > maxRisk[scenarioId] && flow[nextMinArc[i]] > 0) {
          maxRisk[scenarioId] = minCost[i];
        }
      }
      expectedRisk += maxRisk[scenarioId] * scenarioProbs[scenarioId];
      // output[decisionId] = optionSwitch[decisionId % NUM_OPTIONS];
      // output[decisionId] = decisionId;
      // output[decisionId] = minCost[decisionId % 9];
      // output[decisionId] = scenarioProbs[decisionId % NUM_SCENARIOS];
      // output[decisionId] = tmpscenarioProbs[decisionId % NUM_SCENARIOS];
      // output[decisionId] = optionPointer[decisionId % NUM_OPTIONS];
    }
    output[decisionId] = expectedRisk;
    if(scenarioId == 4) {
      std::cout << ", " << output[decisionId] << std::endl;
    }
  }
  std::cout << clock()-startTime;
  return REDP_SUCCESS;
}