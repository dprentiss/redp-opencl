// uncomment to allow CodeXL kernel analysis
#ifndef NUM_NODES
#define NUM_NODES 14
#endif

#ifndef NUM_ARCS
#define NUM_ARCS 53
#endif

#ifndef FLOAT_MAX
#define FLOAT_MAX INFINITY
#endif

#ifndef NUM_OPTIONS
#define NUM_OPTIONS 6
#endif

#ifndef NUM_SCENARIOS
#define NUM_SCENARIOS 4
#endif
/*
*/


float setCost(uint arc, global float *fftt, float *flow, global float *cap, float *arcRiskA, float *arcRiskB)
{
  return arcRiskA[arc] * fftt[arc] 
    * (1.0F + 0.15F*pow((flow[arc]/cap[arc]), 4.0F)) + arcRiskB[arc];
}

float calcCost(uint arc,
               float deltaV,
               global float *fftt,
               float *flow,
               global float *cap, 
               float *arcRiskA,
               float *arcRiskB)
{
  return arcRiskA[arc] * fftt[arc] 
    * (1.0F + 0.15F*pow(((flow[arc]+deltaV)/cap[arc]), 4.0F)) 
    + arcRiskB[arc];
}

float calcCostPrime(uint arc, float deltaV, global float *fftt, float *flow, global float *cap, float *arcRiskA, float *arcRiskB)
{
  return 0.6F * arcRiskA[arc] * fftt[arc] * pow(flow[arc] + deltaV, 3.0F) / pow(cap[arc], 4.0F);
}

float getMaxCap(int arc,
                global uint *options,
                global uint *optionPointer,
                global uint *maxCaps,
                global uint *demand,
                uint sinkNode,
                global uint *tail,
                global uint *head)
{
  int i;
  for(i = 0; i < NUM_OPTIONS; ++i) {
    if(arc == options[optionPointer[i]]) {
      return maxCaps[i];
    }
  }
  if(demand[tail[arc]] && head[arc] == sinkNode) {
    return demand[tail[arc]];
  }
  return 0;
}

int isInNodeArray(uint *array, uint n)
{
  uint i;
  for(i = 0; i < NUM_NODES; ++i) {
    if (array[i] == n) {
      return 1;
    }
  }
  return 0;
}

void updateCost(float *minCost,
                float *maxCost,
                float *cost,
                float *flow,
                uint sinkNode,
                global uint *pointer,
                global uint *tail,
                global uint *head,
                uint *topoOrder,
                uint *nextMinArc,
                uint *nextMaxArc,
                global uint *options,
                global uint *optionPointer,
                global uint *maxCaps,
                global uint *demand)
{
  uint i;
  uint arc;
  uint hasCap;
  float costCand;

  for(i = 0; i < NUM_NODES; i++) {
    minCost[i] = FLOAT_MAX;
    maxCost[i] = 0;
  }

  minCost[sinkNode-1] = 0;
  nextMinArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    arc = pointer[topoOrder[i]-1]-1;
    while(head[arc] == topoOrder[i]) {
      costCand = cost[arc] + minCost[head[arc]-1];
      if(costCand < minCost[tail[arc]-1]) {
        if(hasCap = getMaxCap(arc, options, optionPointer, maxCaps, demand, sinkNode, tail, head)) {
          if(flow[arc] >= hasCap) {
            ++arc;
            continue;
          }
        }
        minCost[tail[arc]-1] = costCand;
        nextMinArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
  }

  for(i = 0; i < NUM_NODES; i++) {
    maxCost[i] = minCost[i];
  }

  maxCost[sinkNode-1] = 0;
  nextMaxArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    arc = pointer[topoOrder[i]-1]-1;
    while(head[arc] == topoOrder[i]) {
      if(cost[arc] >= FLOAT_MAX) {
        ++arc;
        continue;
      }
      costCand = cost[arc] + maxCost[head[arc]-1];
      if(costCand >= maxCost[tail[arc]-1] && flow[arc] > 0) {
        maxCost[tail[arc]-1] = costCand;
        nextMaxArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
  }
}

void getPaths(uint origin,
              uint sinkNode,
              global uint *head,
              uint *minPath,
              uint *maxPath,
              uint *nextMinArc,
              uint *nextMaxArc)
{
  uint i, j, k;
  uint length = 0;
  uint hasCommonNode = 1;

  minPath[0] = origin;
  maxPath[0] = origin;

  while(hasCommonNode) {
    if(minPath[0] == sinkNode - 1 || maxPath[0] == sinkNode - 1) {
      return;
    } else if(nextMaxArc[minPath[0]] == nextMinArc[maxPath[0]]) {
      minPath[0] = head[nextMaxArc[minPath[0]]]-1; 
      maxPath[0] = minPath[0]; 
      continue;
    }
    hasCommonNode = 0;
  }

  while(!hasCommonNode) {
    ++length;
    if(minPath[length - 1] != sinkNode-1 && minPath[length - 1] != NUM_NODES) {
      minPath[length] = head[nextMinArc[minPath[length-1]]]-1;
    } else {
      minPath[length] = NUM_NODES;
    }
    if(maxPath[length - 1] != sinkNode-1 && maxPath[length - 1] != NUM_NODES) {
      maxPath[length] = head[nextMaxArc[maxPath[length-1]]]-1;
    } else {
      maxPath[length] = NUM_NODES;
    }
    for(i = 1; i < length + 1; ++i) {
      for(j = 1; j < length + 1; ++j) {
        if(maxPath[i] == minPath[j]) {
          for(k = j + 1; k < NUM_NODES; ++k) {
            minPath[k] = NUM_NODES;
          }
          for(k = i + 1; k < NUM_NODES; ++k) {
            maxPath[k] = NUM_NODES;
          }
          hasCommonNode = 1;
          return;
        }
      }
    }
  }
}

void shiftFlow(uint sinkNode,
               uint *maxPath,
               uint *minPath,
               global uint *pointer,
               global uint *tail,
               global uint *head,
               uint *topoOrder,
               uint *nextMinArc,
               uint *nextMaxArc,
               float *flow,
               global float *fftt,
               global float *cap,
               float *cost,
               float *minCost,
               float *maxCost,
               float *arcRiskA,
               float *arcRiskB,
               global uint *options,
               global uint *optionPointer,
               global uint *maxCaps,
               global uint *demand)
{
  uint i = 0;
  float nextFlow;
  float max = FLOAT_MAX;
  float shelterCap;
  float root = 0;
  float f = 1.0F;
  float fPrime;
  uint maxIter = 10;
  uint iter = 0;

  // find minimum flow on max path and limit max flow shift
  for(i = 0; i < NUM_NODES; ++i) {
    if(nextMaxArc[maxPath[i]] >= NUM_ARCS
      || maxPath[i] >= NUM_NODES) {
        break;
    }
    nextFlow = flow[nextMaxArc[maxPath[i]]];
    if(nextFlow < max) {
      max = nextFlow;
    }
  }

  // determine if max flow shift is limited by a shelter arc
  for(i = 0; i < NUM_NODES; ++i) {
    if(nextMinArc[minPath[i]] >= NUM_ARCS
      || minPath[i] >= NUM_NODES) {
        break;
    }
    if(shelterCap = getMaxCap(nextMinArc[minPath[i]], options, optionPointer, maxCaps, demand, sinkNode, tail, head)) {
      max = shelterCap - flow[nextMinArc[minPath[i]]];
    }
  }

  root = max;

  while(fabs(f) > 0.0001F) {
    if(iter >= maxIter) {
      break;
    }
    ++iter;
    i = 0;
    f = 0; 
    fPrime = 0; 
    while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
      if(minPath[i+1] < NUM_NODES) {
        f += calcCost(nextMinArc[minPath[i]], root, fftt, flow, cap, arcRiskA, arcRiskB);
        fPrime += calcCostPrime(nextMinArc[minPath[i]], root, fftt, flow, cap, arcRiskA, arcRiskB);
      }
      if(maxPath[i+1] < NUM_NODES) {
        f -= calcCost(nextMaxArc[maxPath[i]], -1*root, fftt, flow, cap, arcRiskA, arcRiskB);
        fPrime += calcCostPrime(nextMaxArc[maxPath[i]], -1*root, fftt, flow, cap, arcRiskA, arcRiskB);
      }
      ++i;
    }
    root -= f/fPrime;
  }
  if(root > max) {
    root = max;
  } else if(root < 0.0001f) {
    root = 0.1 * max;
    //return;
  } 
  i = 0;
  while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
    if(minPath[i+1] < NUM_NODES) {
      flow[nextMinArc[minPath[i]]] += root;
      cost[nextMinArc[minPath[i]]] = setCost(nextMinArc[minPath[i]], fftt, flow, cap, arcRiskA, arcRiskB);
    }
    if(maxPath[i+1] < NUM_NODES) {
      flow[nextMaxArc[maxPath[i]]] -= root;
      cost[nextMaxArc[maxPath[i]]] = setCost(nextMaxArc[maxPath[i]], fftt, flow, cap, arcRiskA, arcRiskB);
    }
    ++i;
  }
  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, options, optionPointer, maxCaps, demand);
}

void equilibrateBush(uint sinkNode, 
                     uint *maxPath,
                     uint *minPath,
                     uint *nextMinArc,
                     uint *nextMaxArc,
                     global uint *pointer,
                     global uint *demand,
                     global uint *tail,
                     global uint *head,
                     uint *topoOrder,
                     float *flow,
                     global float *fftt,
                     global float *cap,
                     float *cost,
                     float *minCost,
                     float *maxCost,
                     float *arcRiskA,
                     float *arcRiskB, 
                     global uint *options,
                     global uint *optionPointer,
                     global uint *maxCaps)
{
  int i;
  int j;
  int eqBush = 0;
  int updated;
  int iter = 0;
  float tol = 0.001f;

  while(!eqBush) {
    updated = 0;
    if(iter > 0 && iter % 100 == 0) {
      tol = tol * 10.f;
    }
    for(i = 0; i < NUM_NODES; ++i) {
      j = topoOrder[i] - 1;
      if(fabs(maxCost[j] - minCost[j]) > tol && demand[j] > 0) {
        getPaths(j, sinkNode, head, minPath, maxPath, nextMinArc, nextMaxArc);
        shiftFlow(sinkNode, maxPath, minPath, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, flow, fftt, cap, cost, minCost, maxCost, arcRiskA, arcRiskB, options, optionPointer, maxCaps, demand);
        updated = 1;
      }
    }
    if(updated == 0) {
      eqBush = 1;
    }
    iter++;
  }
}

void initNetwork(const uint sinkNode,
                 uint *topoOrder,
                 global uint *tail,
                 global uint *head,
                 global uint *pointer,
                 global uint *demand,
                 uint *nextMinArc,
                 uint *nextMaxArc,
                 float *minCost,
                 float *maxCost,
                 float *cost,
                 global float *fftt,
                 global float *cap,
                 float *flow,
                 float *arcRiskA,
                 float *arcRiskB,
                 global uint *options,
                 global uint *optionPointer,
                 global uint *maxCaps)
{
  uint i;
  uint lastI;
  uint j;
  uint k;
  uint l;
  uint tmpNode;
  uint arc;
  float costCand;
  uint nextArc;

  // set ititial topoOrder and min costs
  for(i = 0; i < NUM_NODES; ++i) {
    topoOrder[i] = 0;
    minCost[i] = FLOAT_MAX;
  }
  for(i = 0; i < NUM_ARCS; ++i) {
    cost[i] = FLOAT_MAX;
    flow[i] = 0;
  }
  minCost[sinkNode-1] = 0;
  nextMinArc[sinkNode-1] = NUM_ARCS;
  topoOrder[0] = sinkNode;
  i = 1;
  j = 0;
  while(j < NUM_NODES) {
    lastI = i;
    arc = pointer[topoOrder[j]-1]-1;
    while(head[arc] == topoOrder[j]) {
      // update topoOrder
      if(!isInNodeArray(topoOrder, tail[arc])) {
        topoOrder[i] = tail[arc];
        ++i;
      }
      // propagate mincosts
      costCand = setCost(arc, fftt, flow, cap, arcRiskA, arcRiskB) + minCost[head[arc]-1];
      if(costCand <= minCost[tail[arc]-1] && !getMaxCap(arc, options, optionPointer, maxCaps, demand, sinkNode, tail, head)) {
        minCost[tail[arc]-1] = costCand;
        nextMinArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
    //printStatus();
    // sort new nodes
    for(k = lastI; k < i; ++k) {
      for(l = k + 1; l < i; ++l) {
        if(minCost[l - 1] < minCost[k - 1]) {
          tmpNode = topoOrder[k];
          topoOrder[k] = topoOrder[l];
          topoOrder[l] = tmpNode;
          //printStatus();
        }
      }
    }
    ++j;
  }

  // turn on initial bush
  for(i = 0; i < NUM_ARCS; i++) {
    if(minCost[head[i]-1] == minCost[tail[i]-1]) {
      minCost[tail[i]-1] += 0.00001;
    }
    if(minCost[head[i]-1] < minCost[tail[i]-1]) {
      cost[i] = setCost(i, fftt, flow, cap, arcRiskA, arcRiskB);
    }
  }
  //printStatus();

  // init flows on bush
  nextArc = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    if(demand[i] > 0) {
      nextArc = nextMinArc[i];
      while(nextArc < NUM_ARCS) {
        flow[nextArc] += demand[i];
        nextArc = nextMinArc[head[nextArc]-1];
        // printStatus();
      }
    }
  }

  // init costs on bush
  for(arc = 0; arc < NUM_ARCS; ++arc) {
    if(minCost[head[arc]-1] < minCost[tail[arc]-1] || flow[arc] > 0) {
      cost[arc] = setCost(arc, fftt, flow, cap, arcRiskA, arcRiskB);
    }
  }
  // printStatus();

  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, options, optionPointer, maxCaps, demand);
  // printStatus();
}

int updateBush(uint sinkNode,
               uint *topoOrder,
               global uint tail[],
               global uint head[],
               global uint pointer[],
               uint nextMinArc[],
               uint nextMaxArc[],
               float minCost[],
               float maxCost[],
               float cost[],
               global float *fftt,
               global float *cap,
               float flow[],
               float *arcRiskA,
               float *arcRiskB,
               global uint *options,
               global uint *optionPointer,
               global uint *maxCaps,
               global uint *demand)

{
  uint i;
  uint j;
  uint arc;
  uint isUpdated = 0;

  for(i = 0; i < NUM_ARCS; ++i) {
    if(flow[i] == 0 && cost[i] < FLOAT_MAX) {
      j = pointer[tail[i]-1]-1;
      while(head[j] == tail[i]) {
        if(tail[j] == head[i] && fftt[j] + minCost[head[j]-1] < minCost[tail[j]-1]) {
          cost[i] = FLOAT_MAX;
          cost[j] = setCost(j, fftt, flow, cap, arcRiskA, arcRiskB); 
          isUpdated = 1;
        }
        ++j;
      }
    }
  }

  // set topoOrder
  for(i = 0; i < NUM_NODES; ++i) {
    topoOrder[i] = 0;
  }
  topoOrder[0] = sinkNode;
  i = 1;
  j = 0;
  while(i < NUM_NODES) {
    arc = pointer[topoOrder[j]-1]-1;
    while(head[arc] == topoOrder[j]) {
      // update topoOrder
      if(!isInNodeArray(topoOrder, tail[arc])) {
        topoOrder[i] = tail[arc];
        ++i;
      }
      ++arc;
    }
    ++j;
  }

  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, options, optionPointer, maxCaps, demand);
  return isUpdated;
}

void equilibrateNetwork(uint sinkNode,
                        uint *maxPath,
                        uint *minPath,
                        uint *nextMinArc,
                        uint *nextMaxArc,
                        global uint *pointer,
                        global uint *demand,
                        global uint *tail,
                        global uint *head,
                        uint *topoOrder,
                        float *flow,
                        global float *fftt,
                        global float *cap,
                        float *cost,
                        float *minCost,
                        float *maxCost,
                        float *arcRiskA,
                        float *arcRiskB,
                 global uint *options,
                 global uint *optionPointer,
                 global uint *maxCaps)
{
  uint bushUpdated = 1;
  while(bushUpdated) {
    equilibrateBush(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost, arcRiskA, arcRiskB, options, optionPointer, maxCaps);
    bushUpdated = updateBush(sinkNode, topoOrder, tail, head, pointer, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow, arcRiskA, arcRiskB, options, optionPointer, maxCaps, demand); 
  }
}

float verify(uint sinkNode,
             uint *maxPath,
             uint *minPath,
             uint *nextMinArc,
             uint *nextMaxArc,
             global uint *pointer,
             global uint *demand,
             global uint *tail,
             global uint *head,
             uint *topoOrder,
             float *flow,
             global float *fftt,
             global float *cap,
             float *cost,
             float *minCost,
             float *maxCost,
             float *arcRiskA,
             float *arcRiskB,
             global uint *options,
             global uint *optionPointer,
             global uint *maxCaps)
{
  // demand is met
  // trip conservation
  // non-negativity
  // unsaturated links min cost ~ max costs
  // saturated links < min cost
  // capacities honored
  // no flow on prohibited links
  return 0;
}



__kernel void kernelB(global float * output,
                      global unsigned int * pointer,
                      global unsigned int * demand,
                      global unsigned int * tail,
                      global unsigned int * head,
                      global float * fftt,
                      global float * cap,
                      global unsigned int * optionPointer,
                      global unsigned int * options,
                      global unsigned int * exclusions,
                      global unsigned int * maxCaps,
                      global unsigned int * costs,
                      global float * scenarioProbs,
                      unsigned int budget,
                      unsigned int sinkNode,
                      global float * riskA,
                      global float * riskB)

{
  // private
  const uint decisionIx = get_global_id(0);
  const uint decisionId = decisionIx % 64;
  const uint scenarioId = get_global_id(1);
  uint i;
  uint j;

  // node variables
  float minCost[NUM_NODES];
  float maxCost[NUM_NODES];
  uint nextMaxArc[NUM_NODES];
  uint nextMinArc[NUM_NODES];
  uint minPath[NUM_NODES];
  uint maxPath[NUM_NODES];
  uint topoOrder[NUM_NODES];

  // arc variables
  float flow[NUM_ARCS];
  float cost[NUM_ARCS];
  float arcRiskA[NUM_ARCS];
  float arcRiskB[NUM_ARCS];

  // local (shared by all scenarios for this decision)
  local float maxRisk[NUM_SCENARIOS];
  local bool optionSwitch[NUM_OPTIONS];

  // set optionSwitch for this decision
  if(scenarioId == 0) {
    for(i = 0; i < NUM_OPTIONS; ++i) {
      optionSwitch[NUM_OPTIONS - 1 - i] = (decisionId >> i) & 1;
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);

  // check exclusions for this decision
  for(i = 0; i < NUM_OPTIONS; ++i) {
    for(j = 0; j < NUM_OPTIONS; ++j) {
      if(optionSwitch[i] && optionSwitch[j] && !exclusions[i * NUM_OPTIONS + j]) {
        if(scenarioId == 0) {output[decisionIx] = -1;}
        return;
      }
    }
  }

  // check budget for this decision
  uint totalOptionCost = 0;
  for(i = 0; i < NUM_OPTIONS; ++i) {
    totalOptionCost += optionSwitch[i] * costs[i];
  }
  if(totalOptionCost > budget) {
    if(scenarioId == 0) {output[decisionIx] = -2;}
    return;
  }

  // set risk for this decision and scenario
  for(i = 0; i < NUM_ARCS; ++i) {
    uint headRiskA = riskA[head[i]-1];
    uint tailRiskA = riskA[tail[i]-1];
    uint headRiskB = riskB[head[i]-1];
    uint tailRiskB = riskB[tail[i]-1];
    arcRiskA[i] = headRiskA > tailRiskA ? headRiskA : tailRiskA;
    arcRiskB[i] = headRiskB > tailRiskB ? headRiskB : tailRiskB;
    if(!getMaxCap(i, options, optionPointer, maxCaps, demand, sinkNode, tail, head)) {
      arcRiskB[i] = 0.001;
    }
  }

  // set network parameters for decision
  for(i = 0; i < NUM_OPTIONS; ++i) {
    for(j = optionPointer[i]; j < optionPointer[NUM_OPTIONS] && j < optionPointer[i+1]; ++j) {
      if(optionSwitch[i] == 0) {
        // turn off arc options[j]
        arcRiskB[options[j]] = INFINITY;
      }
    }
  }

  // find equilibrium for this decision and scenario
  initNetwork(sinkNode, topoOrder, tail, head, pointer, demand, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow, arcRiskA, arcRiskB, options, optionPointer, maxCaps);
  equilibrateNetwork(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost, arcRiskA, arcRiskB, options, optionPointer, maxCaps);

  // calculate and return maximum risk for this decision
  maxRisk[scenarioId] = 0;
  for(i = 0; i < NUM_NODES; ++i) {
    if(minCost[i] > maxRisk[scenarioId] && flow[nextMinArc[i]] > 0) {
      maxRisk[scenarioId] = minCost[i];
    }
  }
  barrier(CLK_LOCAL_MEM_FENCE);
  if(scenarioId == 0) {
    float expectedRisk = 0;
    for(i = 0; i < NUM_SCENARIOS; ++i) {
      expectedRisk += maxRisk[i] * scenarioProbs[i];
    }
    // output[decisionId] = optionSwitch[decisionId % NUM_OPTIONS];
    // output[decisionId] = decisionId;
    // output[decisionId] = minCost[decisionId % 9];
    // output[decisionId] = scenarioProbs[decisionId % NUM_SCENARIOS];
    // output[decisionId] = tmpscenarioProbs[decisionId % NUM_SCENARIOS];
    // output[decisionId] = optionPointer[decisionId % NUM_OPTIONS];
    output[decisionIx] = expectedRisk;
  }
}