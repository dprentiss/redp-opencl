#include "math.h"
#include "float.h"

#ifndef NUM_NODES
#define NUM_NODES 14
#endif

#ifndef NUM_ARCS
#define NUM_ARCS 53
#endif

#ifndef FLOAT_MAX
#define FLOAT_MAX CL_INFINITY
#endif

#ifndef NUM_OPTIONS
#define NUM_OPTIONS 6
#endif

#ifndef NUM_SCENARIOS
#define NUM_SCENARIOS 4
#endif


// node variables
int topoOrder[NUM_NODES];
float minCost[NUM_NODES];
float maxCost[NUM_NODES];
int nextMaxArc[NUM_NODES];
int nextMinArc[NUM_NODES];
int minPath[NUM_NODES];
int maxPath[NUM_NODES];

// arc variables
float flow[NUM_ARCS];
float cost[NUM_ARCS];
float arcRiskA[NUM_ARCS];
float arcRiskB[NUM_ARCS];

float maxRisk[NUM_SCENARIOS];
bool optionSwitch[NUM_OPTIONS];

void printStatus()
{
  int i;
  /*
  */
  printf("%3s%8s%8s%8s%11s%10s%11s%10s\n", "", "node", "pointer", "demand", "nextMinArc", "minCost", "nextMaxArc", "maxCost");
  for(i = 0; i < NUM_NODES; i++) {
    printf("%3d%8d%8d%8d%11d%10.2f%11d%10.2f\n", i, i+1, pointer[i], demand[i], nextMinArc[i], minCost[i], nextMaxArc[i], maxCost[i]);
  }
  printf("\n");

  printf("%3s%8s%12s%12s%12s\n", "", "node", "minPath", "maxPath", "topoOrder");
  for(i = 0; i < NUM_NODES; i++) {
    printf("%3d%8d%12d%12d%12d\n", i, i+1, minPath[i], maxPath[i], topoOrder[i]);
  }
  printf("\n");

  printf("%3s%8s%8s%10s%10s\n", "", "head", "tail", "cost", "flow");
  for(i = 0; i < NUM_ARCS; i++) {
    printf("%3d%8d%8d%10.2f%10.2f\n", i, head[i], tail[i], cost[i], flow[i]);
  }
  printf("\n");
}

float getMaxCap(int arc, cl_uint *options, cl_uint *optionPointer) {
  int i;
  for(i = 0; i < NUM_OPTIONS; ++i) {
    if(arc == options[optionPointer[i]]) {
      return maxCaps[i];
    }
  }
  return 0;
}

float setCost(int arc, const float *fftt, float *flow, const float *cap, float *arcRiskA, float *arcRiskB)
{
  return arcRiskA[arc] * fftt[arc] 
  * (1.0F + 0.15F*powf((flow[arc]/cap[arc]), 4.0F)) + arcRiskB[arc];
}

float calcCost(int arc,
               float deltaV,
               const float *fftt,
               float *flow,
               const float *cap, 
               float *arcRiskA,
               float *arcRiskB)
{
  return arcRiskA[arc] * fftt[arc] 
  * (1.0F + 0.15F*powf(((flow[arc]+deltaV)/cap[arc]), 4.0F)) 
    + arcRiskB[arc];
}

float calcCostPrime(int arc, float deltaV, const float *fftt, float *flow, const float *cap, float *arcRiskA, float *arcRiskB)
{
  return 0.6F * arcRiskA[arc] * fftt[arc] * powf(flow[arc] + deltaV, 3.0F) / powf(cap[arc], 4.0F);
}

int isInNodeArray(int *array, int n)
{
  int i;
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
                int sinkNode,
                cl_uint *pointer,
                cl_uint *tail,
                cl_uint *head,
                int *topoOrder,
                int *nextMinArc,
                int *nextMaxArc)
{
  int i;
  int arc;
  int hasCap;
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
        if(hasCap = getMaxCap(arc, options, optionPointer)) {
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
    //nextMaxArc[i] = 0;
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

void getPaths(int origin,
              int sinkNode,
              cl_uint *head,
              int *minPath,
              int *maxPath,
              int *nextMinArc,
              int *nextMaxArc)
{
  int i, j, k;
  int length = 0;
  int hasCommonNode = 1;

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

void shiftFlow(int sinkNode,
               int *maxPath,
               int *minPath,
               cl_uint *pointer,
               cl_uint *tail,
               cl_uint *head,
               int *topoOrder,
               int *nextMinArc,
               int *nextMaxArc,
               float *flow,
               const float *fftt,
               const float *cap,
               float *cost,
               float *minCost,
               float *maxCost)
{
  int i = 0;
  float nextFlow;
  float max = FLOAT_MAX;
  float shelterCap;
  float root = 0;
  float f = 1.0F;
  float fPrime;
  int maxIter = 10;
  int iter = 0;

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
    if(shelterCap = getMaxCap(nextMinArc[minPath[i]], options, optionPointer)) {
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
    root = 0.1f * max;
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
  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
}

void equilibrateBush(int sinkNode, 
                     int *maxPath,
                     int *minPath,
                     int *nextMinArc,
                     int *nextMaxArc,
                     cl_uint *pointer,
                     cl_uint *demand,
                     cl_uint *tail,
                     cl_uint *head,
                     int *topoOrder,
                     float *flow,
                     const float *fftt,
                     const float *cap,
                     float *cost,
                     float *minCost,
                     float *maxCost)
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
        // printStatus();
        shiftFlow(sinkNode, maxPath, minPath, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, flow, fftt, cap, cost, minCost, maxCost);
        // printStatus();
        updated = 1;
      }
    }
    if(updated == 0) {
      eqBush = 1;
    }
    iter++;
  }
}

void initNetwork(const int sinkNode,
                 int *topoOrder,
                 cl_uint *tail,
                 cl_uint *head,
                 cl_uint *pointer,
                 cl_uint *demand,
                 int *nextMinArc,
                 int *nextMaxArc,
                 float *minCost,
                 float *maxCost,
                 float *cost,
                 const float *fftt,
                 const float *cap,
                 float *flow)
{
  int i;
  int lastI;
  int j;
  int k;
  int l;
  int tmpNode;
  int queue;
  int arc;
  float costCand;
  int nextArc;
  float max;

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
      if(costCand <= minCost[tail[arc]-1] && !getMaxCap(arc, options, optionPointer)) {
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
  //printStatus();

  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
  // printStatus();
}

int updateBush(int sinkNode,
               int topoOrder[],
               cl_uint *tail,
               cl_uint *head,
               cl_uint *pointer,
               int nextMinArc[],
               int nextMaxArc[],
               float minCost[],
               float maxCost[],
               float cost[],
               const float *fftt,
               const float *cap,
               float flow[])

{
  int i;
  int j;
  int arc;
  int isUpdated = 0;
  int queue;
  int place;

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

  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
  return isUpdated;
}

int eqStatus(
  float minCost[],
  float maxCost[],
  cl_uint *demand,
  int sinkNode)
{
  int i;
  float nodeBalance[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++) {
    nodeBalance[i] = demand[i];
    if(abs(maxCost[i] - minCost[i]) > .01f) {
      printf("costDiff, node %d: %f", i, abs(maxCost[i] - minCost[i]));
      return 1;
    }
  }
  for(i = 0; i < NUM_ARCS; i++) {
    nodeBalance[head[i]-1] += flow[i];
    nodeBalance[tail[i]-1] -= flow[i];
  }
  for(i = 0; i < NUM_NODES; i++) {
    if(abs(nodeBalance[i]) > .01f && i != sinkNode-1) {
      printf("nodeBalance, node %d: %f", i, abs(maxCost[i] - minCost[i]));
      return 2;
    }
  }
  return 0;
}

void equilibrateNetwork(int sinkNode,
                        int *maxPath,
                        int *minPath,
                        int *nextMinArc,
                        int *nextMaxArc,
                        cl_uint *pointer,
                        cl_uint *demand,
                        cl_uint *tail,
                        cl_uint *head,
                        int *topoOrder,
                        float *flow,
                        const float *fftt,
                        const float *cap,
                        float *cost,
                        float *minCost,
                        float *maxCost)
{
  int bushUpdated = 1;
  while(bushUpdated /*|| eqStatus(minCost, maxCost, demand, sinkNode)*/) {
    equilibrateBush(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
    bushUpdated = updateBush(sinkNode, topoOrder, tail, head, pointer, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow); 
  }
}



int verifyEq(
  float minCost[],
  float maxCost[],
  cl_uint *demand,
  int sinkNode)
{
  int i;
  float nodeBalance[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++) {
    nodeBalance[i] = demand[i];
    if(abs(maxCost[i] - minCost[i]) > .01f) {
      std::cout << abs(maxCost[i] - minCost[i]) << " ";
      return 1;
    }
  }
  for(i = 0; i < NUM_ARCS; i++) {
    nodeBalance[head[i]-1] += flow[i];
    nodeBalance[tail[i]-1] -= flow[i];
  }
  for(i = 0; i < NUM_NODES; i++) {
    if(abs(nodeBalance[i]) > .01f && i != sinkNode-1) {
      std::cout << abs(nodeBalance[i]) << " ";
      return 2;
    }
  }
  return 0;
}

