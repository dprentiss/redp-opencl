/* Implementation of Dial's Algorithm B 
 * algorithmB.h flattened with pointers
 * removed for use in OpenCL kernel
 */

#ifndef B_FLAT
#define B_FLAT

#include <math.h>

float setCost(int arc, float fftt[], float flow[], float cap[])
{
  return fftt[arc] * (1.0F + 0.15F*powf((flow[arc]/cap[arc]), 4.0F));
}

float calcCost(int arc,
               float deltaV,
               float fftt[],
               float flow[],
               float cap[])
{
  return fftt[arc] * (1.0F + 0.15F*powf(((flow[arc]+deltaV)/cap[arc]), 4.0F));
}

float calcCostPrime(int arc, float deltaV, float fftt[], float flow[], float cap[])
{
  return 0.6F * fftt[arc] * powf(flow[arc] + deltaV, 3.0F) / powf(cap[arc], 4.0F);
}

int isInNodeArray(int array[], int n)
{
  int i;
  for(i = 0; i < NUM_NODES; ++i) {
    if (array[i] == n) {
      return 1;
    }
  }
  return 0;
}

void updateCost(float minCost[],
                float maxCost[],
                float cost[],
                float flow[],
                int sinkNode,
                int pointer[],
                int tail[],
                int head[],
                int topoOrder[],
                int nextMinArc[],
                int nextMaxArc[])
{
  int i;
  int arc;
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

void getPaths(int origin,
              int sinkNode,
              int head[],
              int minPath[],
              int maxPath[],
              int nextMinArc[],
              int nextMaxArc[])
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
               int maxPath[],
               int minPath[],
               int pointer[],
               int tail[],
               int head[],
               int topoOrder[],
               int nextMinArc[],
               int nextMaxArc[],
               float flow[],
               float fftt[],
               float cap[],
               float cost[],
               float minCost[],
               float maxCost[])
{
  int i = 0;
  float nextFlow;
  float max = FLOAT_MAX;
  float root = 0;
  float f = 1.0F;
  float fPrime;
  int maxIter = 10;
  int iter = 0;

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
  root = max;

  while(fabs(f) > 0.0001F) {
    if(iter >= maxIter) {
      // printf("******************* Max Iterations\n");
      // printf("maximum iterations: %d %f %f %f %f\n", minPath[0], max, f, root, fPrime);
      // printStatus();
      break;
    }
    ++iter;
    i = 0;
    f = 0; 
    fPrime = 0; 
    while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
      if(minPath[i+1] < NUM_NODES) {
        f += calcCost(nextMinArc[minPath[i]], root, fftt, flow, cap);
        fPrime += calcCostPrime(nextMinArc[minPath[i]], root, fftt, flow, cap);
      }
      if(maxPath[i+1] < NUM_NODES) {
        f -= calcCost(nextMaxArc[maxPath[i]], -1*root, fftt, flow, cap);
        fPrime += calcCostPrime(nextMaxArc[maxPath[i]], -1*root, fftt, flow, cap);
      }
      ++i;
    }
    // printf("%d %f %f %f %f\n", minPath[0], max, f, root, fPrime);
    root -= f/fPrime;
  }
  if(root > max) {
    root = max;
  } else if(root < 0.0001f) {
    root = 0.1 * max;
    // return;
  } 
  i = 0;
  while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
    if(minPath[i+1] < NUM_NODES) {
      flow[nextMinArc[minPath[i]]] += root;
      cost[nextMinArc[minPath[i]]] = setCost(nextMinArc[minPath[i]], fftt, flow, cap);
    }
    if(maxPath[i+1] < NUM_NODES) {
      flow[nextMaxArc[maxPath[i]]] -= root;
      cost[nextMaxArc[maxPath[i]]] = setCost(nextMaxArc[maxPath[i]], fftt, flow, cap);
    }
    ++i;
  }
  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
  /*
  */
}

void
  equilibrateBush(
  int sinkNode, 
  int maxPath[],
  int minPath[],
  int nextMinArc[],
  int nextMaxArc[],
  int pointer[],
  int demand[],
  int tail[],
  int head[],
  int topoOrder[],
  float flow[],
  float fftt[],
  float cap[],
  float cost[],
  float minCost[],
  float maxCost[]
)
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
        shiftFlow(sinkNode, maxPath, minPath, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, flow, fftt, cap, cost, minCost, maxCost);
        updated = 1;
      }
    }
    if(updated == 0) {
      eqBush = 1;
    }
    iter++;
  }
}

void initNetwork(
  int sinkNode,
  int topoOrder[],
  int tail[],
  int head[],
  int pointer[],
  int demand[],
  int *nextMinArc,
  int nextMaxArc[],
  float minCost[],
  float maxCost[],
  float cost[],
  float fftt[],
  float cap[],
  float flow[]
)
{
  int i;
  int place;
  int queue;
  int arc;
  float costCand;
  int nextArc;

  // set topological order of nodes
  for(i = 0; i < NUM_NODES; ++i) {
    topoOrder[i] = 0;
  }
  topoOrder[0] = sinkNode;
  queue = 1;
  for(place = 0; place < NUM_NODES; ++place) {
    for(i=0; i < NUM_ARCS; ++i) {
      if (head[i] == topoOrder[place] && !isInNodeArray(topoOrder, tail[i])) {
        topoOrder[queue] = tail[i];
        ++queue;
      }
    }
  }

  // init cost
  for(i = 0; i < NUM_ARCS; ++i) {
    cost[i] = FLOAT_MAX;
  }

  // init minimum paths
  for(i = 0; i < NUM_NODES; i++) {
    minCost[i] = FLOAT_MAX;
  }
  minCost[sinkNode-1] = 0;
  nextMinArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    arc = pointer[topoOrder[i]-1]-1;
    while(head[arc] == topoOrder[i]) {
      costCand = fftt[arc] + minCost[head[arc]-1];
      if(costCand < minCost[tail[arc]-1]) {
        minCost[tail[arc]-1] = costCand;
        nextMinArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
  }
  // printStatus();

  // init maximum paths
  for(i = 0; i < NUM_NODES; i++) {
    maxCost[i] = minCost[i];
  }

  /*
  nextMaxArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
  arc = pointer[topoOrder[i]-1]-1;
  while(head[arc] == topoOrder[i]) {
  costCand = fftt[arc] + maxCost[head[arc]-1];
  if(costCand > maxCost[tail[arc]-1]) {
  maxCost[tail[arc]-1] = costCand;
  nextMaxArc[tail[arc]-1] = arc;
  // printStatus();
  }
  ++arc;
  }
  }
  */

  // init flows on bush
  nextArc = NUM_ARCS;
  for(i=0; i < NUM_NODES; ++i) {
    nextArc = nextMinArc[i];
    while(nextArc < NUM_ARCS) {
      flow[nextArc] += demand[i];
      nextArc = nextMinArc[head[nextArc]-1];
    }
  }

  // init costs on bush
  for(arc = 0; arc < NUM_ARCS; ++arc) {
    if(minCost[head[arc]-1] < minCost[tail[arc]-1]) {
      cost[arc] = setCost(arc, fftt, flow, cap);
    }
  }

  //printStatus();
  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
  //printStatus();
}

int updateBush(
  int sinkNode,
  int topoOrder[],
  int tail[],
  int head[],
  int pointer[],
  int nextMinArc[],
  int nextMaxArc[],
  float minCost[],
  float maxCost[],
  float cost[],
  float fftt[],
  float cap[],
  float flow[]
)

{
  int i;
  int j;
  int isUpdated = 0;

  for(i = 0; i < NUM_ARCS; ++i) {
    if(flow[i] == 0 && cost[i] < FLOAT_MAX) {
      j = pointer[tail[i]-1]-1;
      while(head[j] == tail[i]) {
        if(tail[j] == head[i] && fftt[j] + minCost[head[j]-1] < minCost[tail[j]-1]) {
          cost[i] = FLOAT_MAX;
          cost[j] = setCost(j, fftt, flow, cap); 
          isUpdated = 1;
          // printStatus();
        }
        ++j;
      }
    }
  }
  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
  return isUpdated;
}

void equilibrateNetwork(
  int sinkNode,
  int maxPath[],
  int minPath[],
  int nextMinArc[],
  int nextMaxArc[],
  int pointer[],
  int demand[],
  int tail[],
  int head[],
  int topoOrder[],
  float flow[],
  float fftt[],
  float cap[],
  float cost[],
  float minCost[],
  float maxCost[]
)
{
  int bushUpdated = 1;
  while(bushUpdated) {
    equilibrateBush(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
    // printStatus();
    bushUpdated = updateBush(sinkNode, topoOrder, tail, head, pointer, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow); 
  }
}

int verifyEq(
  float minCost[],
  float maxCost[],
  int demand[],
  int sinkNode)
{
  int i;
  float nodeBalance[NUM_NODES];
  for(i = 0; i < NUM_NODES; i++) {
    nodeBalance[i] = demand[i];
    if(abs(maxCost[i] - minCost[i]) > .01f) {
      return 1;
    }
  }
  for(i = 0; i < NUM_ARCS; i++) {
    nodeBalance[head[i]-1] += flow[i];
    nodeBalance[tail[i]-1] -= flow[i];
    }
  for(i = 0; i < NUM_NODES; i++) {
    if(abs(nodeBalance[i]) > .01f && i != sinkNode-1) {
      return 2;
    }
  }
  return 0;
}

#endif // B_FLAT