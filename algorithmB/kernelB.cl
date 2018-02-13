#define NUM_NODES 25
#define NUM_ARCS  80
#define FLOAT_MAX 999999

float setCost(int arc, float fftt[], float flow[], float cap[])
{
  return fftt[arc] * (1.0F + 0.15F*powr((flow[arc]/cap[arc]), 4.0F));
}

float calcCost(int arc,
               float deltaV,
               float fftt[],
               float flow[],
               float cap[])
{
  return fftt[arc] * (1.0F + 0.15F*powr(((flow[arc]+deltaV)/cap[arc]), 4.0F));
}

float calcCostPrime(int arc, float deltaV, float fftt[], float flow[], float cap[])
{
  return 0.6F * fftt[arc] * powr(flow[arc] + deltaV, 3.0F) / powr(cap[arc], 4.0F);
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
    maxCost[i] = 0;
  }

  maxCost[sinkNode-1] = 0;
  nextMaxArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    arc = pointer[topoOrder[i]-1]-1;
    while(head[arc] == topoOrder[i]) {
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

  minPath[length] = origin;
  maxPath[length] = origin;

  while(hasCommonNode) {
    if(nextMaxArc[minPath[length]] == nextMinArc[maxPath[length]]) {
      minPath[length] = head[nextMaxArc[minPath[length]]]-1; 
      maxPath[length] = minPath[length]; 
      continue;
    }
    hasCommonNode = 0;
  }

  while(!hasCommonNode) {
    ++length;
    if(minPath[length - 1] != sinkNode-1) {
      minPath[length] = head[nextMinArc[minPath[length-1]]]-1;
    } else {
      minPath[length] = NUM_NODES;
    }
    if(maxPath[length - 1] != sinkNode-1) {
      maxPath[length] = head[nextMaxArc[maxPath[length-1]]]-1;
    } else {
      maxPath[length] = NUM_NODES;
    }
    for(i = 1; i < length + 1; ++i) {
      for(j = 1; j < length + 1; ++j) {
        if(maxPath[i] == minPath[j]) {
          for(k = j; k < NUM_NODES; ++k) {
            minPath[k + 1] = NUM_NODES;
          }
          for(k = i; k < NUM_NODES; ++k) {
            maxPath[k + 1] = NUM_NODES;
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
    if(maxPath[i] >= NUM_NODES - 1) {
      break;
    }
    nextFlow = flow[nextMaxArc[maxPath[i]]];
    if(nextFlow < max) {
      max = nextFlow;
    }
  }
  root = 0;

  while(fabs(f) > 0.001F) {
    if(iter >= maxIter) {
      // printf("******************* Max Iterations\n");
      // printf("%d %f %f %f %f\n", minPath[0], max, f, root, fPrime);
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
    root = 0;
    return;
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
}

void equilibrateBush(int sinkNode, 
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
                     float maxCost[])
{
  int i;
  int eqBush = 0;
  float nextCost = 0;

  while(!eqBush) {
    for(i = 0; i < NUM_NODES; ++i) {
      if(maxCost[i] - minCost[i] > 0.01F && demand[i] > 0) {
        // printStatus();
        getPaths(i, sinkNode, head, minPath, maxPath, nextMinArc, nextMaxArc);
        nextCost = maxCost[i];
        shiftFlow(sinkNode, maxPath, minPath, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc, flow, fftt, cap, cost, minCost, maxCost);
        //updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
        if(maxCost[i] == nextCost) {
          continue;
        }
        //printStatus();
        break;
      }
      if(i == NUM_NODES - 1) {
        eqBush = 1;
      }
    }
  }
}

void initNetwork(int sinkNode,
                 int topoOrder[],
                 int tail[],
                 int head[],
                 int pointer[],
                 int demand[],
                 int nextMinArc[],
                 int nextMaxArc[],
                 float minCost[],
                 float maxCost[],
                 float cost[],
                 float fftt[],
                 float cap[],
                 float flow[])
{
  int i;
  int place = 0;
  int queue = 0;
  int arc;
  float costCand;
  int nextArc;

  // zero flow and cost arrays
  for(i = 0; i < NUM_ARCS; ++i) {
    flow[i] = 0;
    cost[i] = 0;
  }

  // set topological order of nodes
  i = 0;
  topoOrder[queue] = sinkNode;
  ++queue;
  for( ; place < NUM_ARCS; ++place) {
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

  // init maximum paths
  for(i = 0; i < NUM_NODES; i++) {
    maxCost[i] = 0;
  }
  maxCost[sinkNode-1] = 0;
  nextMaxArc[sinkNode-1] = NUM_ARCS;
  for(i = 0; i < NUM_NODES; i++) {
    arc = pointer[topoOrder[i]-1]-1;
    while(head[arc] == topoOrder[i]) {
      costCand = fftt[arc] + maxCost[head[arc]-1];
      if(costCand > maxCost[tail[arc]-1]) {
        maxCost[tail[arc]-1] = costCand;
        nextMaxArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
  }

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

  updateCost(minCost, maxCost, cost, flow, sinkNode, pointer, tail, head, topoOrder, nextMinArc, nextMaxArc);
}

int updateBush(int sinkNode,
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
               float flow[])

{
  int i;
  int j;
  int isUpdated = 0;

  for(i = 0; i < NUM_ARCS; ++i) {
    if(flow[i] == 0 && cost[i] < FLOAT_MAX) {
      j = pointer[tail[i]-1]-1;
      while(head[j] == tail[i]) {
        if(tail[j] == head[i]
        && fftt[j] + minCost[head[j]-1] < minCost[tail[j]-1]) {
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

void testFunc(float flow[])
{
  flow[2] = 666.f;
}

void equilibrateNetwork(int sinkNode,
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
                        float maxCost[])
{
  int bushUpdated = 1;
  while(bushUpdated) {
    equilibrateBush(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
    // printStatus();
    bushUpdated = updateBush(sinkNode, topoOrder, tail, head, pointer, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow); 
  }
}

/* DO NOT BULD. THIS DECLARATION IS FOR CHECKING SYNTAX ONLY */
/*
void kernelB(int * output,
             unsigned int * input,
             unsigned int multiplier)
             */

__kernel void kernelB(__global  unsigned int * output,
                      __global  unsigned int * input,
                      const     unsigned int multiplier)

{
  uint tid = get_global_id(0);

  int sinkNode = 24;
  int pointer[] = {1, 3, 6, 9, 12, 14, 17, 21, 25, 29, 32, 35, 39, 43, 47, 50, 53, 57, 61, 65, 68, 70, 73, 76, 79};

  int demand[] = {50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 0, 50};

  // arc constants
  int tail[] = {2,6,1,3,7,2,4,8,3,5,9,4,10,1,7,11,2,6,8,12,3,7,9,13,4,8,10,14,5,9,15,6,12,16,7,11,13,17,8,12,14,18,9,13,15,19,10,14,20,11,17,21,12,16,18,22,13,17,19,23,14,18,20,24,15,19,25,16,22,17,21,23,18,22,24,19,23,25,20,24};
  int head[] = {1,1,2,2,2,3,3,3,4,4,4,5,5,6,6,6,7,7,7,7,8,8,8,8,9,9,9,9,10,10,10,11,11,11,12,12,12,12,13,13,13,13,14,14,14,14,15,15,15,16,16,16,17,17,17,17,18,18,18,18,19,19,19,19,20,20,20,21,21,22,22,22,23,23,23,24,24,24,25,25};
  float fftt[] = {3.72f,3.72f,3.72f,5.4f,4.02f,5.4f,3.3f,1.37f,3.3f,4.98f,3.09f,4.98f,4.59f,3.72f,4.17f,4.32f,4.02f,4.17f,4.02f,4.11f,1.37f,4.02f,5.55f,1.41f,3.09f,5.55f,3.42f,4.11f,4.59f,3.42f,3.42f,4.32f,1.23f,4.35f,4.11f,1.23f,1.85f,5.49f,1.41f,1.85f,1.89f,2.15f,4.11f,1.89f,1.03f,5.46f,3.42f,1.03f,6.09f,4.35f,3.09f,5.58f,5.49f,3.09f,4.5f,4.02f,2.15f,4.5f,5.94f,1.79f,5.46f,5.94f,2.97f,5.85f,6.09f,2.97f,4.05f,5.58f,5.76f,4.02f,5.76f,4.44f,1.79f,4.44f,5.1f,5.85f,5.1f,3.f,4.05f,3.f};
  float cap[] = {300,300,300,300,300,300,300,200,300,300,300,300,300,300,300,300,300,300,300,300,200,300,300,200,300,300,300,300,300,300,300,300,200,300,300,200,200,300,200,200,200,200,300,200,200,300,300,200,300,300,300,300,300,300,300,300,200,300,300,200,300,300,300,300,300,300,300,300,300,300,300,300,200,300,300,300,300,300,300,300};

  // node variables
  int topoOrder[NUM_NODES];
  float minCost[NUM_NODES];
  float maxCost[NUM_NODES];

  // arc variables
  float flow[NUM_ARCS];
  float cost[NUM_ARCS];
  int nextMaxArc[NUM_NODES];
  int nextMinArc[NUM_NODES];
  int minPath[NUM_NODES];
  int maxPath[NUM_NODES];
  
  initNetwork(sinkNode, topoOrder, tail, head, pointer, demand, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow);
  // equilibrateNetwork(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
  output[tid] = flow[tid];
}