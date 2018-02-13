#ifndef ALGORITHM_B
#define ALGORITHM_B

#include <stdio.h>
#include <math.h>

#define NUM_NODES 5
#define NUM_ARCS  14
#define FLOAT_MAX 999999

// node constants
const int sinkNode = 5;
const int pointer[NUM_NODES] = {1, 3, 6, 9, 13};
const float demand[NUM_NODES] = {50, 50, 50, 50, 0};

// arc constants
const int tail[NUM_ARCS] = {2, 4, 1, 3, 4, 2, 4, 5, 1, 2, 3, 5, 3, 4};
const int head[NUM_ARCS] = {1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5};
const float fftt[NUM_ARCS] = {2, 11, 2, 7, 4, 7, 5, FLOAT_MAX, 11, 4, 5, FLOAT_MAX, 6, 4};
const float cap[NUM_ARCS] = {5, 30, 5, 10, 8, 10, 5, 1, 30, 20, 25, 1, 50, 50};

// node variables
int topoOrder[NUM_NODES];
float minCost[NUM_NODES];
float maxCost[NUM_NODES];

// arc variables
float flow[NUM_ARCS];
float cost[NUM_ARCS];
float costPrime[NUM_ARCS];
int nextMaxArc[NUM_NODES];
int nextMinArc[NUM_NODES];
int minPath[NUM_NODES];
int maxPath[NUM_NODES];

void setCost(int* arc) {
  cost[*arc] = fftt[*arc] * (1.0F + 0.15F*pow((flow[*arc]/cap[*arc]), 4.0F));
}

float calcCost(int* arc, float deltaV) {
  return fftt[*arc] * (1.0F + 0.15F*pow(((flow[*arc]+deltaV)/cap[*arc]), 4.0F));
}

float calcCostPrime(int* arc, float deltaV) {
  return 0.6F * fftt[*arc] * pow(flow[*arc] + deltaV, 3.0F) / pow(cap[*arc], 4.0F);
}

int isInNodeArray(int* a, int n) {
  int i;
  for(i = 0; i < NUM_NODES; ++i) {
    if (a[i] == n) {
      return 1;
    }
  }
  return 0;
}

void setTopoOrder(int sinkNode) {
  int place = 0;
  int queue = 0;
  int i = 0;
  topoOrder[queue] = sinkNode;
  ++queue;
  for(place; place < NUM_ARCS; ++place) {
    for(i=0; i < NUM_ARCS; ++i) {
      if (head[i] == topoOrder[place] && !isInNodeArray(topoOrder, tail[i])) {
        topoOrder[queue] = tail[i];
        ++queue;
      }
    }
  }
}

void setMinCost() {
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
      costCand = fftt[arc] + minCost[head[arc]-1];
      if(costCand < minCost[tail[arc]-1]) {
        minCost[tail[arc]-1] = costCand;
        nextMinArc[tail[arc]-1] = arc;
      }
      ++arc;
    }
  }
}

void setMaxCost() {
  int i;
  int arc;
  float costCand;

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
}

void updateMinCost() {
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
}

void updateMaxCost() {
  int i;
  int arc;
  float costCand;

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

void setArcCost() {
  int arc;
  for(arc = 0; arc < NUM_ARCS; ++arc) {
    if(minCost[head[arc]-1] < minCost[tail[arc]-1]) {
      setCost(&arc);
    }
  }
}

void initBush() {
  int i;
  int nextArc = NUM_ARCS;
  for(i=0; i < NUM_NODES; ++i) {
    nextArc = nextMinArc[i];
    while(nextArc < NUM_ARCS) {
      flow[nextArc] += demand[i];
      nextArc = nextMinArc[head[nextArc]-1];
    }
  }
}

void getPaths(int* origin) {
  int i, j, k;
  int length = 0;
  int hasCommonNode = 1;

  minPath[length] = *origin;
  maxPath[length] = *origin;

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
    }
    if(maxPath[length - 1] != sinkNode-1) {
      maxPath[length] = head[nextMaxArc[maxPath[length-1]]]-1;
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

void printStatus() {
  int i;
  printf("%3s%8s%8s%8s%11s%10s%11s%10s\n", "", "node", "pointer", "demand", "nextMinArc", "minCost", "nextMaxArc", "maxCost");
  for(i = 0; i < NUM_NODES; i++) {
    printf("%3d%8d%8d%8.2f%11d%10.2f%11d%10.2f\n", i, i+1, pointer[i], demand[i], nextMinArc[i], minCost[i], nextMaxArc[i], maxCost[i]);
  }
  printf("\n");

  printf("%3s%8s%12s%12s\n", "", "node", "minPath", "maxPath");
  for(i = 0; i < NUM_NODES; i++) {
    printf("%3d%8d%12d%12d\n", i, i+1, minPath[i], maxPath[i]);
  }
  printf("\n");

  printf("%3s%8s%8s%10s%10s\n", "", "head", "tail", "cost", "flow");
  for(i = 0; i < NUM_ARCS; i++) {
    printf("%3d%8d%8d%10.2f%10.2f\n", i, head[i], tail[i], cost[i], flow[i]);
  }
  printf("\n");
}

void shiftFlow() {
  int i = 0;
  float min = 0;
  float nextFlow;
  float max = FLOAT_MAX;
  float root = 0;
  float f = 1.0F;
  float fPrime;

  for(i = 0; i < NUM_NODES; ++i) {
    if(maxPath[i] >= NUM_NODES - 1) {
      break;
    }
    nextFlow = flow[nextMaxArc[maxPath[i]]];
    if(nextFlow < max) {
      max = nextFlow;
    }
  }
  root = max;

  while(fabs(f) > 0.001F) {
    i = 0;
    f = 0; 
    fPrime = 0; 
    while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
      if(minPath[i+1] < NUM_NODES) {
        f += calcCost(&nextMinArc[minPath[i]], root);
        fPrime += calcCostPrime(&nextMinArc[minPath[i]], root);
      }
      if(maxPath[i+1] < NUM_NODES) {
        f -= calcCost(&nextMaxArc[maxPath[i]], -1*root);
        fPrime += calcCostPrime(&nextMaxArc[maxPath[i]], -1*root);
      }
      ++i;
    }
    //printf("%d %f %f %f %f\n", minPath[0], max, f, root, fPrime);
    root -= f/fPrime;
  }
  if(root > max) {
    root = max;
  //} else if(root < min) {
  //  root = min;
  } 
  i = 0;
  while(minPath[i+1] < NUM_NODES || maxPath[i+1] < NUM_NODES) {
      if(minPath[i+1] < NUM_NODES) {
        flow[nextMinArc[minPath[i]]] += root;
        setCost(&nextMinArc[minPath[i]]);
      }
      if(maxPath[i+1] < NUM_NODES) {
        flow[nextMaxArc[maxPath[i]]] -= root;
        setCost(&nextMaxArc[maxPath[i]]);
      }
      ++i;
  }
  updateMinCost();
  updateMaxCost();
}

void equilibrateBush() {
  int i;
  int eqBush = 0;
  float nextCost = 0;

  while(!eqBush) {
    for(i = 0; i < NUM_NODES; ++i) {
      if(maxCost[i] - minCost[i] > 0.01F && demand[i] > 0) {
        // printStatus();
        getPaths(&i);
        nextCost = maxCost[i];
        shiftFlow();
        updateMaxCost();
        updateMinCost();
        if(maxCost[i] == nextCost) {
          continue;
        }
        // printStatus();
        break;
      }
      if(i == NUM_NODES - 1) {
        eqBush = 1;
      }
    }
  }
}

void initNetwork() {
  int i;

  for(i = 0; i < NUM_ARCS; ++i) {
    cost[i] = FLOAT_MAX;
  }

  setTopoOrder(sinkNode);
  setMinCost();
  setMaxCost();
  initBush();
  setArcCost();
  updateMinCost();
  updateMaxCost();
}

int updateBush() {
  int i;
  int j;
  int isUpdated = 0;

  for(i = 0; i < NUM_ARCS; ++i) {
    if(flow[i] == 0 && cost[i] < FLOAT_MAX) {
      j = pointer[tail[i]-1]-1;
      while(head[j] == tail[i]) {
        if(tail[j] == head[i]
        && fftt[j] + minCost[tail[j]-1] > minCost[head[j]-1]) {
          cost[i] = FLOAT_MAX;
          setCost(&j); 
          isUpdated = 1;
          break;
        }
      ++j;
      }
    }
  }
  updateMaxCost();
  updateMinCost();
  return isUpdated;
}

void equilibrateNetwork() {
  int bushUpdated = 1;
  while(bushUpdated) {
    equilibrateBush();
    bushUpdated = updateBush(); 
  }
}

#endif // ALGORITHM_B