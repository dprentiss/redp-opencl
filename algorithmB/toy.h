#define NUM_NODES 5
#define NUM_ARCS  14
#define FLOAT_MAX 999999

// node constants
const int sinkNode = 5;
const int pointer[NUM_NODES] = {1, 3, 6, 9, 13};
const float demand[NUM_NODES] = {0, 50, 50, 50, 0};

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