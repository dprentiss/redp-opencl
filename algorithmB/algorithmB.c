#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include "dial.h"


void printStatus()
{
  int i;
  printf("%3s%8s%8s%8s%11s%10s%11s%10s\n", "", "node", "pointer", "demand", "nextMinArc", "minCost", "nextMaxArc", "maxCost");
  for(i = 0; i < NUM_NODES; i++) {
    printf("%3d%8d%8d%8d%11d%10.2f%11d%10.2f\n", i, i+1, pointer[i], demand[i], nextMinArc[i], minCost[i], nextMaxArc[i], maxCost[i]);
  }
  printf("\n");

  /*
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
  */
}

#include "bflat.h"

int main(int argc, const char* argv[])
{
  int i;
  int j;
  // file input
  // usage is redpcl kernelFile networkFile
  if ( argc != 2 ) /* argc should be 2 for correct execution */ {
    printf( "\nusage: redcl kernelFile networkFile\n");
  }
  else {
    FILE *file = fopen( argv[1], "r" );
    if ( file == 0 ) {
      printf( "Could not open file\n" );
    }
    else {
      int x;
      while  ( ( x = fgetc( file ) ) != EOF ) {
        printf( "%c", x );
      }
      fclose( file );
    }
  }
  for(i = 0; i < 80; ++i) {
    for(j = 0; j < 80; ++j) {
      flow[j] = 0;
      cost[j] = 0;
    }
    sinkNode = (i % 25) + 1;
    initNetwork(sinkNode, topoOrder, tail, head, pointer, demand, nextMinArc, nextMaxArc, minCost, maxCost, cost, fftt, cap, flow);
    equilibrateNetwork(sinkNode, maxPath, minPath, nextMinArc, nextMaxArc, pointer, demand, tail, head, topoOrder, flow, fftt, cap, cost, minCost, maxCost);
    printf("sinkNode: %d, verified: %d, minCost[0]: %d\n", sinkNode, verifyEq(minCost, maxCost, demand, sinkNode), (int) minCost[0]);
    // printStatus();
  }
}