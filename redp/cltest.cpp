#include <stdio.h>
#include <CL/cl.h>
 
int main(int argc, char **argv){
    cl_platform_id test;
    cl_uint num;
    cl_uint ok = 1;
    clGetPlatformIDs(ok, &test, &num);

    printf("%d, %d\n", test, num);
 
    return 0;
}