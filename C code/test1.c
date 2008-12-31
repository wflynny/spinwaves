#include <stdlib.h>
#include <stdio.h>

typedef float InteractionMatrix[3][3];

float energy(InteractionMatrix m)
{
  int j,k;
  printf("m=%p\n",m);
  for (j=0; j < 3; j++)
    for (k=0; k<3; k++)
        printf("m[%d,%d] = %f\n",j,k,m[j][k]);
  return 0;
}

void caller()
{
  int i;
  float *d;
  InteractionMatrix *jMatrices;
  d = (float*)malloc(sizeof(float)*3*3*15);
  for (i = 0; i < 15*9; i++) d[i] = i;
  jMatrices = (InteractionMatrix *)(d);
  printf("d=%p,jMatrices=%p\n",d,jMatrices);
  energy(jMatrices[3]);
}


int main(int argc, char *argv[])
{
   caller();
   return 0;
}
