#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
//#include "C:/random/dSFMT-src-1.3/dSFMT.h"

float randf()
{
    return rand()/(((float)RAND_MAX)+1);  //Since RAND_MAX is only 32000 this only gives that many unique values
}

float randFloat(float min, float max)  //uses randfloat() to generate a random float between min and max
{
    if (min>max)
    {
        return randf()*(min-max)+max;
    }
    else
    {
        return randf()*(max-min)+min;
    }
}

struct interaction
{
       int otherAtomIndex;  //may need to be long if the list gets big
       int jMatrixIndex;
};

struct atom
{
//       float s[3];   //this notation is cuasing problems when i try to assign it a value for soe reason
//       float pos[3];
       int numInteractions;
//       struct interaction interactions[];
         float *s;
         float *pos;
         struct interaction* interactions;
};
typedef float InteractionMatrix[3][3];

float matMultiply(float s1[3], InteractionMatrix jMat, float s2[3])
{
      //multiply s1 by jMat
      float x1 = (jMat[0][0] * s1[0]) + (jMat[0][1] * s1[1]) + (jMat[0][2] * s1[2]);
      float y1 = (jMat[1][0] * s1[0]) + (jMat[1][1] * s1[1]) + (jMat[1][2] * s1[2]);
      float z1 = (jMat[2][0] * s1[0]) + (jMat[2][1] * s1[1]) + (jMat[2][2] * s1[2]);

      //return the dot product
      return ((x1*s2[0]) + (y1*s2[1]) + (z1*s1[2]));
}


/* I'm hoping that python will do this malloc business for me
void caller()
{
  float *d;
  InteractionMatrix *jMatrices;
  d = (float*)malloc(sizeof(float)*3*3*15);
  jMatrices = (InteractionMatrix*)(d);
  Energy(atoms, a, jMatrices)
}
*/

float Energy(struct atom atoms[], struct atom a, float s1[3], InteractionMatrix* jMatrices)
{
      float E = 0;
      float *s2;
      int j;
      int index;
      for(j = 0; j < a.numInteractions; j++)
      {
              s2 = atoms[a.interactions[j].otherAtomIndex].s;
              index = a.interactions[j].jMatrixIndex;
              E -= matMultiply(s1, jMatrices[index], s2);
      }
      return E;
}

void flipSpins(struct atom atoms[], int numAtoms, float jMatrices[][3][3], float T)
{
     int index;
     float oldE;
     float newE;
     for(index = 0; index < numAtoms; index++)
     {
             float newS[3] = {0,0,randFloat(-1, 1)};
             oldE = Energy(atoms, atoms[index], atoms[index].s, jMatrices);
             newE = Energy(atoms, atoms[index], newS, jMatrices);
             if(newE < oldE)
             {
                     atoms[index].s = newS;
             }
             else
             {
                 float deltE = newE - oldE;
                 float probChange = exp(-deltE/T);
                 if(randf() < probChange)
                 {
                            atoms[index].s = newS;
                 }
             }
     }
}


//numAtoms may need to be switched to long for long lists of atoms
void simulate(struct atom atoms[], int numAtoms, float jMatrices[][3][3], int k, float maxTemp, float minTemp, float tFactor)
{
     //set the seed for the random numbers
     time_t *time1;
     
     srand(time(time1));

     int i;
     float T = maxTemp;
     while(T < minTemp)
     {
             for(i = 0; i < k; i++)
             {
                     flipSpins(atoms, numAtoms, jMatrices, T);
             }
             T = T*tFactor;
     }
}


int main(int argc, char *argv[])
{
    int k = 1000;
    float tMax = 15;
    float tMin = .01;
    float tFactor = .9;
    /*time(NULL)
    long seed = 1;
    srand48(seed); */
    int random = rand();


    printf("%d",random);
    system("PAUSE");
    return 0;
}
