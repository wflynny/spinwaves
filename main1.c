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
//       float s[3];   //this notation is cuasing problems when i try to assign it a value for some reason
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

//test
/*struct atomArray
{
       struct atom atomArr[];
}

*/

/*
typedef struct {
  double *Spp, *Smm, *Spm, *Smp ;
  double *Epp, *Emm, *Epm, *Emp ;
  double *R ;
} PBoutdata ;
*/

/*
typedef struct
{
       int num;
       int ints[];
} testStruct;

int basicTest(testStruct  *struc)
{
    int i;
    printf("%d\n", struc->num);
    for( i = 0; i < struc->num; i++)
    { 
        printf("i=%d struc.ints=%d\n",i, struc->ints[i]);
    }
    return struc->num + 1;
}

typedef struct
{
        int num;
        testStruct *structs[];
}ArrayStruct;

void test2(ArrayStruct *struc)
{
     int i;
     int j;
     for( j = 0; j < struc->num; j++)
     {
              for( i = 0; i < struc->structs[j]->num; i++)
              { 
                printf("i=%d struc.ints=%d\n",i, struc->structs[j]->ints[i]);
              }
     }
}

int test(struct atom atoms[], int numAtoms)
{
//    struct atom *atoms = atomArray.atomArr;
    int i = 0;
    int j;
    struct atom a;
    for(i = 0; i < numAtoms; i++)
    {
      a = atoms[i];
      printf("%d) NumInteractions: %d\n",i,a.numInteractions);
      printf("%d) s[0]: %d\n",i,a.s[0]);
      printf("%d) s[1]: %d\n",i,a.s[1]);
      printf("%d) s[2]: %d\n",i,a.s[2]);
      printf("%d) pos[0]: %d\n",i,a.pos[0]);
      printf("%d) pos[1]: %d\n",i,a.pos[1]);
      printf("%d) pos[2]: %d\n",i,a.pos[2]);
      for(j = 0; j < a.numInteractions; j++)
      {
            printf("%d,%d) NumInteractions: %d\n",i,j, a.interactions[j].jMatrixIndex);
            printf("%d,%d) NumInteractions: %d\n",i,j, a.interactions[j].otherAtomIndex);
      }
    }
    
    return 0;
}

*/
/*
typedef struct
{
        int length;
        int interactions[];
}interactionsStruct;

typedef struct
{
        int number;
        float sArray[][3][3];
        interactionsStruct interactionsArray[];
}atomArrayStruct;  

void test3(int k , float tMax, float tMin, float tFactor, atomArrayStruct atoms)
{
     return
}
*/

typedef float Spin[3];

//Global list of atoms so that the atoms can be added one at a time
//Atom *atoms;  //not necessary if python keeps the address


typedef struct
{
        int *interaction_matrix;
        int *nbr_list;
        float *spin;
//        Spin spin;
} Atom;
        
Atom* new_atom_list(int n)
{
     Atom *a;
     a=malloc(n*sizeof(Atom));
     return a;//return to calling python method to add the atoms
}

void del_atom(Atom *p) {free(p);}

void set_atom(Atom *p, int k, int *mat, int matLength, int* neighbors, int nbrLength, Spin spin)
{
     p[k].interaction_matrix = mat;
     p[k].nbr_list = neighbors;
     p[k].spin = spin;
}

void atomTest(Atom *p, int num)
{
     int length;
     length = sizeof(*p)/sizeof(Atom);
     printf("length: %d\n", length);  //just a test, does not know the amount allocated to it with malloc
     
     int i;
     for (i = 0; i < num; i++)
     {
         printf("atom %d) mat[0]:%d mat[1]:%d mat[2]:%d\n", i, p[i].interaction_matrix[0], p[i].interaction_matrix[1], p[i].interaction_matrix[2]);
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
