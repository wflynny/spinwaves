#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "dSFMT.h"
/*
float randf()
{
    return rand()/(((float)RAND_MAX)+1);  //Since RAND_MAX is only 32000 this only gives that many unique values
}
*/
typedef float Spin[3];
typedef float InteractionMatrix[3][3];

static dsfmt_t *dsfmt_state;

float randf()
{
      return (float)dsfmt_genrand_close_open(dsfmt_state);
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

void initializeRandState()
{
     dsfmt_state = (dsfmt_t*)malloc(sizeof(dsfmt_t));
     dsfmt_init_gen_rand(dsfmt_state, time(NULL));
}

void randSpin(Spin spin)
{
       //Marsaglia Method
       float Sx, Sy, Sz;
       float r1, r2;
       float C1, C2, Csqr;
       Csqr = 1;
       while( Csqr >= 1)
       {
           r1 = randf();
           r2 = randf();
           C1 = 1-(2*r1);
           C2 = 1-(2*r2);
           Csqr = pow(C1, 2) + pow(C2, 2);
           if(Csqr < 1)
           {
                   Sx = 2 * C1 * sqrt(1-Csqr);
                   Sy = 2 * C2 * sqrt(1-Csqr);
                   Sz = 1 - 2*Csqr;
           }
       }
            
       spin[0] = Sx;
       spin[1] = Sy;
       spin[2] = Sz;
//       printf("Spin: (%f, %f, %f) ", Sx, Sy, Sz);
       return;
}

typedef struct
{
        int numInteractions;
        int *interaction_matrix;
        int *nbr_list;
//        float *spin;
        Spin spin;
} Atom;
        

InteractionMatrix* new_jMatrix_list(int n)
{
     InteractionMatrix *m;
     m = (InteractionMatrix*)malloc(n*sizeof(InteractionMatrix));
     return m;//return pointer to python
}

void add_matrix(InteractionMatrix *m, int index, float j11, float j12, float j13, float j21, float j22, float j23, float j31, float j32, float j33)
{
     m[index][0][0] = j11;
     m[index][0][1] = j12;
     m[index][0][2] = j13;
     m[index][1][0] = j21;
     m[index][1][1] = j22;
     m[index][1][2] = j23;
     m[index][2][0] = j31;
     m[index][2][1] = j32;
     m[index][2][2] = j33;
//     printf("%f %f %f %f %f %f %f %f %f\n", j11, j12, j13, j21, j22, j23, j31, j32, j33);
}     

void del_jMat(InteractionMatrix *m){free(m);}

Atom* new_atom_list(int n)
{
     Atom *a;
     a=(Atom*)malloc(n*sizeof(Atom));
     return a;//return to calling python method to add the atoms
}

void del_atom(Atom *p) {free(p);}

void set_atom(Atom *p, int k, int *mat, int* neighbors, int numInteractions, Spin spin)
{
     p[k].interaction_matrix = mat;
     p[k].nbr_list = neighbors;
     
     p[k].spin[0] = spin[0];
     p[k].spin[1] = spin[1];
     p[k].spin[2] = spin[2];
     p[k].numInteractions = numInteractions;
     
     //test
//     printf("mat* = %d   mat[0],[1],[2] = (%d,%d,%d)   p[k].mat* = %d p[k].mat[0],[1],[2] = (%d,%d, %d)\n", mat, mat[0], mat[1], mat[2], p[k].interaction_matrix, p[k].interaction_matrix[0],p[k].interaction_matrix[1], p[k].interaction_matrix[2]);
}

void atomTest(Atom *p, int num)
{
     int i;
     int j;
     for(i = 0; i < num; i++)
     {
           printf("%d) Spin: (%d,%d,%d) interactions: %d( ", i, p[i].spin[0], p[i].spin[1], p[i].spin[2], p[i].nbr_list);
           for(j = 0; j < p[i].numInteractions; j++)
           {
                 printf("%d ", p[i].nbr_list[j]);
           }
           printf(") jMats: %d( ", p[i].interaction_matrix);
           for(j = 0; j < p[i].numInteractions; j++)
           {
                 printf("%d ", p[i].interaction_matrix[j]);
           }
           printf(")\n");
//     system("PAUSE");
     }
}

void matrixTest(InteractionMatrix *m, int num)
{
     int i, j, k;
     for(i = 0; i < num ; i++)
     {
           printf("%d) %f %f %f %f %f %f %f %f %f\n", i,m[i][0][0], m[i][0][1], m[i][0][2], m[i][1][0], m[i][1][1], m[i][1][2], m[i][2][0], m[i][2][1], m[i][2][2]);
     }
}

/*
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
*/

float matCalc(float s1[3], InteractionMatrix jMat, float s2[3])
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



float Energy(Atom *atoms, Atom *a, float s1[3], InteractionMatrix* jMatrices)
{
      float E = 0;
      float *s2;
      int j, index;
      for(j = 0; j < a->numInteractions; j++)
      {
              s2 = atoms[a->nbr_list[j]].spin;
              index = a->interaction_matrix[j];
              E -= matCalc(s1, jMatrices[index], s2);
      }
      return E;
}

void flipSpins(Atom *atoms, int numAtoms, InteractionMatrix *jMatrices, float T)
{
     int index;
     float oldE, newE;
     Spin newS;
     for(index = 0; index < numAtoms; index++)
     {
             randSpin(newS);
             oldE = Energy(atoms, atoms+index, atoms[index].spin, jMatrices);
             newE = Energy(atoms, atoms+index, newS, jMatrices);
             if(newE < oldE)
             {
                     //atoms[index].spin = newS;
                     atoms[index].spin[0] = newS[0];
                     atoms[index].spin[1] = newS[1];
                     atoms[index].spin[2] = newS[2];
//                     printf("Spin assigned = (%f, %f, %f)\n", atoms[index].spin[0], atoms[index].spin[1], atoms[index].spin[2]);
             }
             else
             {
                 float deltE = newE - oldE;
                 float probChange = exp(-deltE/T);
                 if(randf() < probChange)
                 {
                            //atoms[index].spin = newS;
                     atoms[index].spin[0] = newS[0];
                     atoms[index].spin[1] = newS[1];
                     atoms[index].spin[2] = newS[2];
//                     printf("Spin assigned = (%f, %f, %f)\n", atoms[index].spin[0], atoms[index].spin[1], atoms[index].spin[2]);
                 }
             }
//             printf("here");
//             system("PAUSE");
             
     }
}

float* getSpins(Atom *atoms, int index)
{
       return atoms[index].spin;
}


//numAtoms may need to be switched to long for long lists of atoms
void simulate(Atom *atoms, int numAtoms, InteractionMatrix *jMatrices, int k, float maxTemp, float minTemp, float tFactor)
{
     //set the seed for the random numbers
//     time_t *time1;
//     srand(123123);
//     printf("here");
      initializeRandState();
     
     int i;
     float T = maxTemp;
     while(T > minTemp)
     {
             for(i = 0; i < k; i++)
             {
//                     printf("flipping spins");
                     flipSpins(atoms, numAtoms, jMatrices, T);
             }
             T = T*tFactor;
     }
     
     //for test purposes
     float avgMag = 0;
     for(i = 0; i < numAtoms; i++)
     {
           avgMag += atoms[i].spin[2];
//           printf("%d) %f (%f, %f, %f)\n", i, avgMag, atoms[i].spin[0], atoms[i].spin[1], atoms[i].spin[2]);
     }
     avgMag = avgMag/numAtoms;
     printf("average Magnetization: %f\n",avgMag);
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


int main(int argc, char *argv[])
{
//    time_t mytime = time(0);
    int i, j;
    dsfmt_state = (dsfmt_t*)malloc(sizeof(dsfmt_t));
    dsfmt_init_gen_rand(dsfmt_state, time(NULL));
    float sum = 0;
    int num = 2000000;
    for(j = 0; j < 20; j++)
    {
        for(i = 0; i < num; i++)
        {
    //        printf("num: %f\n",randFloat(-1,1));
              sum += randFloat(-1,1);
        }
              sum = sum/num;
              printf("%f", sum);
            system("PAUSE");
    }
    
    return 0;
}
