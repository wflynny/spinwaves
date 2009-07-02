#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "dSFMT.h"

typedef float Spin[3];
typedef float Anisotropy[3];
typedef float InteractionMatrix[3][3];

static dsfmt_t dsfmt_state;

float randf()
{
      return (float)dsfmt_genrand_close_open(&dsfmt_state);
}

float randFloat(float min, float max)  //uses randf() to generate a random float between min and max
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
     dsfmt_init_gen_rand(&dsfmt_state, time(NULL));
}


void randSpin(Spin spin, float magnitude)
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
            
       spin[0] = Sx * magnitude;
       spin[1] = Sy * magnitude;
       spin[2] = Sz * magnitude;
       return;
}

typedef struct
{
        int numInteractions;
        int *interaction_matrix;
        int *nbr_list;
//        float *spin;
        Spin spin;
        Anisotropy anisotropy;
        float spinMag;
} Atom;
        
void del_jMat(InteractionMatrix *m)
{
     printf("%x %d", m, m);
      int i, j;
      printf("deleting\n");
      for(i = 0; i < 3; i++)
      {
            for(j = 0; j < 3; j++)
            {
                  printf("%d %d %f \n", i,j,m[0][i][j]);
                  }
                  }
     free(m);
}

InteractionMatrix* new_jMatrix_list(int n, int *success)
{
     InteractionMatrix *m;
     m = (InteractionMatrix*)malloc(n*sizeof(InteractionMatrix));
//     printf("jmatrix start %d  end: %d\n", m, m+n);
 //    system("PAUSE");
     if(m == NULL)
     {
          *success = 0;
     }
     else
     {
         *success = 1;
     }
//     del_jMat(m);
//     printf("deleted\n");
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

void del_atom(Atom *p) 
{
     free(p);
     //for(i = numAtoms-1;  i >0; i--)
     //{ 
     //      printf("%d\n", i);
     //      free(p + i);
    // }
}

Atom* new_atom_list(int n, int *success)
{
     Atom *a;
     a=(Atom*)malloc(n*sizeof(Atom));
//     printf("N = %d  start = %x end = %x\n", n, a, a+n);
 //    system("PAUSE");
     if(a == NULL)
     {
          *success = 0;
     }
     else
     {
         *success = 1;
     }
     return a;//return to calling python method to add the atoms
}


void set_atom(Atom *p, int k, Anisotropy anisotropy, int *mat, int* neighbors, int numInteractions, Spin spin, float spinMag)
{

     int i;

     p[k].interaction_matrix = mat;
     p[k].nbr_list = neighbors;
     p[k].spin[0] = spin[0];
     p[k].spin[1] = spin[1];
     p[k].spin[2] = spin[2];
     p[k].anisotropy[0] = anisotropy[0];
     p[k].anisotropy[1] = anisotropy[1];
     p[k].anisotropy[2] = anisotropy[2];
     p[k].numInteractions = numInteractions;
     p[k].spinMag = spinMag;
     
}

void getSpin(Atom *p, int index, Spin spin)
{
       spin[0] = p[index].spin[0];
       spin[1] = p[index].spin[1];
       spin[2] = p[index].spin[2];
       return;
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



float matCalc(float *s1, InteractionMatrix jMat, float s2[3])
{

      //multiply s1 by jMat
      float x1 = (jMat[0][0] * s1[0]) + (jMat[0][1] * s1[1]) + (jMat[0][2] * s1[2]);
      float y1 = (jMat[1][0] * s1[0]) + (jMat[1][1] * s1[1]) + (jMat[1][2] * s1[2]);
      float z1 = (jMat[2][0] * s1[0]) + (jMat[2][1] * s1[1]) + (jMat[2][2] * s1[2]);

      //return the dot product
      return ((x1*s2[0]) + (y1*s2[1]) + (z1*s2[2]));
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
//      printf("s1: %x (%f, %f, %f)  atom* = %x\n", s1, s1[0],s1[1], s1[2], a);
      for(j = 0; j < a->numInteractions; j++)
      {
              s2 = atoms[a->nbr_list[j]].spin;
//              printf("s2: %x (%f, %f, %f)", s2, s2[0], s2[1], s2[2]);
              index = a->interaction_matrix[j];
//              printf("here2");
                //E = -s1*Jij*s2 - Dx*Sx^2 - Dy*Sy^2 - DzSz^2
              E -= matCalc(s1, jMatrices[index], s2) + a->anisotropy[0]*pow(s1[0],2) + a->anisotropy[1]*pow(s1[1],2) + a->anisotropy[2]*pow(s1[2],2);
//              anisotropyE = a->anisotropy[0]*pow(s1[0],2) + a->anisotropy[1]*pow(s1[1],2) + a->anisotropy[2]*pow(s1[2],2);
//              printf("E = %f\n", E);
//              E -= anisotropyE;
//              printf("anisotropy: (%f, %f, %f)  s: (%f, %f, %f) E = %f\n", a->anisotropy[0], a->anisotropy[1], a->anisotropy[2],s1[0],s1[1],s1[2],E);
//              system("PAUSE");
//              printf("here3\n");
      }
//      printf("ok E = %f\n", E);
      
      return E;
}

void flipSpins(Atom *atoms, int numAtoms, InteractionMatrix *jMatrices, float T, int *flipped)
{
     int index;
     float oldE, newE;
     Spin newS;
     for(index = 0; index < numAtoms; index++)
     {
             randSpin(newS, atoms[index].spinMag);
             oldE = Energy(atoms, atoms+index, atoms[index].spin, jMatrices);

             newE = Energy(atoms, atoms+index, newS, jMatrices);

             if(newE < oldE)
             {
                     //atoms[index].spin = newS;
                     atoms[index].spin[0] = newS[0];
                     atoms[index].spin[1] = newS[1];
                     atoms[index].spin[2] = newS[2];
                     (*flipped)++;
//                     printf("Spin assigned = (%f, %f, %f)\n", atoms[index].spin[0], atoms[index].spin[1], atoms[index].spin[2]);
             }
             else
             {
                 float deltE = newE - oldE;
//                 printf("deltE : %f  ", deltE);
                 float probChange = exp(-deltE/T);
//                 printf("probChange : %f\n", probChange);
                 if(randf() < probChange)
                 {
                            //atoms[index].spin = newS;
                     atoms[index].spin[0] = newS[0];
                     atoms[index].spin[1] = newS[1];
                     atoms[index].spin[2] = newS[2];
                     (*flipped)++;
//                     printf("Spin assigned = (%f, %f, %f)\n", atoms[index].spin[0], atoms[index].spin[1], atoms[index].spin[2]);
                 }
             }
//             printf("here");
//             system("PAUSE");
             
     }
}

/*
float* getSpins(Atom *atoms, int index)
{
       return atoms[index].spin;
}
*/

//numAtoms may need to be switched to long for long lists of atoms
void simulate(Atom *atoms, int numAtoms, InteractionMatrix *jMatrices, int k, float maxTemp, float minTemp, float tFactor)
{
     float fracFlipped;
     int i,j, flipped;
     
     for(i = 0; i < numAtoms; i++)
     {
           printf("Anisotropy: (%f, %f, %f)\n", atoms[i].anisotropy[0], atoms[i].anisotropy[1], atoms[i].anisotropy[2]);
     }   
     
     
     printf("\nnumAtoms = %d\n", numAtoms);
     //set the seed for the random numbers
     initializeRandState();
     
     float T = maxTemp;
     while(T > minTemp)
     {
             flipped = 0;
             for(i = 0; i < k; i++)
             {
                     flipSpins(atoms, numAtoms, jMatrices, T, &flipped);
             }
             fracFlipped = ((float)flipped)/(k * numAtoms);
             printf("Temperature: %f  Fraction flipped: %f\n", T, fracFlipped);
             T = T*tFactor;
     }
 
     return;
}



int main(int argc, char *argv[])
{
//    time_t mytime = time(0);
    int i, j;
    initializeRandState();
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
