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
/*
void randSpin(Spin spin)
{
     float random = randf();
     if(random < 0.5)
     {
      spin[2] = -1;
     }
     else
     {
         spin[2] = 1;
     }
//     printf("Spin: %f\n", spin[2]);
}
*/


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
        Anisotropy anisotropy;
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


void set_atom(Atom *p, int k, Anisotropy anisotropy, int *mat, int* neighbors, int numInteractions, Spin spin)
{
//     printf("spin: (%f, %f, %f) ", spin[0],spin[1], spin[2]);
//     printf("%d) %x  %x  %x\n", k, p, p + k, p + k + 1);
       int i;
/*       printf("(%d)\n", k);
       for(i =0; i < numInteractions;  i++)
       {
             printf("%d)(%d, %d)\n",i, neighbors[i], mat[i]);
       }
*/
//       system("PAUSE");
     p[k].interaction_matrix = mat;
     p[k].nbr_list = neighbors;
//     printf("here1\n");
     p[k].spin[0] = spin[0];
//     printf("here2\n");
     p[k].spin[1] = spin[1];
//     printf("here3\n");
     p[k].spin[2] = spin[2];
     p[k].anisotropy[0] = anisotropy[0];
     p[k].anisotropy[1] = anisotropy[1];
     p[k].anisotropy[2] = anisotropy[2];
//     printf("anisotropy: (%f, %f, %f)\n", anisotropy[0], anisotropy[1], anisotropy[2]);
//     system("PAUSE");
//     printf("here4\n");
     p[k].numInteractions = numInteractions;
     
/*     printf("%d)", k);
     int i;
     for(i = 0; i < numInteractions; i++)
     {
           printf("%d %d\n", mat[i], neighbors[i]);
     }
     system("PAUSE");*/
     
/*     int i, j;
     printf("K = %d\n", k);
     for(i = 0; i < k+1; i++)
     {
           printf("\n%d) spin: (%f, %f, %f)", i, p[i].spin[0], p[i].spin[1], p[i].spin[2]);
           for(j = 0; j < p[i].numInteractions; j++)
           {
                 printf("(%d, %d) ", p[i].nbr_list[j], p[i].interaction_matrix[j]);
           }
           system("PAuSE");
     }*/
     
     //test
//     printf("mat* = %d   mat[0],[1],[2] = (%d,%d,%d)   p[k].mat* = %d p[k].mat[0],[1],[2] = (%d,%d, %d)\n", mat, mat[0], mat[1], mat[2], p[k].interaction_matrix, p[k].interaction_matrix[0],p[k].interaction_matrix[1], p[k].interaction_matrix[2]);
}

void getSpin(Atom *p, int index, Spin spin)
{
//       printf("(%f, %f, %f)\n", p[index].spin[0], p[index].spin[1], p[index].spin[2]);
       spin[0] = p[index].spin[0];
       spin[1] = p[index].spin[1];
       spin[2] = p[index].spin[2];
       return;
}
/* implemented in python instead
float getMagnetization(Atom *atoms, int numAtoms)
{
      //Returns the length of the sum of all the spins
      float avgMag = 0;
      float Xsum = 0;
      float Ysum = 0;
      float Zsum = 0;
      int i;
      for(i = 0; i < numAtoms; i++)
      {
            Xsum += atoms[i].spin[0];
            Ysum += atoms[i].spin[1];
            Zsum += atoms[i].spin[2];
      }
      avgMag = sqrt(pow(Xsum,2) + pow(Ysum,2) + pow(Zsum,2));
      return avgMag;
}
*/
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

float matCalc(float *s1, InteractionMatrix jMat, float s2[3])
{
 //     printf("%f, %f, %f\n", s1[0], s1[1], s1[2]);
//      system("PAUSE");
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
//     printf("flipping spins");
     for(index = 0; index < numAtoms; index++)
     {
//             printf("atom %d) s = (%f,%f,%f) \n", index, atoms[index].spin[0], atoms[index].spin[1], atoms[index].spin[2]);
             randSpin(newS);
             oldE = Energy(atoms, atoms+index, atoms[index].spin, jMatrices);
//             printf("oldE : %f\n", oldE);

//             printf("news = (%f, %f, %f)\n", newS[0], newS[1], newS[2]);
             newE = Energy(atoms, atoms+index, newS, jMatrices);
//             printf("newE : %f\n", newE);
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
     
/*     printf("C side interaction matrix list %d\n", jMatrices);
     for (i=0; i<3; i++)
     for (j=0; j<3; j++){
         printf("%d %d Jij=%f\n",i,j,jMatrices[0][i][j]);
         }
*/    
     
     
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
 


     
     //for test purposes
     float avgMag = 0;
     for(i = 0; i < numAtoms; i++)
     {
           avgMag += atoms[i].spin[2];
//           printf("%d) %f (%f, %f, %f)\n", i, avgMag, atoms[i].spin[0], atoms[i].spin[1], atoms[i].spin[2]);
     }
     avgMag = avgMag/numAtoms;
     printf("average Magnetization: %f\n",avgMag);

     return;
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