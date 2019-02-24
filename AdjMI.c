#include <stdarg.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include<assert.h>
#include<MutualInformation.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <unistd.h>
#include <err.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
//#include<MIToolbox/MutualInformation.h>

//#define maxSNPCount 2000
//#define acc         1300

//#define DEBUG 1
#define PARALLEL 1

//This is for controling the parallel part. If commented the parallel doesn't work

const unsigned long maxSNPCount=11000000;
const unsigned long acc=1300;
const int           maxCores=40;
const int           ASCII_LIMIT=200;


#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)

/* Dependency: 
https://github.com/Craigacp/MIToolbox/
   make, then make install as root, had to libMIToolbox.so to /usr/lib64

compile with:
change #include<MIToolsbox/MutualInformation.h>
gcc -O3 AdjMI.c -lMIToolbox -lm -o AdjMI.exe

##Golem
module add lib/MIToolbox
gcc -g AdjMI.c -I/apps/lib/MIToolbox/include -lMIToolbox -lm -o AdjMI.exe

#Draco
module load gcc
export MIToolbox=~/MIToolbox
gcc -g AdjMI.c -I$MIToolbox -L$MIToolbox -lMIToolbox -lm -o AdjMI.exe -Wl,-rpath $MIToolbox

#Minerva
module add MIToolbox/gcc
gcc -g AdjMI.c -I/cluster/apps/MIToolbox/gcc/20161213/include/MIToolbox -L/cluster/apps/MIToolbox/gcc/20161213/lib -lMIToolbox -lm -o AdjMI.exe
//The include and the libraries are set up differenly (the organization of the directiories)

The for loops of this program are in many occasions written to be used in biallelic SNPs.
This program is meant to be run in parallel using share memory to read the input file.
Generates the intersection between pairs avoiding Ns. Needs a minimum of accession to create the pairs (Threshold).
Only prints out Adjusted Mutual Information that is over 0.9. It Doesn't take into account polymorphic sites which 
its minor allele frequeny is lower than the 5% of the length of the polymorphic site. Only works for biallelics.

*/

//-----------------------------------------------------------------Functions------------------------------------------------------//

//Never count longer than the length of the array!! NO num_elements+1
int max_array(int *a, int num_elements) {

   int g, max=a[0]; 

   for (g=0; g<num_elements; g++){
      if (a[g]>max) { 
         max=a[g]; 
      } 
   } 
   return max;
}

int min_array(int *b, int num_elements) {

   int g, min=b[0]; 

   for (g=0; g<num_elements; g++){
      if (b[g]<min) { 
         min=b[g]; 
      } 
   } 
   return min;
}

int array_length(int first, int last) { 
   int l; 
   l=(last-first)+1; 
   return l;
}

void array_seq(int first, int last, int array[]) { 
    int l,e; 
    int start=first; 
    l=(last-first)+1; 

    for(e=0;e<l;e++) { 
       array[e]=start++; 
    } 
}

//#define max_array(ret, a, num_elements) {int g, max=a[0]; for (g=0; g<num_elements+1; g++){if (a[g]>max) { max=a[g]; } } ret=max;}
//#define min_array(ret, b, num_elements) {int g, min=b[0]; for (g=0; g<num_elements+1; g++){if (b[g]<min) { min=b[g]; } } ret=min;}
//#define array_length(ret,first,last) { int l; l=(last-first)+1; ret=l;}
//#define array_seq(first, last, array) { int l,e; int start=first; l=(last-first)+1; for(e=0;e<l;e++) { array[e]=start++; } }

//It is better for debugging to use the complete functions, because the debagger goes step by step
static void
check (int test, const char * message, ...)
{
    if (test) {
        va_list args;
        va_start (args, message);
        vfprintf (stderr, message, args);
        va_end (args);
        fprintf (stderr, "\n");
        exit (EXIT_FAILURE);
    }
}

int Intersection(int arr1[], int arr2[], int m, int n, int array[])
{
   //printf("Intersection ");
   int i=0, j=0, index=0;

   while (i < m && j < n) {

      if (arr1[i] < arr2[j])
      	i++;
      else if (arr2[j] < arr1[i])
      	j++;
      else{
        array[index++]=arr2[j++];
        i++;
      }
   }

   return index;
}



//Function to check if a value is in an array// //This function didn't make much sense. Use instate: lessminor1= ( entropy1[0] <= minorallele || entropy1[1] <= minorallele ) ? 1 : 0;
int exist(int val, int *arr, int size){

    int i,answer=0;
    for (i=0; i < size; i++) {
        if (arr[i]<=val){
	    answer++;
	}
    }
	return answer;
}

//Function to calculate AdjMI 
void  ChildProcess(int pid, int start,int stop,char (*alleles)[acc],char inputFileName[300],int *pos, int *chr,int numberSNPs,int Threshold)//tengo que pasarle a la function todo loque ya haya definido 
{
   FILE    *out;
   char    outputFileName[200];
   int     n1,n2,i,j,k,g,f,h,s,z; 
   int     end,X1,X2,X3,aa,bb,nij;

   int     frequency1_A, frequency1_T, frequency1_C, frequency1_G;
   int     frequency2_A, frequency2_T, frequency2_C, frequency2_G;

   int     nomlength,demlength,nom1length,nom2length,nom3length,dem1length,dem2length,dem3length;
   int     *nom,*dem;

   double  miself1=0,miself2=0;

   int     *array1,*array2;
   int     *v1,*v2;
   double  *V1,*V2;
   int     t;
   int     entropy1[2],entropy2[2];

   double  *rp1;
   double  mi;
   double  Ha=0.0,Hb=0.0, AMI=0.0,NMI=0.0,EMI=0.0,EMI_bound=0.0;
   double  *logNij;

   int indexarray1, indexarray2;
/* mehlert
   double  E3[2][2];
   double  EPLNP[2][2];
   double  bound[2][2];
*/
   int     count;
   int     *array;

   double  multi;
   double  minorallele;

   double  CC;

//   int frequency1[ASCII_LIMIT],frequency2[ASCII_LIMIT];

//   int n3;
   int numberOfDiffAlleles1;
   int numberOfDiffAlleles2;

   int lessminor1;
   int lessminor2;
   int both;

   double prod;

   int *nom1;
   int *nom2;
   int *nom3;
   int *dem1;
   int *dem2;
   int *dem3;

    time_t timer;
    char buffer[26];
    struct tm* tm_info;

   time(&timer);
   tm_info = localtime(&timer);
   strftime(buffer, 26, "%Y_%m_%d %H:%M:%S", tm_info);

   v1     = calloc(acc,sizeof(int));
   v2     = calloc(acc,sizeof(int));
   V1     = calloc(acc,sizeof(double));
   V2     = calloc(acc,sizeof(double));

   array1 = calloc(acc,sizeof(int));
   array2 = calloc(acc,sizeof(int));
   array  = calloc(2*acc,sizeof(int));


   sprintf(outputFileName,"%s_%i_%i_AdjMIresults.txt",inputFileName,start,stop);
   out=fopen(outputFileName,"w");

#ifdef DEBUG
   fprintf(out,"DEBUG %s ChildProcess() start \n",buffer);
#endif

   fprintf(out,"EMI\tMI\tAMI\tNMI\tEntropyA\tEntropyB\tCHR_A\tPOS_A\tFreqA1\tFreqA2\tCHR_B\tPOS_B\tFreqB1\tFreqB2\tLength\tAccCode\t\n"); fflush(out);

   for(i=start;i<stop;i++) {

#ifdef DEBUG
      time(&timer);
      tm_info = localtime(&timer);

      strftime(buffer, 26, "%Y_%m_%d %H:%M:%S", tm_info);

      fprintf(out,"DEBUG %s i=%i\n",buffer,i); fflush(out);
#endif

      // only interval
      // i = line
      // k = pos in line
      
      n1     = strlen(alleles[i]);

      indexarray1 = 0;

      for(k=0; k<n1; k++) {

         if( alleles[i][k] != 78 ){//ASCII code for 'N' 

            v1[k]=alleles[i][k];
       	    array1[indexarray1++]=k;

//            printf("DEBUG alleles[%i][%i]=%c v1[%i]=%i \n",i,k,alleles[i][k],k,v1[k]);
	 }	
      }

      for(j=i+1; j<numberSNPs; j++) {

//          printf("DEBUG %i\n",j);
         
//I need to clean these every time to avoid corrupted numbers. Like 'free' for statics  
          double E3[2][2]    = {{0,0},{0,0}};
          double EPLNP[2][2] = {{0,0},{0,0}};
          double bound[2][2] = {{0,0},{0,0}};
 
//         printf("DEBUG j=%i\n",j);

         n2=strlen(alleles[j]);

         if(n1==n2) {

	    indexarray2=0;

	    for(k=0;k<n2;k++) {

               if( alleles[j][k] != 78 ){

       	          v2[k]=alleles[j][k];
       		  array2[indexarray2++]=k;

	       }
	    }

	    Intersection(array1,array2,indexarray1,indexarray2,array);//con este salvo el array de resultado sin problemas de alloc

	    count = Intersection(array1,array2,indexarray1,indexarray2,array);//con este salvo la longitud del arreglo

	    //fprintf(out,"\nIntersections: ");
	    //int xx;for(xx=0;xx<count;xx++){fprintf(out,"%i ",array[xx]);}

	    multi       = (double)5/100;//Cuando se divide debe obligatoriamente el denominador ser double
	    minorallele = count*multi;

	    //fprintf(out,"\n\ncount: %i minoralle: %e\n Vectores ",count,minorallele);

	    if(count>Threshold) {

               frequency1_A=0;
               frequency1_T=0;
               frequency1_C=0;
               frequency1_G=0;
               frequency2_A=0;
               frequency2_T=0;
               frequency2_C=0;
               frequency2_G=0;
//I was wondering why V1 and V2 don't have to be allocated inside the loop to work properly, I mean to avoid acumulating values from the last iterations. 
//The key is here I visit the new array only with the length of count, meaning that even if the array has some addtional values from the last run I wont visit them to fill V1/V2.
//These values have been replace in this iteration and then there is no a problem. The problem comes when I visit more position that this iteration haven't been able to fill.
//In that cases the array (when no reinitialized) will have the values from the last iteration. 
	       for(z=0;z<count;z++) {

                  V1[z]=v1[array[z]];
                  V2[z]=v2[array[z]];

                  //fprintf(out," %e  ",V2[z]);
                  //fprintf(out," %e  ",V1[z]);  
// mehlert		  frequency1[(int)V1[z]]++;
// mehlert                frequency2[(int)V2[z]]++;

//                  printf("DEBUG V1[%i]=%e V2[%i]=%e \n",z,V1[z],z,V2[z]);
//This is a much better way to do the counting.
                  switch ( (int)V1[z] ){
                   case 65: frequency1_A++; break;
                   case 84: frequency1_T++; break;
                   case 67: frequency1_C++; break;
                   case 71: frequency1_G++;
                  }

                  switch ( (int)V2[z] ){
                   case 65: frequency2_A++; break;
                   case 84: frequency2_T++; break;
                   case 67: frequency2_C++; break;
                   case 71: frequency2_G++; 
                  }


               }

//               n3=127;
               numberOfDiffAlleles1=0;
               numberOfDiffAlleles2=0;

/* mehlert               for (l=0;l<n3;l++) {

                  if(frequency1[l]>0) { 
                     entropy1[numberOfDiffAlleles1++]=frequency1[l];
                  }
               }
*/
//Clean or reinitialize. Avoid corrected numbers. Also is importan for additive values '++', next loop could keep adding if not reinitialized
                 entropy1[0]=0; entropy1[1]=0;
                 entropy2[0]=0; entropy2[1]=0;

                 if (frequency1_A) entropy1[numberOfDiffAlleles1++]=frequency1_A;
                 if (frequency1_T) entropy1[numberOfDiffAlleles1++]=frequency1_T;
                 if (frequency1_C) entropy1[numberOfDiffAlleles1++]=frequency1_C;
                 if (frequency1_G) entropy1[numberOfDiffAlleles1++]=frequency1_G;

                 if (frequency2_A) entropy2[numberOfDiffAlleles2++]=frequency2_A;
                 if (frequency2_T) entropy2[numberOfDiffAlleles2++]=frequency2_T;
                 if (frequency2_C) entropy2[numberOfDiffAlleles2++]=frequency2_C;
                 if (frequency2_G) entropy2[numberOfDiffAlleles2++]=frequency2_G;
                  
                  //int xx;for(xx=0;xx<3;xx++){fprintf(out,"position %i entropy %i",xx,entropy1[xx]);}



 /* mehlert                 for (l=0;l<n3;l++) {

                     if(frequency2[l]>0) {
                        entropy2[numberOfDiffAlleles2++]=frequency2[l];
                        //fprintf(out,"\n\n 2: item %i position %i\n\n",frequency2[l],numberOfDiffAlleles2);
                     }
                  }
*/

//                  lessminor1 = exist((int)minorallele,entropy1,2);
//                  lessminor2 = exist((int)minorallele,entropy2,2);

                  lessminor1= ( entropy1[0] <= minorallele || entropy1[1] <= minorallele ) ? 1 : 0;
                  lessminor2= ( entropy2[0] <= minorallele || entropy2[1] <= minorallele ) ? 1 : 0;

                  //fprintf(out,"\n\nV1:MINOR %i FA1 %i FA2 %i IsFA %i\n",(int)minorallele,entropy1[0],entropy1[1],lessminor1);
                  //fprintf(out,"\n\nV2:MINOR %i FA1 %i FA2 %i IsFA %i\n",(int)minorallele,entropy2[0],entropy2[1],lessminor2);

                  both=lessminor1+lessminor2;
                  //fprintf(out,"\n\n BOTH %i\n\n",both);

		  if(both<1) {
                     Ha = -(((double)entropy1[0]/count)*(log((double)entropy1[0]/count)/log(2))+((double)entropy1[1]/count)*(log((double)entropy1[1]/count))/log(2));
                     Hb = -(((double)entropy2[0]/count)*(log((double)entropy2[0]/count)/log(2))+((double)entropy2[1]/count)*(log((double)entropy2[1]/count))/log(2));

                     

                     mi = calculateMutualInformation(V1,V2,count);

                     miself1=calculateMutualInformation(V1,V1,count);
                     miself2=calculateMutualInformation(V2,V2,count);
                     //fprintf(out,"i:%i j:%i  %e = %e  %e = %e \n  %i   %i  %i   %i\n",i,j,Ha,miself1,Hb,miself2,entropy1[0],entropy2[0],entropy1[1],entropy2[1]); fflush(out);
                     int dx=entropy1[0]*entropy2[0];
                     int dy=entropy1[0]*entropy2[1];
                     int dw=entropy1[1]*entropy2[0];
                     int dz=entropy1[1]*entropy2[1];
                     int AB[2][2]= { {dx,dy},{dw,dz} };

                     for (s=0;s<2;s++) {
                        for (g=0;g<2;g++) {
                           E3[s][g]=(double)AB[s][g]/pow(count,2)*( ( log( (  (double)AB[s][g]  /  pow( count,2 ) ) ) / log(2)));
                        }
                     }

                     int MaxA = max_array(entropy1,2);
                     int MaxB = max_array(entropy2,2);
                     int Nij[2]={MaxA, MaxB};
                     end  = min_array(Nij,1);

                     logNij=calloc(end+2,sizeof(double));

                     for (t=0;t<end+1;t++) {
                        logNij[0]=0;
                        logNij[t+1]=(log(((double)(t+1)/count)))/log(2);
                     }

                     for(h=0;h<2;h++) {

                        for(f=0;f<2;f++) {

                           int q=((entropy1[h]+entropy2[f])-count);
                           int p[2]={1,q};

                           nij=max_array(p,2);

                           int r=count-entropy1[h]-entropy2[f]+nij;
                           int o[2]={nij,r};

                           if(o[0]>o[1]) {
                              int tmp=o[0]; o[0]=o[1]; o[1]=tmp;
                           }

                           aa=entropy1[h]-nij+1;
                           bb=entropy2[f]-nij+1;
                           X1=o[1]+1;
                           X3=count-entropy1[h]+1; 
                           X2=count-entropy2[f]; 

                           if(X2>o[1]) {
                              nom1length = array_length(aa,entropy1[h]);
                              nom2length = array_length(bb,entropy2[f]);
                              nom3length = array_length(X1,X2);
                              dem1length = array_length(X3,count);
                              dem2length = array_length(1,o[0]);

//When a function arr[] is equal to *arr. Both are pointers. Do no make statics arrays that by the time of compiling their size is unknown. Better use calloc in these cases
                              nom1 = calloc(nom1length,sizeof(int));
                              nom2 = calloc(nom2length,sizeof(int));
                              nom3 = calloc(nom3length,sizeof(int));
                              dem1 = calloc(dem1length,sizeof(int));
                              dem2 = calloc(dem2length,sizeof(int));

                              array_seq(aa,entropy1[h],nom1);
                              array_seq(bb,entropy2[f],nom2);
                              array_seq(X1,X2,nom3);
                              array_seq(X3,count,dem1);
                              array_seq(1,o[0],dem2);

                              nomlength=nom1length+nom2length+nom3length;
                              demlength=dem1length+dem2length;

                              nom=calloc(nomlength,sizeof(int));
                              dem=calloc(demlength,sizeof(int));

                              memcpy(nom,nom1,nom1length * sizeof(int));
// mehlert                              memcpy(&nom[nom1length], nom2, nom2length * sizeof(int));
                              memcpy(nom+nom1length, nom2, nom2length * sizeof(int));
// mehlert                              memcpy(&nom[nom1length + nom2length], nom3, nom3length * sizeof(int));
                              memcpy(nom+nom1length + nom2length, nom3, nom3length * sizeof(int));
                              // int m; for(m=0;m<nom1length+nom2length+nom3length;m++) {printf("-%i",nom[m]);} printf("\n");
                              memcpy(dem,dem1,dem1length * sizeof(int));
// mehlert                    memcpy(&dem[dem1length], dem2, dem2length * sizeof(int));
                              memcpy(dem+dem1length, dem2, dem2length * sizeof(int));
//Using & is like visiting all your friends on the way home
                              free(nom1);
                              free(nom2);
                              free(nom3);
                              free(dem1);
                              free(dem2);

                           } else {

                              nom1length = array_length(aa,entropy1[h]);
                              nom2length = array_length(bb,entropy2[f]);
                              dem1length = array_length(X3,count);
                              dem2length = array_length((X2+1),o[1]);
                              dem3length = array_length(1,o[0]);

                              nom1 = calloc(nom1length,sizeof(int));
                              nom2 = calloc(nom2length,sizeof(int));
                              dem1 = calloc(dem1length,sizeof(int));
                              dem2 = calloc(dem2length,sizeof(int));
                              dem3 = calloc(dem3length,sizeof(int));

                              array_seq(aa,entropy1[h],nom1);
                              array_seq(bb,entropy2[f],nom2); 
                              array_seq(X3,count,dem1);
                              array_seq((X2+1),o[1],dem2);
                              array_seq(1,o[0],dem3); 

                              nomlength=nom1length+nom2length;
                              demlength=dem1length+dem2length+dem3length;

                              nom=calloc(nomlength,sizeof(int));
                              dem=calloc(demlength,sizeof(int));

                              memcpy(nom,nom1,nom1length * sizeof(int));
// mehlert                    memcpy(&nom[nom1length], nom2, nom2length * sizeof(int));
                              memcpy(nom+nom1length, nom2, nom2length * sizeof(int));
                              //  int m; for(m=0;m<nom1length+nom2length;m++) {printf("-%i",nom[m]);} printf("\n");
                              memcpy(dem,dem1,dem1length * sizeof(int)); 
// mehlert                    memcpy(&dem[dem1length], dem2, dem2length * sizeof(int));
                              memcpy(dem+dem1length, dem2, dem2length * sizeof(int));
// mehlert                    memcpy(&dem[dem1length + dem2length], dem3, dem3length * sizeof(int));
                              memcpy(dem+dem1length+dem2length, dem3, dem3length * sizeof(int));

                              free(nom1);
                              free(nom2);
                              free(dem1);
                              free(dem2);
                              free(dem3);
                           }

                           int d;

                           rp1=calloc(nomlength+2,sizeof(double));

                           prod = 1.0; 

                           for(d=0;d<nomlength;d++){

                              if(nom[d]>0){
                                 rp1[d]=(double)nom[d]/dem[d];
                                 prod= prod*rp1[d];
                                 //printf("\n Div %e Prod %e %i %i\n ",rp1[d],prod,nom[d],dem[d]);
                              }
                           }	

                           //This needs to be double because the multiplication at some point is so big that a float can't handle. 
                           prod=prod/count;

                           int nij1=entropy2[f]+entropy1[h]-count;
                           int nij2_array[2]={entropy2[f],entropy1[h]};
                           int nij1_array[2]={1,nij1};

                           int primero;
                           primero = max_array(nij1_array,2); 

                           int ultimo;
                           ultimo = min_array(nij2_array,1);

                           int nij_length;
                           nij_length = array_length(primero,ultimo);  // is the correct order

                           int nij_seq[nij_length];
                           array_seq(primero,ultimo,nij_seq);

                           int x;

                           for(x=0;x<nij_length;x++) {
                              //printf("EPLNP %e prod %e logNij %e nij %i x %i\n",EPLNP[h][f],prod,logNij[(int)nij_seq[x]],nij_seq[x],x);
                              EPLNP[h][f]=(EPLNP[h][f]+nij_seq[x]*logNij[nij_seq[x]]*prod);

                              //printf("EPLNP %e prod %e logNij %e nij %i\n",EPLNP[h][f],prod,logNij[(int)nij_seq[x]],nij_seq[x]);
                              prod=(double)(prod*(entropy2[f]-nij_seq[x])*(entropy1[h]-nij_seq[x]))/(double)(nij_seq[x]+1)/(count-entropy2[f]-entropy1[h]+nij_seq[x]+1);    

                              //printf("EPLNP %e prod %e logNij %e nij %i   %i %i %i %i %i\n",EPLNP[h][f],prod,logNij[(int)nij_seq[x]],nij_seq[x],h,f,x,entropy2[h],entropy1[f]);
                           }

                           prod=1.0;
                           CC=(double)((double)((double)(count*(entropy2[f]-1)*(entropy1[h]-1))/(double)entropy2[f])/entropy1[h])/((count-1))+(((double)count/entropy2[f])/entropy1[h]);
                           //bound[h][f]=(float)( ( (double)entropy2[f]*(double)entropy1[h] ) / (count*count)* (log(CC)/ log(2)) );

                           int uu =  entropy2[f] * entropy1[h] ;
                           double ll1,ll2;
                           ll1 = log(2);
                           ll2 = log(CC);
                           double  ll =  (double)ll1 / ll2;
                           double  dd = count * count * ll;
                           bound[h][f]= (double)uu / dd;
                           //EMI_bound=EMI_bound+bound[h][f];
                           //EMI=EMI+(EPLNP[h][f]-E3[h][f]);
                           //printf("\n\n EMI %e EPLNP %e E3 %e\n\n",EMI,EPLNP[h][f],E3[h][f]);

                           free(nom);
                           free(dem);
                           free(rp1);
                        } 
                     }

                     //fprintf(out,"\n\n%e %e %e %e %e %e %e %e\n\n",EPLNP[0][0],E3[0][0],EPLNP[0][1],E3[0][1],EPLNP[1][1],E3[1][1],EPLNP[1][0],E3[1][0]);
                     EMI_bound=bound[0][0]+bound[1][0]+bound[0][1]+bound[1][1];
                     EMI=((EPLNP[0][0])-(E3[0][0]))+((EPLNP[0][1])-(E3[0][1]))+((EPLNP[1][1])-(E3[1][1]))+((EPLNP[1][0])-(E3[1][0]));

                     AMI=(double)(mi-EMI)/(sqrt((double)(Ha*Hb))-EMI);
                     NMI=(double)mi/sqrt((double)(Ha*Hb));

                     if (abs(EMI)>EMI_bound){
// DEBUG                        printf("The EMI is small: EMI < %f, setting AMI=NMI\n",EMI_bound);
                        AMI=NMI;
                     }

                      if(AMI>0.9) {
//                     if(AMI>0) {
                        fprintf(out,"%e\t%e\t%e\t%e\t%e\t%e\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\t",EMI,mi,AMI,NMI,Ha,Hb,chr[i],pos[i],(int)entropy1[0],(int)entropy1[1],chr[j],pos[j],(int)entropy2[0],(int)entropy2[1],(int)count);
                        int m; for(m=0;m<count;m++) {fprintf(out,"%i,",array[m]);} fprintf(out,"\n");
                        //printf("MI:%e\tAMI:%e\tNMI:%e\tEMIe=%e\tEMIfloat=%f\t%i\t%i\t%i\t%i\n",mi,AMI,NMI,EMI,EMI,chr[i],pos[i],chr[j],pos[j]);
                        //fprintf(out,"%i\t%i\t%e\t%e\t%i\t%i\t%e\t%e\t%e\n",entropy1[0],entropy1[1],Ha,miself1,entropy2[0],entropy2[1],Hb,miself2,EMI);/
                        //fprintf(out,"%i\t%i\t%e\t%i\t%i\t%e\t%e\n",entropy1[0],entropy1[1],Ha,entropy2[0],entropy2[1],Hb,EMI);
                        //fprintf(out,"%i\t%i\t%e\t%i\t%i\t%e\t%e\n",entropy1[1],entropy1[0],Ha,entropy2[1],entropy2[0],Hb,EMI);
                     } 

                     free(logNij);

		  } // if(both<1)

            } // if(count>Threshold)

         } // if(n1==n2)


      } // for(j=i+1; j<numberSNPs; j++)


   } // for(i=start;i<stop;i++)

   free(v1);
   free(v2);
   free(V1);
   free(V2);

   free(array1);
   free(array2);
   free(array);

#ifdef DEBUG
   time(&timer);
   tm_info = localtime(&timer);
   strftime(buffer, 26, "%Y_%m_%d %H:%M:%S", tm_info);

   fprintf(out,"DEBUG %s ChildProcess() return \n",buffer); fflush(out);
#endif

} // ChildProcess()


	void fill(char **p,int rows, int cols)
	  {
		memset(p, (char)' ', cols*rows);
	   }

//-----------------------------------------------------------main------------------------------------------------------//


int main(int argc, char **argv)
{
        FILE    *in;
        char    *line;
        char    inputFileName[300],Intervals[300],outputName[300];
        int     *chr,*pos;
        int     numberSNPs;
        int     linealloc=1500;
        int     *interval_start,*interval_end,Threshold;
        char (*alleles)[acc];

//-----------------------------------------------------------------Reading Arguments--------------------------------------------------------------------//

	if(argc==4) {

      	   sprintf(inputFileName,"%s",argv[1]);
      	   sprintf(Intervals,"%s",    argv[2]);

	   Threshold=atoi(argv[3]);

	   printf("\nInput File = %s\n",inputFileName);
    	   printf("\nIntervals = %s\n",Intervals);
	   printf("\nThreshold = %i\n",Threshold);

	} else {
      	   printf("Wrong number of parameters!\n");
      	   printf("run with: <program> <inputFile> <Intervals> <Threshold> \n");
      	   printf("output will be written to file inputFileName_AdjAMIresults.txt.\n");
	   exit(0);
    	}

//-----------------------------------------------------------------Reading The Intervals--------------------------------------------------------------------//

	FILE *inter_file;
        char *buf1;
	int   index1=0;

	buf1= calloc(1500,sizeof(char));

        inter_file     = fopen(Intervals,"r");
	interval_start = calloc(maxCores,sizeof(int));
	interval_end   = calloc(maxCores,sizeof(int));

        while (fgets(buf1,1500,inter_file)!=NULL) {

           sscanf(buf1,"%i %i",&interval_start[index1],&interval_end[index1]);
           //printf("Intervals table %i\t%i\n",interval_start[index1],interval_end[index1]);
           index1++;

        }

        fclose(inter_file);

//---------------------------------------------------------------------------Reading input File----------------------------------------------------------//
	//Share memory

	int         fd;
   	 /* Information about the file. */
    	struct      stat s1;
    	int         status;
//   	size_t      size;
    
    	/* Open the file for reading. */
    	fd = open (inputFileName, O_RDONLY);
    	check (fd < 0, "open %s failed: %s", inputFileName, strerror (errno));

    	/* Get the size of the file. */
    	status = fstat (fd, & s1);
    	check (status < 0, "stat %s failed: %s", inputFileName, strerror (errno));
//    	size = s1.st_size;

	long unsigned int ARRAY_SIZE_alleles = maxSNPCount*(acc+1);
	long unsigned int ARRAY_SIZE_pos     = maxSNPCount*sizeof(int);
	long unsigned int ARRAY_SIZE_chr     = maxSNPCount*sizeof(int);

	/*el numero del segundo array debe ser menor que el ARRAY_SIZE para el mmap, sino no funciona.*/
//	char (*alleles)[maxSNPCount]=(char (*)[maxSNPCount]) mmap(0,ARRAY_SIZE_alleles, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED,-1, 0); /*fd reemplaza -1*/

        printf("ARRAY_SIZE_alleles=%lu \n",ARRAY_SIZE_alleles);
        printf("ARRAY_SIZE_pos=%lu     \n",ARRAY_SIZE_pos);
        printf("ARRAY_SIZE_chr=%lu     \n",ARRAY_SIZE_chr);


	alleles = (char (*)[acc]) mmap(0,ARRAY_SIZE_alleles, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0); /*fd reemplaza -1*/
	if (alleles == MAP_FAILED)
       	   handle_error("mmap alleles");

	pos     = (int *)         mmap(0,ARRAY_SIZE_pos,     PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (pos == MAP_FAILED)
           handle_error("mmap pos");

	chr     = (int *)         mmap(0,ARRAY_SIZE_chr,     PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (chr == MAP_FAILED)
           handle_error("mmap chr");

  	in   = fopen(inputFileName,"r"); 

 	line = calloc(linealloc,sizeof(char));

  	numberSNPs=0;

       	printf(" InputFile '%s' read started! \n",inputFileName);
        while (fgets(line,linealloc, in)!=NULL){

	//     printf("%s\n",line);
    	   if(strlen(line)>3) {

       	      sscanf(line,"%i %i %s",&chr[numberSNPs],&pos[numberSNPs],alleles[numberSNPs]);

//       	      printf(" inputFile %i :  %i %i %s\n",numberSNPs,chr[numberSNPs],pos[numberSNPs],alleles[numberSNPs]);
       	      numberSNPs++;
           }

       	}
       	printf(" InputFile '%s' read finished! \n",inputFileName);

	fclose(in);

//---------------------------------------------------------Fork: Parallel jobs over the same input file---------------------------------------------------//

	int x = 0;

        pid_t pid,wpid;

	fflush(stdout);
// mehlert	fflush(stdin);
	/*	//Shared memory. Fork makes a copy of all variables for each job, this causes overfloating.
	int shmid = shmget(IPC_PRIVATE, sizeof(char[maxSNPCount][maxSNPCount]),IPC_CREAT | 0666);
	char (*shmPtr)[maxSNPCount] = shmat(shmid, 0, 0);
	fill(shmPtr);*/


        pid=0;

        printf("Starting %u processes\n",index1); fflush(stdout);

#ifdef PARALLEL
	for(x = 0; x < index1 ; x++) {

    	   pid = fork();

   	   if(pid < 0) {
              printf("Error");
              exit(EXIT_FAILURE);

    	   } else if (pid == 0) {

              // ChildProcess
            // printf("Child (%d): %d\n",                            x + 1, getpid());

	      printf("Child (%d): %d 		interval: %i %i  numberSNPs: %i started\n", x + 1, getpid(),interval_start[x],interval_end[x],numberSNPs); fflush(stdout);
	      sprintf(outputName,"%s_%i_%i_AdjMIresults.txt",inputFileName,interval_start[x],interval_end[x]);
		printf("%s outputfile created!\n",outputName);fflush(stdout);

#endif
	      ChildProcess(getpid(), interval_start[x],interval_end[x],alleles,inputFileName,pos,chr,numberSNPs,Threshold);
               
#ifdef PARALLEL
	      printf("Child (%d): %d 		interval: %i %i  finished\n", x + 1, getpid(),interval_start[x],interval_end[x]); fflush(stdout);
              exit(EXIT_SUCCESS); 
  	   }
	}
#endif

 	munmap(alleles,ARRAY_SIZE_alleles);
 	munmap(alleles,acc);
	munmap(chr,ARRAY_SIZE_chr);
	munmap(pos,ARRAY_SIZE_pos);

	free(line);
 	free(buf1);
	free(interval_start);
	free(interval_end);

	// print child status
  	char(* status_print)[200] = (char (*)[200]) malloc(200*index1);

  	int  idx=0;
        int  i;

  	while ( (wpid = wait(&status)) > 0) {
    		sprintf(status_print[idx],"Exit status of %d was %d (%s)\n", (int)wpid, status, (status == EXIT_SUCCESS) ? "EXIT_SUCCESS" : "EXIT_FAILURE");
        	idx++;
  	}

	for ( i=0; i<idx; i++){
       		printf("%s",status_print[i]);
	} 	

   return 0;
}
