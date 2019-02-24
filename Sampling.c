#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>
#include <assert.h>


/* compile with

gcc -O3 Sampling.c -o Sampling.exe

Dist    Freq
10000   13746129

chr     pos     snp     bin
1       73      73      100033


*/


main(int argc, char **argv)
{
 int maxSNPCount=2000000;
 int maxdistCount=12000;
 int maxCount=42000000;
 int binwidth;

  FILE *out;
  char line[1000],outputFileName[200],inputFileName[200],Distances[200];
  int  *SNP1,*Chr,*Pos,*bin,i,k;
  unsigned long int *Dist,*Count;


  if(argc==3)
    {
      sprintf(inputFileName,"%s",argv[1]);
      sprintf(Distances,"%s",argv[2]);
      sprintf(outputFileName,"%s_%s_Sampled.txt",inputFileName,Distances);
    }
  else
    {
      printf("Wrong number of parameters!\n");
      printf("run with: <program> <SNPsunique> <DistancesProb>\n");
      printf("output will be written to file inputFile_Sampled.txt\n");
      exit(0);
    }

	out=fopen(outputFileName,"w");
	FILE *in;
        char *buf1;
        int  index1=0;

        buf1= calloc(1500,sizeof(char));

        in     = fopen(inputFileName,"r");
	SNP1=calloc(maxSNPCount,sizeof(int));
  	bin=calloc(maxSNPCount,sizeof(int));
  	Pos=calloc(maxSNPCount,sizeof(int)); 
        Chr=calloc(maxSNPCount,sizeof(int));

        while (fgets(buf1,1500,in)!=NULL) {

           sscanf(buf1,"%i %i %i %i",&Chr[index1],&Pos[index1],&SNP1[index1],&bin[index1]);
           //printf("SNP:%i\tChr:%i\tBin:%i\tIndex:%i\n",SNP1[index1],Chr[index1],bin[index1],index1);
           index1++;

        }

        fclose(in);



        FILE *in2;
        char *buf2;
        int   index2=0;

        buf2= calloc(1500,sizeof(char));

        in2     = fopen(Distances,"r");
        Dist=calloc(maxdistCount,sizeof(unsigned long int));
        Count=calloc(maxdistCount,sizeof(unsigned long int));

        while (fgets(buf2,1500,in2)!=NULL) {

           sscanf(buf2,"%lu %lu",&Dist[index2],&Count[index2]);
           //printf("Dist:%lu\tCount:%lu\tIndex:%i\n",Dist[index2],Count[index2],index2);
           index2++;

        }

        fclose(in2);

	binwidth=Dist[0]-1;
//	printf("\n%i\n",binwidth); //

	int randSNP1;
	int randSNP2,sucess=0;
	unsigned long int randDistance;
	int *randomSNPs1,*randomSNPs2,*randomChr1,*randomChr2,*randombins1,*randombins2,*randomPos1,*randomPos2;

       randomSNPs1=calloc(maxCount,sizeof(int));
       randomSNPs2=calloc(maxCount,sizeof(int));
       randomChr1=calloc(maxCount,sizeof(int));
       randomChr2=calloc(maxCount,sizeof(int));
       randombins1=calloc(maxCount,sizeof(int));
       randombins2=calloc(maxCount,sizeof(int));
       randomPos1=calloc(maxCount,sizeof(int));
       randomPos2=calloc(maxCount,sizeof(int));


	int index3=0;
	unsigned long int index4;

	for(i=0;i<index2;i++){
		index4=0;

		if(Count[i]==0)continue;
//		 printf("1index4=%lu i:%i dist:%lu count:%lu\n",index4,i,Dist[i],Count[i]);
		while(index4 < Count[i]) {
		//	printf("while with %i Dist:%lu\n",i,Dist[i]);
    			randSNP1 = rand()% index1;
			randSNP2 = rand()% index1;
			randDistance=abs(SNP1[randSNP1]-SNP1[randSNP2]);
			int start=abs(Dist[i]-binwidth); 
			int end =Dist[i];
//		printf(" start: %i end:%i\n",start,end);

				if( randDistance > start && randDistance < end){
				 index3++;
				 index4++;
				 randomSNPs1[index3]=SNP1[randSNP1];
				 randomSNPs2[index3]=SNP1[randSNP2];
				 randomChr1[index3]=Chr[randSNP1];
				 randomChr2[index3]=Chr[randSNP2];
				 randombins1[index3]=bin[randSNP1];
				 randombins2[index3]=bin[randSNP2];
				 randomPos1[index3]=Pos[randSNP1];
				 randomPos2[index3]=Pos[randSNP2];

		                  }
		}



	}



		int x;

        for (x=0;x<index3;x++)
        {
     fprintf(out,"%i\t%i\t%i\t%i\t%i\t%i\t%i\t%i\n",randomChr1[x],randomPos1[x],randomSNPs1[x],randomChr2[x],randomPos2[x],randomSNPs2[x],randombins1[x],randombins2[x]);
         }

  fclose(out);

	free(SNP1);
	free(Pos);
	free(Chr);
	free(Dist);
	free(Count);
        free(randomSNPs1);
        free(randomSNPs2);
        free(randomChr1);
        free(randomChr2);
	free(randomPos1);
        free(randomPos2);
        free(randombins1);
        free(randombins2);
	free(bin);
	free(buf1);
	free(buf2);

 }
