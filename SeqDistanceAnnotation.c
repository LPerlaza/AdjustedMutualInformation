#include <stdarg.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
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

/*

This annotates ....

Compilation
 gcc -O3 SeqDistanceAnnotation.c -lm -o SeqDistanceAnnotation.exe

*/



//----------------------------------------------------GLobal variables-------------------------------------------------------------//
const unsigned long maxCountResults=1000000;
const unsigned long maxSNPCount=31000000;
const unsigned long acc=1300;
const int           maxCores=40;
const int           ASCII_LIMIT=200;
const unsigned long codelength=4600;



//-----------------------------------------------------------error handeler------------------------------------------------------//
#define handle_error(msg) \
           do { perror(msg); exit(EXIT_FAILURE); } while (0)


//#define PARALLEL 1
           
//-----------------------------------------------------------Functions------------------------------------------------------//
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


float CalculateMean(double data[],int length)
{
 double sum = 0.0;
  int i;
    for(i=0; i<length; ++i){
        sum +=data[i];
    }
    return((float)sum/length);
}


float CalculateMedian(double x[],int n) {
    float temp;
    int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(x[j] < x[i]) {
                // swap elements
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }
    }

    if(n%2==0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((x[n/2] + x[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return x[n/2];
    }
}


float CalculateSD(double data[],int length)
{
    int i;
    float mean,standardDeviation=0;
    mean=CalculateMean(data,length);

    for(i=0; i<length; ++i) 
    {
    standardDeviation += pow(data[i] - mean, 2);
     }
    return sqrt(standardDeviation/(float)length);
}

void  ChildProcess(int pid, int start,int stop,char (*alleles)[acc],double (*accessions)[acc],char (*CodeString)[codelength],char AdjMIResults[300],int *pos, int *pos1,int *pos2, int *chr,int *chr1, int *chr2, int numberSNPs)
{

 FILE    *out;
 char    outputFileName[200];
 char    *Code;
 int     i,k,j,l,y,x,w,z;
 int     *V1,*V2;
 int     *AcCode;


  AcCode=calloc(acc,sizeof(int));
  Code=calloc(acc,sizeof(char));
  V1= calloc(acc,sizeof(int));
  V2= calloc(acc,sizeof(int));

   sprintf(outputFileName,"%s_%i_%i_SeqDistAnnotated.txt",AdjMIResults,start,stop);
   out=fopen(outputFileName,"w");
   fprintf(out,"Chr1\tPos1\tChr2\tPos2\tmeang1\tmeang2\tmediang1\tmediang2\tsdg1\tsdg2\n"); fflush(out);

   for(i=start;i<stop;i++) {
	for (j=0; j<1135; j++) {
			V1[j]=alleles[pos1[i]-1][j];
				//printf("%d ",V1[j]);
				printf("%i %i %i %i %i %i %s \n",k,i,pos1[i],pos[k],chr1[i],chr[k],alleles[k]);
			} 


//        printf("\nCode: %s\n",CodeString[i]);
//despues de esto no lo puedo imprimir porque el strtok modifica la original
	int indx=0;
	Code =strtok(CodeString[i], ",");
	while (Code != NULL) {
               	AcCode[indx++]= atoi(Code);
                //printf("%d ",AcCode[indx]);
                Code = strtok (NULL, ",");
        }
      	int sizeofAcCode=indx;

	int prev = V1[AcCode[0]];
        int count1 =0;
	int count2 =0;

	int group1[1135];
	int group2[1135];

	 for(l=0; l<sizeofAcCode; l++) {
	 //printf("Code split %i %d %i\t",l,AcCode[l],V1[AcCode[l]]);
		if (V1[AcCode[l]] == prev)
		{
		 group1[count1++]=AcCode[l];
		 //printf("letra: %i acc: %i group1: %i  index: %i \n", V1[AcCode[l]],AcCode[l],group1[count1],count1);
		}else {
		group2[count2++]=AcCode[l];
		}
	}
//ojo aqui con hacer bien el conteo si no me aparecen Ns
	
	double subSeqDist1[count1*count1];
	int countA=0;

	for(x=0;x<count1;x++){
	//printf("letra: %i  group1: %i  index: %i \n", V1[group1[i]],group1[i],i);
		for(w=0;w<count1;w++){
		subSeqDist1[countA++]=accessions[group1[x]][group1[w]];
		}
	}

	double subSeqDist2[count2*count2];
	int countB=0;
	for(z=0;z<count2;z++){
	//printf("letra: %i  group2: %i  index: %i \n", V1[group2[i]],group2[i],i);}
		 for(y=0;y<count2;y++){
                subSeqDist2[countB++]=accessions[group1[z]][group1[y]];
                }
        }

	float mean1=CalculateMean(subSeqDist1,countA);
	float median1=CalculateMedian(subSeqDist1,countA);
	float sd1=CalculateSD(subSeqDist1,countA);
	float mean2=CalculateMean(subSeqDist1,countA);
        float median2=CalculateMedian(subSeqDist1,countA);
        float sd2=CalculateSD(subSeqDist1,countA);

	fprintf(out,"%i\t%i\t%i\t%i\t%f\t%f\t%f\t%f\t%f\t%f\n",chr1[i],pos1[i],chr2[i],pos2[i],mean1,median1,sd1,mean2,median2,sd2); fflush(out);
	//mean=sum(subSeqDist2)/countB; Escribir una funcion para esto

  }
}


//-----------------------------------------------------------main------------------------------------------------------//


int main(int argc, char **argv)
{
        FILE    *in,*in1,*in2;
        char    *line,*line1,*line2;
        char    inputFileName[300],SeqDist[300],Intervals[300],AdjMIResults[300],outputName[300];
        int     *chr,*pos,*chr1,*pos1,*chr2,*pos2;
        int     numberSNPs,numberAcc,numberResults;
        int     linealloc=5000;
        int     *interval_start,*interval_end;
        char    (*alleles)[acc];
        double		(*accessions)[acc];
        char	(*CodeString)[codelength];
        int     i,j;

//-----------------------------------------------------------------Reading Arguments--------------------------------------------------------------------//

			if(argc==5) {

      	   sprintf(inputFileName,"%s",argv[1]);
      	   sprintf(SeqDist,"%s", argv[2]);
      	   sprintf(AdjMIResults,"%s", argv[3]);
		   sprintf(Intervals,"%s", argv[4]);
	   	   
	   		printf("\nInput File = %s\n",inputFileName);
	   		printf("\nSeqDist = %s\n",SeqDist);
	   		printf("\nAdjMIResults = %s\n",AdjMIResults);
	   		printf("\nIntervals = %s\n",Intervals);
	  
			} else {
      	   	printf("Wrong number of parameters!\n");
      	   	printf("run with: <program> <inputFile> <SeqDist.Matrix> <AdjMIResults> <Intervals>\n");
      	   	printf("output will be written to file inputFileName_SeqDistAnnotated.txt.\n");
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
        //   printf("Intervals table %i\t%i\n",interval_start[index1],interval_end[index1]);
           index1++;

        }

        fclose(inter_file);
	free(buf1);
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

    //    printf("ARRAY_SIZE_alleles=%lu \n",ARRAY_SIZE_alleles);
     //   printf("ARRAY_SIZE_pos=%lu     \n",ARRAY_SIZE_pos);
      //  printf("ARRAY_SIZE_chr=%lu     \n",ARRAY_SIZE_chr);


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

	//	printf(" inputFile %i :  %i %i %s\n",numberSNPs,chr[numberSNPs],pos[numberSNPs],alleles[numberSNPs]);
       	      numberSNPs++;
           }

       	}
       	printf(" InputFile '%s' read finished! \n",inputFileName);

	fclose(in);
	free(line);
//-----------------------------------------------------------------Reading The Distance Sequence--------------------------------------------------------------------//

//Share memory
	double data;
	int         fd1;
   	 /* Information about the file. */
    	struct      stat s11;
    	int         status1;
//   	size_t      size;
    
    	/* Open the file for reading. */
    	fd1 = open (SeqDist, O_RDONLY);
    	check (fd1 < 0, "open %s failed: %s",SeqDist, strerror (errno));

    	/* Get the size of the file. */
    	status1 = fstat (fd1, & s11);
    	check (status1 < 0, "stat %s failed: %s", SeqDist, strerror (errno));
//    	size = s1.st_size;

	long unsigned int ARRAY_SIZE_acc = (acc*sizeof(double))*(acc);

	/*el numero del segundo array debe ser menor que el ARRAY_SIZE para el mmap, sino no funciona.*/
//	char (*alleles)[maxSNPCount]=(char (*)[maxSNPCount]) mmap(0,ARRAY_SIZE_alleles, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED,-1, 0); /*fd reemplaza -1*/

  //     printf("ARRAY_SIZE_acc=%lu \n",ARRAY_SIZE_acc);
       

	accessions = (double (*)[acc]) mmap(0,ARRAY_SIZE_acc, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0); /*fd reemplaza -1*/
	if (accessions == MAP_FAILED)
       	   handle_error("mmap accessions");

  	in1= fopen(SeqDist,"r"); 
 	//line1 = calloc(linealloc,sizeof(char));
       	printf(" SeqDist '%s' read started! \n",SeqDist);


	 //while (fgets(line1,linealloc, in1)!=NULL){

    	   //	if(strlen(line1)>3) {
				 	for(i = 0; i < 1135; i++)
 					 	{	
      						for(j = 0; j < 1135; j++) 
      						{
      							if(fscanf(in1, "%lf",&data)==1){
							accessions[i][j]=data;}
							else{printf("Failed to read integer.\n");}
							//printf("%i %i %e\n",i,j,accessions[i][j]);
							}
           				}
      			//}
       //}
       	
       	printf("SeqDist '%s' read finished! \n",SeqDist);

	fclose(in1);

//-----------------------------------------------------------------Reading AdjMIResults--------------------------------------------------------------------//

//Share memory

	int         fd2;
   	 /* Information about the file. */
    	struct      stat s111;
    	int         status2;
//   	size_t      size;
    
    	/* Open the file for reading. */
    	fd2 = open (AdjMIResults, O_RDONLY);
    	check (fd2 < 0, "open %s failed: %s",AdjMIResults, strerror (errno));

    	/* Get the size of the file. */
    	status2 = fstat (fd2, & s111);
    	check (status2 < 0, "stat %s failed: %s", AdjMIResults, strerror (errno));
//    	size = s1.st_size;
/// possible longest code is having a vector fro 0:1135 separated by "," 4570 characters: codelength 


	long unsigned int ARRAY_SIZE_Code = maxCountResults*(codelength);
	long unsigned int ARRAY_SIZE_pos1     = maxCountResults*sizeof(int);
	long unsigned int ARRAY_SIZE_chr1     = maxCountResults*sizeof(int);
	long unsigned int ARRAY_SIZE_pos2     = maxCountResults*sizeof(int);
	long unsigned int ARRAY_SIZE_chr2     = maxCountResults*sizeof(int);


	/*el numero del segundo array debe ser menor que el ARRAY_SIZE para el mmap, sino no funciona.*/
//	char (*alleles)[maxSNPCount]=(char (*)[maxSNPCount]) mmap(0,ARRAY_SIZE_alleles, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED,-1, 0); /*fd reemplaza -1*/

//        printf("ARRAY_SIZE_Code=%lu \n",ARRAY_SIZE_Code);
//        printf("ARRAY_SIZE_pos1=%lu     \n",ARRAY_SIZE_pos1);
//        printf("ARRAY_SIZE_chr1=%lu     \n",ARRAY_SIZE_chr1);
//		printf("ARRAY_SIZE_pos2=%lu     \n",ARRAY_SIZE_pos2);
 //       printf("ARRAY_SIZE_chr2=%lu     \n",ARRAY_SIZE_chr2);

	CodeString = (char (*)[codelength]) mmap(0,ARRAY_SIZE_Code, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0); /*fd reemplaza -1*/
	if (CodeString == MAP_FAILED)
       	   handle_error("mmap CodeString");

	pos1     = (int *)         mmap(0,ARRAY_SIZE_pos1, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (pos1 == MAP_FAILED)
           handle_error("mmap pos1");

	chr1     = (int *)         mmap(0,ARRAY_SIZE_chr1,PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (chr1 == MAP_FAILED)
           handle_error("mmap chr1");
           
           
    	pos2     = (int *)         mmap(0,ARRAY_SIZE_pos2, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (pos2 == MAP_FAILED)
           handle_error("mmap pos2");

	chr2     = (int *)         mmap(0,ARRAY_SIZE_chr2, PROT_READ | PROT_WRITE, MAP_ANON | MAP_SHARED, -1, 0);
	if (chr2 == MAP_FAILED)
           handle_error("mmap chr2");
           

  	in2   = fopen(AdjMIResults,"r"); 

 	line2 = calloc(linealloc,sizeof(char));

  	numberResults=0;

       	printf("AdjMIResults '%s' read started! \n",AdjMIResults);
        while (fgets(line2,linealloc, in2)!=NULL){

	//     printf("%s\n",line);
    	   if(strlen(line2)>3) {

       	      sscanf(line2,"%i %i %i %i %s",&chr1[numberResults],&pos1[numberResults],&chr2[numberResults],&pos2[numberResults],CodeString[numberResults]);

 	     // printf(" inputFile %i : %i %i  %i %i %s\n",numberResults,chr1[numberResults],pos1[numberResults],chr2[numberResults],pos2[numberResults],CodeString[numberResults]);
       	      numberResults++;
           }

       	}
       	printf(" AdjMIResults '%s' read finished! \n",AdjMIResults);

	fclose(in2);
	free(line1);

//---------------------------------------------------------Fork: Parallel jobs over the same input file---------------------------------------------------//

	int x = 0;

        pid_t pid,wpid;
        pid=0;


#ifdef PARALLEL
	for(x = 0; x < index1 ; x++) {
			pid = fork();
				if(pid < 0) {
              	printf("Error");
              	exit(EXIT_FAILURE);
    	   		} else if (pid == 0) {

	      		printf("Child (%d): %d 		interval: %i %i  numberSNPs: %i started\n", x + 1, getpid(),interval_start[x],interval_end[x],numberSNPs); fflush(stdout);

				#endif

	      		ChildProcess(getpid(), interval_start[x],interval_end[x],alleles,accessions,CodeString,AdjMIResults,pos,pos1,pos2,chr,chr1,chr2,numberSNPs);
               		printf("%s outputfile created!\n",outputName);fflush(stdout);
		
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
	munmap(chr,ARRAY_SIZE_chr1);
	munmap(pos,ARRAY_SIZE_pos1);
	munmap(chr,ARRAY_SIZE_chr2);
	munmap(pos,ARRAY_SIZE_pos2);
	munmap(accessions,ARRAY_SIZE_acc);
	munmap(accessions,acc);
	
	free(interval_start);
	free(interval_end);

	// print child status
  	char(* status_print)[200] = (char (*)[200]) malloc(200*index1);

  	int  idx=0;
       
  	while ( (wpid = wait(&status)) > 0) {
    		sprintf(status_print[idx],"Exit status of %d was %d (%s)\n", (int)wpid, status, (status == EXIT_SUCCESS) ? "EXIT_SUCCESS" : "EXIT_FAILURE");
        	idx++;
  	}

	for ( i=0; i<idx; i++){
       		printf("%s",status_print[i]);
	} 	


/**/

   return 0;
}






