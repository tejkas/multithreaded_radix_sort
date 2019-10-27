/********************************************************
 * Kernels to be optimized for the CS:APP Performance Lab
 ********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include <pthread.h>
#include <semaphore.h>
#include <algorithm>
#include <string.h>
#include "helper.h"

pthread_t threads[32];

/* 
 * Please fill in the following team struct 
 */
user_t user = {
    (char*)  "11111111",            /* UID */
    (char*)  "abc",          /* Full name */
    (char*)  "abc@acd.com",  /* Email */
};

//If you worked with anyone in this project, please indicate that here:
// I worked with  "Myles Byrne, Bryan Le ".

// Of course, your entire file should be hand-written by you.  You are free to
// look at tutorials and academic literature for radix-sort based sorting.  You
// are not allowed to look at or copy from implementations on online code repositories
// like github.


//  You should modify only this file, but feel free to modify however you like!


/*
 * setup - This function runs one time, and will not be timed.
 *         You can do whatever initialization you need here, but
 *         it is not required -- don't use if you don't want to.  ^_^
 */
void setup() {
  //..

  // So, in my experiments, it take fewer cycles for each run if we 
  // waste some time doing absolute nothing in particular
  // at the begining of the program.  It might be because of Intel's Turbo
  // mode and DVFS somehow??  TBH IDK, but I would leave this in if I were
  // you.
  for(int i = 0; i < 50000;++i) {
     do_nothing(i);
  }
}

/******************************************************
 * Your different versions of the singlethread kernel go here
 ******************************************************/
bool kvp_compare(kvp lhs, kvp rhs) { 
  return lhs.key < rhs.key; 
}

/*
 * naive_singlethread - The naive baseline version of singlethread 
 */
char naive_singlethread_descr[] = "naive_singlethread: Naive baseline implementation";
void naive_singlethread(int dim, kvp *src, kvp *dst) 
{
    //This is the built-in stable sort if you want to try it
    //memcpy(dst, src, dim*sizeof(kvp));
    //std::stable_sort(dst, dst+dim, kvp_compare);
    //return;

    int log_radix=8; //Radix of radix-sort is 2^8
    int iters=(sizeof(unsigned int)*8/log_radix);

    // 256 buckets for 2^8 bins, one for each iteration 
    unsigned long long buckets[256+1][iters];
    unsigned long long sum[256+1][iters];

    for(int iter = 0; iter < iters; ++iter) {
      for(int i = 0; i < bucket_size(log_radix); ++i) {
        buckets[i][iter]=0;
        sum[i][iter]=0;
      }
 
      
      //1. Generate the bucket count
      for(int i = 0; i < dim; ++i) {
        int index = gen_shift(src[i].key,iter*log_radix,
                              (bucket_size(log_radix)-1))+1;
        buckets[index][iter]++;
      }

      //2. Perform scan
      for(int i = 1; i < bucket_size(log_radix); ++i) {
        sum[i][iter] = buckets[i][iter] + sum[i-1][iter];
      }

      //3. Move Data items
      for(int i = 0; i < dim; ++i) {
        int index = gen_shift(src[i].key,iter*log_radix,
                              bucket_size(log_radix)-1);
        int out_index = sum[index][iter];
	 move_kvp(dst,src,i,out_index);

	//	dst[out_index].key = src[i].key;
	//	dst[out_index].value = src[i].value;
        sum[index][iter]++;
      }

      // Move dest back to source
      for(int i = 0; i < dim; ++i) {
	 move_kvp(src,dst,i,i);

	//	dst[i].key = src[i].key;
	//	dst[i].value = src[i].value;
      }

    }
}


/*
 * singlethread - Your current working version of singlethread. 
 * IMPORTANT: This is the version you will be graded on
 */


char singlethread_descr[] = "singlethread: Current working version";
void singlethread(int dim, kvp *src, kvp *dst) 
{
   

   //NEW PART
     kvp *temp;
    int log_radix=8; //Radix of radix-sort is 2^8
    int iters=(sizeof(unsigned int)*8/log_radix);

    // 256 buckets for 2^8 bins, one for each iteration 
    // unsigned long long buckets[256+1][iters];
    // unsigned long long sum[256+1][iters];

    unsigned long long buckets[iters][256+1];
    unsigned long long sum[iters][256+1];

    

    size_t bucketsLength = sizeof(buckets[0][0])*257*iters;
    size_t sumLength = sizeof(sum[0][0])*257*iters;

    int bucketSize = bucket_size(log_radix);

    //experimental implementation here
    
    // memset(buckets, 0, bucketsLength);
    //  memset(sum, 0, sumLength);

      



    
    // for(int iter = 0; iter < iters; ++iter) {
      // for(int i = 0; i < bucketSize; ++i) {
      // buckets[i][iter]=0;
      // sum[i][iter]=0;
      // }

      memset(buckets, 0, bucketsLength);
      memset(sum, 0, sumLength);
      
      
      
      int shift0 = 0*log_radix;
      int shift1 = 1*log_radix;
      int shift2 = 2*log_radix;
      int shift3 = 3*log_radix;
      int mask  = (bucketSize-1);
      //1. Generate the bucket count
      for(int i = 0; i < dim; ++i) {
	// int index = gen_shift(src[i].key,shift,
                              //(bucketSize-1))+1;
	int index0 = ((src[i].key >> shift0) & mask) + 1;
        buckets[0][index0]++;
	int index1 = ((src[i].key >> shift1) & mask) + 1;
        buckets[1][index1]++;
	int index2 = ((src[i].key >> shift2) & mask) + 1;
        buckets[2][index2]++;
	int index3 = ((src[i].key >> shift3) & mask) + 1;
        buckets[3][index3]++;
	
      }

      //2. Perform scan
      for(int i = 1; i < bucketSize; ++i) {
        sum[0][i] = buckets[0][i] + sum[0][i-1];
	sum[1][i] = buckets[1][i] + sum[1][i-1];
	sum[2][i] = buckets[2][i] + sum[2][i-1];
	sum[3][i] = buckets[3][i] + sum[3][i-1];
      }

   for(int iter = 0; iter < iters; ++iter){

     int shift = iter*log_radix;
	
      //3. Move Data items
      for(int i = 0; i < dim; i+=8) {
	// int index = gen_shift(src[i].key,iter*log_radix,
	//	      bucke


	//new

	int index0 = ((src[i].key >> shift) & mask);
        int out_index0 = sum[iter][index0];

	//		int index1 = ((src[i+1].key >> shift) & mask);
	//	int out_index1 = sum[iter][index1];
	/*
	int index2 = ((src[i+2].key >> shift) & mask);
	int out_index2 = sum[iter][index2];

	int index3 = ((src[i+3].key >> shift) & mask);
	int out_index3 = sum[iter][index3];
	*/
	//    move_kvp(dst,src,i,out_index);
	//    move_kvp(dst,src,i,out_index);
	
	//    move_kvp(dst,src,i,out_index);
	dst[out_index0] = src[i];
		//	dst[out_index0]  = src[i];
        sum[iter][index0]++;

	int index1 = ((src[i+1].key >> shift) & mask);
	int out_index1 = sum[iter][index1];

	dst[out_index1] = src[i+1];
	sum[iter][index1]++;

	int index2 = ((src[i+2].key >> shift) & mask);
	int out_index2 = sum[iter][index2];

	dst[out_index2] = src[i+2];
	sum[iter][index2]++;

	int index3 = ((src[i+3].key >> shift) & mask);
	int out_index3 = sum[iter][index3];

	dst[out_index3] = src[i+3];
	sum[iter][index3]++;

	int index4 = ((src[i+4].key >> shift) & mask);
	int out_index4 = sum[iter][index4];

	dst[out_index4] = src[i+4];
	sum[iter][index4]++;

	int index5 = ((src[i+5].key >> shift) & mask);
	int out_index5 = sum[iter][index5];

	dst[out_index5] = src[i+5];
	sum[iter][index5]++;

	int index6 = ((src[i+6].key >> shift) & mask);
	int out_index6 = sum[iter][index6];

	dst[out_index6] = src[i+6];
	sum[iter][index6]++;

	int index7 = ((src[i+7].key >> shift) & mask);
	int out_index7 = sum[iter][index7];

	dst[out_index7] = src[i+7];
	sum[iter][index7]++;

		      
	
	//	dst[out_index1] = src[i+1];
	//	dst[out_index1]  = src[i+1];
	// sum[iter][index1]++;
	/*
	dst[out_index2] = src[i+2];
	//	dst[out_index2]  = src[i+2];
        sum[iter][index2]++;

	dst[out_index3] = src[i+3];
	//	dst[out_index3]  = src[i+3];
        sum[iter][index3]++;
	*/

	
      }


           if(iter != iters-1)
      	{
	  
	  temp  = dst;
	  dst = src;
	  src = temp;
	  
	 
	  
	  	}
	   else
	  {
	  
	

	    //       Move dest back to source
		         for(int i = 0; i < dim; i+=8) {
	// move_kvp(src,dst,i,i);

			   src[i] = dst[i];
			   src[i+1] = dst[i+1];
			   src[i+2] = dst[i+2];
			   src[i+3] = dst[i+3];
			   src[i+4] = dst[i+4];
			   src[i+5] = dst[i+5];
			   src[i+6] = dst[i+6];
			   src[i+7] = dst[i+7];    


	 // (*src).key = (*dst).key;
	 // (*src).value = (*dst).value;
	 // src++;
	 // dst++;
		  	}

	 	}
      
      

    }



   
}


char singlethread_testdescr2[] = "singlethread: new test 2 ";
void test2_singlethread(int dim, kvp *src, kvp *dst) 
{
  kvp *temp;
    int log_radix=8; //Radix of radix-sort is 2^8
    int iters=(sizeof(unsigned int)*8/log_radix);

    // 256 buckets for 2^8 bins, one for each iteration 
    // unsigned long long buckets[256+1][iters];
    // unsigned long long sum[256+1][iters];

    unsigned long long buckets[iters][256+1];
    unsigned long long sum[iters][256+1];

    

    size_t bucketsLength = sizeof(buckets[0][0])*257*iters;
    size_t sumLength = sizeof(sum[0][0])*257*iters;

    int bucketSize = bucket_size(log_radix);

    //experimental implementation here
    
    // memset(buckets, 0, bucketsLength);
    //  memset(sum, 0, sumLength);

      



    
    // for(int iter = 0; iter < iters; ++iter) {
      // for(int i = 0; i < bucketSize; ++i) {
      // buckets[i][iter]=0;
      // sum[i][iter]=0;
      // }

      memset(buckets, 0, bucketsLength);
      memset(sum, 0, sumLength);
      
      
      
      int shift0 = 0*log_radix;
      int shift1 = 1*log_radix;
      int shift2 = 2*log_radix;
      int shift3 = 3*log_radix;
      int mask  = (bucketSize-1);
      //1. Generate the bucket count
      for(int i = 0; i < dim; ++i) {
	// int index = gen_shift(src[i].key,shift,
                              //(bucketSize-1))+1;
	int index0 = ((src[i].key >> shift0) & mask) + 1;
        buckets[0][index0]++;
	int index1 = ((src[i].key >> shift1) & mask) + 1;
        buckets[1][index1]++;
	int index2 = ((src[i].key >> shift2) & mask) + 1;
        buckets[2][index2]++;
	int index3 = ((src[i].key >> shift3) & mask) + 1;
        buckets[3][index3]++;
	
      }

      //2. Perform scan
      for(int i = 1; i < bucketSize; ++i) {
        sum[0][i] = buckets[0][i] + sum[0][i-1];
	sum[1][i] = buckets[1][i] + sum[1][i-1];
	sum[2][i] = buckets[2][i] + sum[2][i-1];
	sum[3][i] = buckets[3][i] + sum[3][i-1];
      }

   for(int iter = 0; iter < iters; ++iter){

     int shift = iter*log_radix;
	
      //3. Move Data items
      for(int i = 0; i < dim; i++) {
	// int index = gen_shift(src[i].key,iter*log_radix,
	//	      bucke


	//new

	int index0 = ((src[i].key >> shift) & mask);
        int out_index0 = sum[iter][index0];

	/*	int index1 = ((src[i+1].key >> shift) & mask);
	int out_index1 = sum[iter][index1];

	int index2 = ((src[i+2].key >> shift) & mask);
	int out_index2 = sum[iter][index2];

	int index3 = ((src[i+3].key >> shift) & mask);
	int out_index3 = sum[iter][index3];
	*/
	//    move_kvp(dst,src,i,out_index);
	//    move_kvp(dst,src,i,out_index);
	
	//    move_kvp(dst,src,i,out_index);
	dst[out_index0] = src[i];
		//	dst[out_index0]  = src[i];
        sum[iter][index0]++;
	/*
	dst[out_index1] = src[i+1];
	//	dst[out_index1]  = src[i+1];
        sum[iter][index1]++;

	dst[out_index2] = src[i+2];
	//	dst[out_index2]  = src[i+2];
        sum[iter][index2]++;

	dst[out_index3] = src[i+3];
	//	dst[out_index3]  = src[i+3];
        sum[iter][index3]++;
	*/

	
      }


           if(iter != iters-1)
      	{
	  
	  temp  = dst;
	  dst = src;
	  src = temp;
	  
	 
	  
	  	}
	   else
	  {
	  
	

	    //       Move dest back to source
		         for(int i = 0; i < dim; ++i) {
	// move_kvp(src,dst,i,i);

			   src[i] = dst[i];
			   


	 // (*src).key = (*dst).key;
	 // (*src).value = (*dst).value;
	 // src++;
	 // dst++;
		  	}

	 	}
      
      

    }
 
}




/********************************************************************* 
 * register_singlethread_functions - Register all of your different versions
 *     of the singlethread kernel with the driver by calling the
 *     add_singlethread_function() for each test function.  When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_singlethread_functions() {
    add_singlethread_function(&naive_singlethread, naive_singlethread_descr);
    add_singlethread_function(&singlethread, singlethread_descr);
    /* ... Register additional test functions here */
    //  add_singlethread_function(&test_singlethread, singlethread_testdescr);
    add_singlethread_function(&test2_singlethread, singlethread_testdescr2);
}



// ----------------------- do multi-threaded versions here ------------------


// I'm kind of cheating to pass global variables to my thread function
// There are nicer ways to do this, but w/e
int gdim;
kvp* gsrc;
kvp* gdst;

void *do_sort(void* threadid) {
   long tid = (long)threadid;
   //printf("Hello from thread %d\n", tid);
   int dim = gdim;
   kvp* src = gsrc;
   kvp* dst = gdst;
  
   naive_singlethread(dim,src,dst);   
   return 0;
}


/* 
 * naive_multithread - The naive baseline version of multithread 
 */
char naive_multithread_descr[] = "naive_multithread: Naive baseline implementation";
void naive_multithread(int dim, kvp *src, kvp *dst) 
{
    gdim=dim;
    gsrc=src;
    gdst=dst;

    //call one thread to do our work -- I'm sure that will make things go faster
    pthread_create(&threads[0], NULL, do_sort, (void *)0 /*tid*/);

    void** ret=0;
    pthread_join(threads[0],ret);
}


/* 
 * multithread - Your current working version of multithread
 * IMPORTANT: This is the version you will be graded on
 */

#define NUMTHREADS 16

//global iters
int giters;

//global array
unsigned long long tag[NUMTHREADS+1][256+1];

//function declaration
void *multithreaded_algorithm(void * threadNum);


//barrier declaration
pthread_barrier_t mybarrier;


//global radix
int gradix;

char multithread_descr[] = "multithread: Current working version";
void multithread(int dim, kvp *src, kvp *dst) 
{

  if(dim <= 65536)
    {
      singlethread(dim, src, dst);
      return;
    }
  // fprintf(stdout, "STARTED THE FUNCTION! \n");
  int log_radix = 8;
  gradix = log_radix;
  gdim = dim;
  gsrc = src;
  gdst = dst;
  int iters = (sizeof(unsigned int)*8/log_radix);
  giters = iters;
  int thread_id[NUMTHREADS];

  pthread_barrier_init(&mybarrier,NULL,NUMTHREADS);

  //create threads here for step 1
  for(int i = 0; i < NUMTHREADS; i++)
    {
      thread_id[i] = i;
      // fprintf(stdout, "gonna create some threads!");
      int ret = pthread_create(&threads[i], NULL, multithreaded_algorithm, &thread_id[i]);
      
      if(ret)
	     {
	      fprintf(stderr, "Thread Creation Error");
    	  exit(-1);
       }

      //  fprintf(stdout, "thread %d  created! \n", i);
    }

  //Free the threads!
    for(int i = 0; i < NUMTHREADS; i++)
    {
      int ret = pthread_join(threads[i], NULL);
      if(ret)
      {
        fprintf(stderr, "Thread freeing error");
        exit(-1);
      }
    }

    pthread_barrier_destroy(&mybarrier);
  
			       
    
}

void *multithreaded_algorithm(void * threadNum)
{

  // fprintf(stdout, "algorithm entered! \n");
  int tagSize = sizeof(tag[0][0])*(NUMTHREADS+1)*(257);
  int mask = bucket_size(gradix)-1; //255 for 8
  int threadID = *(int*)threadNum;
  int tstart =(threadID)*(gdim/NUMTHREADS);
  int tmax = (threadID + 1)*(gdim/NUMTHREADS);

   
   
  for(int iter = 0; iter < giters; iter++)
    {
      // fprintf(stdout, "ITER LOOP ENTERED! \n");
      //zero out arrays
      memset(tag, 0, tagSize);

      pthread_barrier_wait(&mybarrier);

      // fprintf(stdout, "Past the first barrier! \n");
      //Histo phase
      int shift = iter*gradix;

      // fprintf(stdout, "tstart for this thread %d  is: %d \n",threadID, tstart);
      // fprintf(stdout, "tmax for this thread %d is: %d \n", threadID, tmax);
      //  int mask = bucket_size(gradix)-1; //255 for 8
      // int threadID = *(int*)threadNum;
      // int t =(threadID)*(gdim/4);
       //   int tmax = (threadID + 1)*(0.25)*gdim;
      for(int t = tstart; t < tmax; t++)
	   {
	     // fprintf(stdout, "histo creation beginning \n");
	     int index = ((gsrc[t].key >> shift) & mask);
	     tag[threadID+1][index]++;
	     // fprintf(stdout, "%X \n", tag[threadID+1][index]);
	   }

      pthread_barrier_wait(&mybarrier);

      
      //prefix sum scan (only one thread allowed, it is singlethreaded)
            if(threadID == 1)
        { 
          
          for(int n = 0; n < 257; n++)
  	     {
  	       for(int m = 1; m < NUMTHREADS+1; m++)
  	        {
  	      
  	         tag[m][n] = tag[m][n] + tag[m-1][n];
  	         if(m == NUMTHREADS && n!= 256)
  		        {
  		          tag[0][n+1] = tag[m][n];
  		        }
  	  
  	        }
          }
        }   


      pthread_barrier_wait(&mybarrier);

      //distribute the data in parallel
      for(int t = tstart; t < tmax; t+=4)
      {
        int index = ((gsrc[t].key >> shift) & mask);
        int offset = tag[threadID][index];

       	gdst[offset] = gsrc[t];
	
        tag[threadID][index]++;

        int index1 = ((gsrc[t+1].key >> shift) & mask);
        int offset1 = tag[threadID][index1];

       	gdst[offset1] = gsrc[t+1];
	
        tag[threadID][index1]++;

	 int index2 = ((gsrc[t+2].key >> shift) & mask);
        int offset2 = tag[threadID][index2];

       	gdst[offset2] = gsrc[t+2];
	
        tag[threadID][index2]++;

	 int index3 = ((gsrc[t+3].key >> shift) & mask);
        int offset3 = tag[threadID][index3];

       	gdst[offset3] = gsrc[t+3];
	
        tag[threadID][index3]++;

	
      }

      pthread_barrier_wait(&mybarrier);

      //copy dst to src
      for(int t = tstart; t < tmax; t+=8)
      {
        gsrc[t] = gdst[t];
	gsrc[t+1] = gdst[t+1];
	gsrc[t+2] = gdst[t+2];
	gsrc[t+3] = gdst[t+3];
        gsrc[t+4] = gdst[t+4];
	gsrc[t+5] = gdst[t+5];
	gsrc[t+6] = gdst[t+6];
	gsrc[t+7] = gdst[t+7];
      
        
      }

      pthread_barrier_wait(&mybarrier);

      
     }
      
  return NULL;
      
}

unsigned long long bag[257][NUMTHREADS+1];

void* multithreadedalgo(void * threadNum);

char multithread_descr22[] = "multithread: Current working version 222222";
void multithread2(int dim, kvp *src, kvp *dst) 
{

  if(dim <= 65536)
    {
      singlethread(dim, src, dst);
      return;
    }
  // fprintf(stdout, "STARTED THE FUNCTION! \n");
  int log_radix = 8;
  gradix = log_radix;
  gdim = dim;
  gsrc = src;
  gdst = dst;
  int iters = (sizeof(unsigned int)*8/log_radix);
  giters = iters;
  int thread_id[NUMTHREADS];

  pthread_barrier_init(&mybarrier,NULL,NUMTHREADS);

  //create threads here for step 1
  for(int i = 0; i < NUMTHREADS; i++)
    {
      thread_id[i] = i;
      // fprintf(stdout, "gonna create some threads!");
      int ret = pthread_create(&threads[i], NULL, multithreadedalgo, &thread_id[i]);
      
      if(ret)
	     {
	      fprintf(stderr, "Thread Creation Error");
    	  exit(-1);
       }

      //  fprintf(stdout, "thread %d  created! \n", i);
    }

  //Free the threads!
    for(int i = 0; i < NUMTHREADS; i++)
    {
      int ret = pthread_join(threads[i], NULL);
      if(ret)
      {
        fprintf(stderr, "Thread freeing error");
        exit(-1);
      }
    }

    pthread_barrier_destroy(&mybarrier);
  
			       
    
}

void *multithreadedalgo(void * threadNum)
{

  // fprintf(stdout, "algorithm entered! \n");
  int tagSize = sizeof(tag[0][0])*(NUMTHREADS+1)*(257);
  int mask = bucket_size(gradix)-1; //255 for 8
  int threadID = *(int*)threadNum;
  int tstart =(threadID)*(gdim/NUMTHREADS);
  int tmax = (threadID + 1)*(gdim/NUMTHREADS);

   
   
  for(int iter = 0; iter < giters; iter++)
    {
      // fprintf(stdout, "ITER LOOP ENTERED! \n");
      //zero out arrays
      memset(bag, 0, tagSize);

      pthread_barrier_wait(&mybarrier);

      // fprintf(stdout, "Past the first barrier! \n");
      //Histo phase
      int shift = iter*gradix;

      // fprintf(stdout, "tstart for this thread %d  is: %d \n",threadID, tstart);
      // fprintf(stdout, "tmax for this thread %d is: %d \n", threadID, tmax);
      //  int mask = bucket_size(gradix)-1; //255 for 8
      // int threadID = *(int*)threadNum;
      // int t =(threadID)*(gdim/4);
       //   int tmax = (threadID + 1)*(0.25)*gdim;
      for(int t = tstart; t < tmax; t++)
	   {
	     // fprintf(stdout, "histo creation beginning \n");
	     int index = ((gsrc[t].key >> shift) & mask);
	     bag[index][threadID+1]++;
	     // fprintf(stdout, "%X \n", tag[threadID+1][index]);
	   }

      pthread_barrier_wait(&mybarrier);

      
      //prefix sum scan (only one thread allowed, it is singlethreaded)
            if(threadID == 1)
        { 
          
          for(int n = 0; n < 257; n++)
  	     {
  	       for(int m = 1; m < NUMTHREADS+1; m++)
  	        {
  	      
  	         bag[n][m] = bag[n][m] + bag[n][m-1];
  	         if(m == NUMTHREADS && n!= 256)
  		        {
  		          bag[n+1][0] = bag[n][m];
  		        }
  	  
  	        }
          }
        }   


      pthread_barrier_wait(&mybarrier);

      //distribute the data in parallel
      for(int t = tstart; t < tmax; t+=4)
      {
        int index = ((gsrc[t].key >> shift) & mask);
        int offset = bag[index][threadID];

       	gdst[offset] = gsrc[t];
	
        bag[index][threadID]++;

        int index1 = ((gsrc[t+1].key >> shift) & mask);
        int offset1 = bag[index1][threadID];

       	gdst[offset1] = gsrc[t+1];
	
        bag[index1][threadID]++;

	 int index2 = ((gsrc[t+2].key >> shift) & mask);
        int offset2 = bag[index2][threadID];

       	gdst[offset2] = gsrc[t+2];
	
        bag[index2][threadID]++;

	 int index3 = ((gsrc[t+3].key >> shift) & mask);
        int offset3 = bag[index3][threadID];

       	gdst[offset3] = gsrc[t+3];
	
        bag[index3][threadID]++;

	
      }

      pthread_barrier_wait(&mybarrier);

      //copy dst to src
      for(int t = tstart; t < tmax; t+=4)
      {
        gsrc[t] = gdst[t];
	gsrc[t+1] = gdst[t+1];
	gsrc[t+2] = gdst[t+2];
	gsrc[t+3] = gdst[t+3];
        
      }

      pthread_barrier_wait(&mybarrier);

      
     }
      
  return NULL;
      
}

/*********************************************************************
 * register_multithread_functions - Register all of your different versions
 *     of the multithread kernel with the driver by calling the
 *     add_multithread_function() for each test function. When you run the
 *     driver program, it will test and report the performance of each
 *     registered test function.  
 *********************************************************************/

void register_multithread_functions() 
{
    add_multithread_function(&naive_multithread, naive_multithread_descr);   
    add_multithread_function(&multithread, multithread_descr);   
    /* ... Register additional test functions here */
    //    add_multithread_function(&multithread2, multithread_descr22);
}


