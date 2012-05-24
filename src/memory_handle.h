#include "R.h"
#ifndef MEMORY_HANDLE__
#define MEMORY_HANDLE__
/*this can be used to handle the memory allocation independent of R by chaning
 myalloc back to malloc, etc*/

/*Since R has used the marco, there is no neat way for me to get around it, 
  though macro should be avoided as much as possible*/

#define mymalloc(n,t) Calloc(n,t) 
#define mycalloc(n,t) Calloc(n,t)
#define myfree(p) Free(p)

//If we want to have R indepedent code,
/* The above needs to be modified as 
#define mymalloc(n,t) (t*) malloc( (size_t)(n)*sizeof(t)), 
#deifne mycalloc(n,t) (t*) calloc( (size_t)(n), sizeof(t))
#define myfree(p) free(p);*/
#endif
