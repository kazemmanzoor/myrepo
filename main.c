#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Added for seconds() 
#include <time.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>


// Includes header of my library
#include "mylib.h"

#define SIZE 32

double seconds()
/* Return elapsed time on second from January 1, 1970 */
{
    struct timeval tmp;
    double sec;
    gettimeofday( &tmp, (struct timezone *)0 );
    sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;
    return sec;
}


void print_matrix( double *A ){

  int i, j;

  fprintf( stdout, "\n" );
  
  for( i = 0; i < SIZE; i++ ) {
    for( j = 0; j < SIZE; j++ ) {

      fprintf( stdout, "%.3g \t", A[ i * SIZE + j ] );

    }
    fprintf( stdout, "\n" );
  }

}


int main ( int * argc, char * argv[] ){

  double *A, *B, *BT, *C1, *C2;
  int i, j, k;
  double t1, t2;
  
  A = (double *) malloc( SIZE * SIZE * sizeof( double ) );
  B = (double *) malloc( SIZE * SIZE * sizeof( double ) );
  BT = (double *) malloc( SIZE * SIZE * sizeof( double ) );

  // Matrix C1 stores A * B (MM Multiplication)
  C1 = (double *) malloc( SIZE * SIZE * sizeof( double ) );

  // Matrix C2 stores A * Transpose(B) (MM Multiplication)
  C2 = (double *) malloc( SIZE * SIZE * sizeof( double ) );
  
  for( i = 0; i < SIZE; i++ ) {
    for( j = 0; j < SIZE; j++ ) {
      A[ i * SIZE + j ] = (double) ( rand() % 100 + 1);
      B[ i * SIZE + j ] = (double) ( rand() % 100 + 1);
      C1[ i * SIZE + j ] = 0.0;
      C2[ i * SIZE + j ] = 0.0;
    }
  }

  t1 = seconds();
  matmul( A, B, C1, SIZE );

#ifdef __DEBUG
  if( SIZE < 20 ) print_matrix( B );
#endif
  
  transpose( B, SIZE );

#ifdef __DEBUG
  if( SIZE < 20 ) print_matrix( B );
#endif
  
  matmul( A, B, C2, SIZE );

  t2 = seconds();
  
  fprintf( stdout, "\n All operations where correctly performed " );
  fprintf( stdout, "\n Time to solution (seconds): %.3g \n\n", t2 - t1 );

  return 0;

}

