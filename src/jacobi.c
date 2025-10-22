////////////////////////////////////////////////////////////////////////////////
// jacobi.c --- TP2 : resolución de un sistema lineal por el método de Jacobi
//
// Auteur          : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Créé le         : Mon Aug 26 14:56:32 2013
// Dern. mod. par  : Jeremie Gaidamour (CNRS/IDRIS) <gaidamou@idris.fr>
// Dern. mod. le   : Mon Aug 26 14:56:32 2013
////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>


// Dimensión por defecto y el tamaño de las matrices
#ifndef VAL_N
#define VAL_N 1201
#endif
#ifndef VAL_D
#define VAL_D 800
#endif

// Inicialización aleatoria
void random_number(double* array, int size) {
  for(int i=0; i<size; i++) {
    // Genera un número entre 0 y 1
    double r = (double)rand() / (double)(RAND_MAX - 1);
    array[i] = r;
  }
}

int main() {

  int n=VAL_N, diag=VAL_D;
  int i, j, iteration=0;
  double a[n][n];
  double x[n], x_courant[n], b[n];
  double norme;

  struct timeval t_elapsed_0, t_elapsed_1;
  double t_elapsed;

  clock_t t_cpu_0, t_cpu_1;
  double  t_cpu;


  // Inicialización de la segunda matriz
  srand(421); // aseguramos reproducibilidad
  random_number(&a[0][0], n*n);
  random_number(&b[0], n);

  // Aumentamos el valor de la diagonal
  for(int i=0; i<n; i++) {
    a[i][i] += diag;
  }

  // Damos una solución inicial
  for(int i=0; i<n; i++) {
    x[i] = 1.0;
  }

  // Iniciamos reloj del cálculo inicial
  t_cpu_0 = clock();

  // Tiempo transcurrido de referencia
  gettimeofday(&t_elapsed_0, NULL);

  // Solución por el método de Jacobi
  while(1) {
    iteration++;

    for(int i=0; i<n; i++) {
      x_courant[i] = 0;
      for(j=0; j<i; j++) {
	x_courant[i] += a[j][i]*x[j];
      }
      for(j=i+1; j<n; j++) {
	x_courant[i] += a[j][i]*x[j];
      }
      x_courant[i] = (b[i] - x_courant[i])/a[i][i];
    }

    // Prueba de convergencia
    {
      double absmax = 0;
      for(int i=0; i<n; i++) {
	double curr = fabs(x[i] - x_courant[i]);
	if (curr > absmax)
	  absmax = curr;
      }
      norme = absmax / n;
    }
    if( (norme <= DBL_EPSILON) || (iteration >= n) ) break;

    // Copiamos de x temporal a x
    memcpy(x, x_courant, n*sizeof(double));
  }

  // Finalización de tiempo transcurrido
  gettimeofday(&t_elapsed_1, NULL);
  t_elapsed = (t_elapsed_1.tv_sec - t_elapsed_0.tv_sec) + (t_elapsed_1.tv_usec - t_elapsed_0.tv_usec) / (double)1000000;

  // Tiempo de CPU de cálculo final
  t_cpu_1 = clock();
  t_cpu = (t_cpu_1 - t_cpu_0) / (double)CLOCKS_PER_SEC;

  // Impresión de resultados
  fprintf(stdout, "\n\n"
	  "   Tamaño de la matriz : %5d\n"
	  "   Iteraciones         : %4d\n"
	  "   Norma               : %10.3E\n"
	  "   Tiempo transcurrido : %10.3E sec.\n"
	  "   Tiempo CPU          : %10.3E sec.\n",
	  n, iteration, norme, t_elapsed, t_cpu
	  );

  return EXIT_SUCCESS;
}
