#include <stdio.h>
#include <unistd.h> //para sleep()
#include <time.h>   //para el random
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

/**
 * @brief Archivo en C con la solución a la tarea programada 1
 *
 *      En este archivo se especifica la solución a la tarea programada 1.
 * Continuar con la descripción del archivo aquí.
 *
 * @author  Geovanny Cordero Valverde   (correo@ucr.ac.cr)
 * @author  Carlos Delgado Rojas        (carlos.delgadorojas@ucr.ac.cr)
 */


//  *** Firma de los métodos usados en el programa ***

void tp(int,char**);




//  *** Método de ejecución principal ***

/**
 * Método de ejecución principal del programa
 * @param argc
 * @param argv
 * @return 0 si el proceso es exitoso, -1 en caso contrario
 */
int main(int argc, char** argv){

    tp(argc, argv);

    return 0;
}



//  *** Implementación de los métodos usados ***

/**
 * @brief Método que llena e imprime las matrices A y B
 *
 * 		Este método llena las matrices A y B, con los números entre 0 y 5 para A
 * y los números entre 0 y 2 para B. Luego, imprime A y B para verificar el resultado.
 *
 * @param n dimensión de las matrices (siempre son nxn)
 * @param A direción en memoria de la matriz A
 * @param B direción en memoria de la matriz B
 */
void llenarMatrices(int n, int* A, int* B){
	int i=0;
	for(; i < (n*n) ; i++){
        A[i] = rand() % 6; // llena la matriz A con números entre 0 y 5
    }

	i = 0;
    for(; i< (n*n) ; i++){
        B[i] = rand() % 3; // llena la matriz B con números entre 0 y 2
    }

	// imprime A
	printf("A =");
	i=0;
	for (; i < n*n; i++) {
		if(!(i%n)){
			printf("\n\t");
		}
		printf("%i ", A[i]);
	}

    // imprime B
	printf("\nB =");
	i=0;
	for (; i < n*n; i++) {
		if(!(i%n)){
			printf("\n\t");
		}
		printf("%i ", B[i]);
	}
}

/**
 *
 * @param ac
 * @param av
 */
void tp(int ac, char** av){

    int myId; // identificador del proceso actual (#)
    int numProcs; // cantidad de procesos creados
    int nameLen; // nombre del nodo del clúster en el corre el proceso actual
    int n; // dimensión de la matriz (siempre es de tamaño n x n)
    int nPrimos; // cantidad de números primos en M

    double i_ttime;
    double f_ttime;
    double i_time;
    double f_time;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

	int* A; // Matriz con números aleatorios entre 0 y 5
	int* B; // Matriz con números aleatorios entre 0 y 2
	int* M; // Matriz M = A * B
	int* P; // Vector donde P[i] = cantidad de # primos en la columna i de M
	int* C; // Matriz donde C[i,j] = M[i,j] + M[i,j-1] + M[i-1,j] + M[i,j+1] + M[i+1,j]
	

	//se inicia el trabajo con MPI
    MPI_Init(&ac, &av); 

	//guarda el tiempo actual, en segundos
    i_ttime = MPI_Wtime(); 

	//obtiene el id del proceso actual
    MPI_Comm_rank(MPI_COMM_WORLD, &myId);

	//guarda la cantidad de procesos 
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

	//guarda el nombre de la computadora en la que corre el proceso actual y el tamaño
    MPI_Get_processor_name(processor_name, &nameLen);

	printf("Proceso %i de %i corriendo en %s \n", myId, numProcs, processor_name);
	
	/* Barrera de sincronizacion.
	   Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(myId == 0) {
		printf("Ingrese la dimensión de su matriz: ");
		scanf("%i", &n);

		//guarda el tiempo actual, en segundos
		i_time = MPI_Wtime();

		srand(time(NULL));

		//llenado de matrices		
		A = (int *) malloc(n * n * sizeof(int)); // reserva de memoria para una matriz de enteros n x n
		B = (int *) malloc(n * n * sizeof(int));

		llenarMatrices(n,A,B);

	}
	/* BRETE */
	
	/*
	para calcular M se puede hacer de esta manera
	for(int i=0; i<n; i++){
    	    for(int k=0; k<n; k++){
	    	for(int j=0;j<n;j++){
		    m[i][k] += a[i][j] * b[j][k];
    		}
    	    }
    	}
    	*/
	
	/* BRETE */
		
	/* Barrera de sincronizacion.
	   Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	//guarda el tiempo actual, en segundos
    i_time = MPI_Wtime(); 
	
	if(myId == 0){
		printf("\nResultados finales:\n");
		printf("	Valor de n = %i \n", n);
		printf("	Número total de procesos que corrieron: %i \n", numProcs);
		printf("	Total de valores primos en M: %i \n", nPrimos);
			    
	    /*El tiempo que tardó desde que ya el usuario comunicó sus valores hasta 
		  antes de que se desplieguen resultados en pantalla y se escriban los 
		  archivos de texto. */
		printf("	Tiempo total de \"procesamiento\": %f segundos. \n", f_time - i_time);
		
		/* escribir en los archivos */
		
		/*Si n ≤ 100, se despliegan en pantalla también: A, B, M, P y C, de manera
		  que se puedan distinguir fácilmente las filas y las columnas en el caso
	      de las matrices. Si n > 100 se copia en un archivo de texto cada uno de 
	      estos arreglos. */
		//if(dim <100){
		//if(0){
//			printf("\n***	A	****\n");
//			for(i=0;i<dim;i++){
//				for(j=0;j<dim;j++){
//					printf(" %i ", A[i][j]);
//				}
//				printf("\n");
//			}
//
//			printf("\n\n***	B	****\n");
//			for(i=0;i<dim;i++){
//				for(j=0;j<dim;j++){
//					printf(" %i ", B[i][j]);
//				}
//				printf("\n");
//			}
//
//			//imprimir M
//			printf("\n\n***	M	****\n");
//			for(i=0;i<dim;i++){
//				for(j=0;j<dim;j++){
//					printf(" %i ", M[i][j]);
//				}
//				printf("\n");
//			}
//			//imprimir P
//			printf("\n\n***	P	****\n");
//			for(i=0;i<dim-i;i++){
//				printf("%i, ", P[i]);
//			}
//			printf("%i. \n", P[i]);
//
//
//			//imprimir C
//			printf("\n\n***	C	****\n");
//			for(i=0;i<dim;i++){
//				for(j=0;j<dim;j++){
//					printf(" %i ", C[i][j]);
//				}
//				printf("\n");
//			}
//		}
//	}
	
	//toma el tiempo al momento del final de ejecucion
    f_time = MPI_Wtime();
    printf("	Tiempo total: %f segundos. \n", f_ttime - i_ttime);
    
    MPI_Finalize();
}
}

