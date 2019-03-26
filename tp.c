#include <stdio.h>
//para sleep()
#include <unistd.h> 
//para el random
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

void tp(int,char**);

int main(int argc, char** argv){

    tp(argc, argv);

    return 0;
}


void tp(int ac, char** av){

    int i, j, myid, numprocs, namelen, dim, nprimos;
    double i_ttime, f_ttime, i_time, f_time;;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
        
	//creación de matrices
	int a[100][100];
	int b[100][100];
	int m[100][100];
	int p[100];
	int c[100][100];
	

	//se inicia el trabajo con MPI
    MPI_Init(&ac, &av); 

	//guarda el tiempo actual, en segundos
    i_ttime = MPI_Wtime(); 

	//obtiene el id del proceso actual
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	//guarda la cantidad de procesos 
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	//guarda el nombre de la computadora en la que corre el proceso actual y el tamaño
    MPI_Get_processor_name(processor_name, &namelen); 

	printf("Proceso %i de %i corriendo en %s \n", myid, numprocs, processor_name);
	
	/* Barrera de sincronizacion.
	   Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	if(myid == 0){
		printf("Ingrese la dimensión de su matriz: ");
		scanf("%i", &dim);
		
		//guarda el tiempo actual, en segundos
	    i_time = MPI_Wtime(); 

		srand(time(NULL));
		
		//llenado de matrices		
		for(i=0;i<dim;i++){
			for(j=0;j<dim;j++){
				//números entre 0 y 2
				b[i][j] = rand() % 6;
			}
		}
				
		for(i=0;i<dim;i++){
			for(j=0;j<dim;j++){
				//números entre 0 y 2
				a[i][j] = rand() % 3;
			}
		}

	}
	
	/* BRETE HARDCORE */
	
	/* BRETE HARDCORE */
		
	/* Barrera de sincronizacion.
	   Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
	MPI_Barrier(MPI_COMM_WORLD);
	
	//guarda el tiempo actual, en segundos
    i_time = MPI_Wtime(); 
	
	if(myid == 0){
		printf("\nResultados finales:\n");
		printf("	Valor de n = %i \n", dim);
		printf("	Número total de procesos que corrieron: %i \n", numprocs);
		printf("	Total de valores primos en M: %i \n", nprimos);
			    
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
		if(0){	
			printf("\n***	A	****\n");
			for(i=0;i<dim;i++){
				for(j=0;j<dim;j++){
					printf(" %i ", a[i][j]);
				}
				printf("\n");
			}
				
			printf("\n\n***	B	****\n");
			for(i=0;i<dim;i++){
				for(j=0;j<dim;j++){
					printf(" %i ", b[i][j]);
				}
				printf("\n");
			}
			
			//imprimir M
			printf("\n\n***	M	****\n");
			for(i=0;i<dim;i++){
				for(j=0;j<dim;j++){
					printf(" %i ", m[i][j]);
				}
				printf("\n");
			}
			//imprimir P
			printf("\n\n***	P	****\n");
			for(i=0;i<dim-i;i++){
				printf("%i, ", p[i]);
			}
			printf("%i. \n", p[i]);
		
		
			//imprimir C
			printf("\n\n***	C	****\n");
			for(i=0;i<dim;i++){
				for(j=0;j<dim;j++){
					printf(" %i ", c[i][j]);
				}
				printf("\n");
			}
		}
	}	
	
	//toma el tiempo al momento del final de ejecucion
    f_time = MPI_Wtime();
    printf("	Tiempo total: %f segundos. \n", f_ttime - i_ttime);
    
    MPI_Finalize();
}

