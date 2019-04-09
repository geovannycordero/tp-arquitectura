#include <stdio.h>
#include <unistd.h> //para sleep()
#include <time.h>   //para el random
#include <stdlib.h>
#include <math.h>
#include "mpi.h"

/**
 * @brief Archivo en C con la solución a la tarea programada 1
 *
 * En este archivo se especifica la solución a la tarea programada 1.
 * Continuar con la descripción del archivo aquí.
 *
 * @author  Geovanny Cordero Valverde   (geovanny.corderovalverde@ucr.ac.cr)
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
            printf("\n");
        }
        printf("\t%i ", A[i]);
    }

    // imprime B
    printf("\nB =");
    i=0;
    for (; i < n*n; i++) {
        if(!(i%n)){
            printf("\n");
        }
        printf("\t%i ", B[i]);
    }
}

/**
 * @brief Método utilizado por cada proceso para calcular su parte de M
 *
 *      Agregar una descripción más amplia aquí
 * También aprovechar esta línea ...
 *
 * @param Mx Buffer para el cálculo parcial de M
 * @param Ax Buffer con el subconjunto de filas de la matriz A
 * @param B Matriz B
 * @param filas cantidad de filas correspondientes al proceso
 * @param columnas las n columnas
 */
void calcularM(int* Mx, int* Ax, int* B, int filas, int columnas){
    int fila=0;
    for (; fila < filas; ++fila) {
        int columna = 0;
        for(; columna < columnas; ++columna){
            int offset = fila * columnas + columna; // índice de acceso de la matriz Mx(fila,columna) = Mx[offset]
            Mx[offset] = 0; // limpia la basura en memoria
            int x = 0; //los n elementos de cada fila
            for(; x < columnas; ++x){
                Mx[offset] = Mx[offset] + Ax[fila * columnas + x] * B[x * columnas + columna];
            }
        }
    }
}

/**
 * @brief
 *
 * @param num
 * @return
 */
int esPrimo(int num){
    int i;
    if(num<=3 && num!=0){
        return 1; // num = 0 o 1 o 2 o 3
    } else{
        i = 2;
        for(; i <= (int)sqrt(num) ; i++){
            if(!(num%i)){ // si num%i = 0 entonces no es primo
                return 0;
            }
        }
        return 1; // es primo
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
    int tp; // cantidad de números primos en M (suma de todos los tpx)
    int tpx; // cantidad de números detectados por cada proceso
    int i, j; //contadores genéricos

	double i_time, f_time;
    double i_ttime, f_ttime;
    
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    int* A; // Matriz con números aleatorios entre 0 y 5
    int* Ax; // porción de la matriz A que le corresponde a cada proceso

    int* B; // Matriz con números aleatorios entre 0 y 2

    int* M; // Matriz M = A * B
    int* Mx; // porción de la matriz M de resultados
    int* My; // porción de la matriz M ampliada (para calcular C)

    int* P; // Vector donde P[i] = cantidad de # primos en la columna i de M
    int* Px; // porción del vector P para cada proceso

    int* C; // Matriz donde C[i,j] = M[i,j] + M[i,j-1] + M[i-1,j] + M[i,j+1] + M[i+1,j]
    int* Cx; // porción de la matriz C que le corresponde calcular a cada proceso


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

    // Proceso 0 se encarga de crear las matrices A y B
    if(myId == 0) {
        printf("Ingrese la dimensión de su matriz: ");
        scanf("%i", &n);

        //guarda el tiempo actual, en segundos
        i_time = MPI_Wtime();

        srand(time(NULL));

        //reserva de memoria
        A = (int *) malloc(n * n * sizeof(int)); // matriz de enteros n x n
        B = (int *) malloc(n * n * sizeof(int));
        M = (int *) malloc(n * n * sizeof(int));

        P = (int *) malloc(n * sizeof(int));     // vector de n entradas
        C = (int *) malloc(n * n * sizeof(int));

        // verifica que efectivamente se haya reservado la memoria necesaria
        if (A==NULL || B==NULL || M == NULL || P==NULL || C==NULL){
            fprintf(stderr, "ERROR: La aplicación no pudo reservar memoria para las matrices, por lo que se ha cerrado\n");
            exit(EXIT_FAILURE);
        } else{
            //llena las matrices con valores aleatorios
            llenarMatrices(n,A,B);
        }
    }

    // propagación del valor n definido por el usuario
    MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);

    //reserva de memoria
    Ax = (int*) malloc(n/numProcs * n * sizeof(int));
    Mx = (int*) malloc(n/numProcs * n * sizeof(int));
    Px = (int*) calloc(n , sizeof(int));
    My = (int*) calloc(((n/numProcs)+2) * n, sizeof(int)); // matriz My, con dos filas y columnas extra
    Cx = (int*) calloc(n/numProcs * n , sizeof(int));

    if(myId!=0){
        B = (int *) malloc(n * n * sizeof(int));
        C = (int *) malloc(n * n * sizeof(int));
    }

    // propagación de la totalidad de la matriz B
    MPI_Bcast(B,(n*n),MPI_INT,0,MPI_COMM_WORLD);

    int filas = n/numProcs;
    int columnas = n;


    // Divide A en Ax
    MPI_Scatter(A,filas*columnas,MPI_INT,
                Ax,filas*columnas,MPI_INT,
                0,MPI_COMM_WORLD);

    // cada proceso calcula su parte de M (Mx)
    calcularM(Mx,Ax,B,filas,columnas);

    // Envía Mx para que se transforme en M
    MPI_Gather(Mx,filas*columnas,MPI_INT,
               M,filas*columnas,MPI_INT,
               0,MPI_COMM_WORLD);

    /* Barrera de sincronización.
       Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
    MPI_Barrier(MPI_COMM_WORLD);

    // Inicia el recorrido único de la matriz Mx

    if(myId == 0){
        // imprime M
        printf("\nM=");
        int i=0;
        for (; i < n*n; i++) {
            if(!(i%n)){
                printf("\n");
            }
            printf("\t%i ", M[i]);
        }
        printf("\n");
    }

    int* cuantos; //filas que le tocan a cada proceso
    int* inicio; //fila de inicio de cada proceso

    cuantos = malloc(numProcs * sizeof(int));
    inicio = malloc(numProcs * sizeof(int));

    if(myId == 0){
        j = 0;
        for( ; j < numProcs; ++j){
            if(j == 0){
                cuantos[j] = (filas + 1) * n;
                inicio[j] = 0;
            }
            else if(j == numProcs-1){
                cuantos[j] = (filas + 1) * n;
                inicio[j] = (filas * j - 1) * n;
            }
            else{
                cuantos[j] = (filas + 2) * n;
                inicio[j] = (filas * j - 1) * n;
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Scatterv(M, cuantos, inicio, MPI_INT, My, (filas+2)*n, MPI_INT, 0, MPI_COMM_WORLD);

    //setea el inicio del recorrido de la matriz My
    if(myId==0){
        i = 0;
    } else {
        i = n;
    }
    int des = 0;
    if(myId != 0){
        des = (filas*n)+n; //cantidad de desplazamientos
    } else{
        des = filas*n; //cantidad de desplazamientos
    }

    j = 0; // índice de escritura sobre Cx

    for( ; i < des; i++){
        if(myId==0 & i<n){ // primera fila
            if((i%n)==0){ //está sobre la columna 0
                Cx[j] = My[i] + My[i+1] + My[i+n];
                j=j+1;
            } else if ( i%n == n-1 ) { //está sobre la columna n-1
                Cx[j] = My[i] + My[i-1] + My[i+n];
                j=j+1;
            } else{ //entrada normal
                Cx[j] = My[i] + My[i-1] + My[i+1] + My[i+n];
                j=j+1;
            }
        } else{//resto de los procesos
            if((i%n)==0){ //está sobre la columna 0
                Cx[j] = My[i] + My[i-n] + My[i+1] + My[i+n];
                j=j+1;
            } else if ( i%n == n-1 ) { //está sobre la columna n-1
                Cx[j] = My[i] + My[i-n] + My[i-1] + My[i+n];
                j=j+1;
            } else{ //entrada normal
                Cx[j] = My[i] + My[i-n] + My[i-1] + My[i+1] + My[i+n];
                j=j+1;
            }
        }
        if(esPrimo(My[i])){
            Px[i%n] = Px[i%n] + 1; //suma a la columna del vector correspondiente
            tpx = tpx + 1; // suma al contador de números primos
        }
    }

    /* Barrera de sincronización.
       Hasta que todos los procesos alcancen este llamado ninguno puede proseguir.*/
    MPI_Barrier(MPI_COMM_WORLD);

    // Envía Cx para que se transforme en C
    MPI_Gather(Cx,filas*columnas,MPI_INT,
               C,filas*columnas,MPI_INT,
               0,MPI_COMM_WORLD);

    if(myId == 0){
        printf("\nC=");
        i = 0;
        for (; i < n*n; i++) {
            if(!(i%n)){
                printf("\n");
            }
            printf("\t%i ", C[i]);
        }
        printf("\n");
    }

    MPI_Reduce(&tpx,&tp,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); // envía los datos al proceso ROOT para tp
    MPI_Reduce(Px,P,n,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD); // envía los datos al proceso ROOT para vector P

    //guarda el tiempo actual, en segundos
    f_time = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);
    
    //

    if(myId == 0){
        // imprime el vector P para checkear los resultados
        int a = 0;
        printf("\n");
        for(;a<n;++a){
            printf("P[%i]: %i, ", a, P[a]);
        }

        printf("\nResultados finales:\n");
        printf("\tValor de n = %i \n", n);
        printf("\tNúmero total de procesos que corrieron: %i \n", numProcs);
        printf("\tTotal de valores primos en M: %i \n", tp);

        /*El tiempo que tardó desde que ya el usuario comunicó sus valores hasta
          antes de que se desplieguen resultados en pantalla y se escriban los
          archivos de texto. */
          printf("	Tiempo parcial: %f segundos. \n", f_time - i_time);

        //toma el tiempo al momento del final de ejecucion
        f_ttime = MPI_Wtime();
        printf("	Tiempo total: %f segundos. \n", f_ttime - i_ttime);

    }

    MPI_Finalize();
}

