#include <stdio.h>
#include <stdlib.h>

int calculo(int k ){
    int f = k, c = k;

    int *m1, *m2;

    m1 = (int*) malloc(f * c * sizeof(int));
    m2 = (int*)	malloc(f * c * sizeof(int));
    
    int i = 0, j;
    for( ; i < f; ++i){
	j = 0;
	for( ; j < c; ++j){
	    m1[i * j] = 0;
	}
    }

    m1[1 * 1] = 66;

    i = 0;
    for( ; i < f; ++i){
        j = 0;
        for( ; j < c; ++j){
            printf(" %i ", m1[i * j]);
        }
	printf("\n");
    }

    free(m1);

    return 0;
}


int 
main(int argc, char** argv){
	calculo(atoi(argv[1]));
	printf("\nThe End!\n");	
	return 0;
}
