#include <stdio.h>


int prime()
{
    int n, i, flag = 0;

    printf("Enter a positive integer: ");
    scanf("%d", &n);

    for(i = 2; i <= n/2; ++i)
    {
        // condition for nonprime number
        if(n%i == 0)
        {
            flag = 1;
        }
    }
    
    return flag;
}


int main(){
	if(prime()){
		printf("Es primo");
	}
	else{
		printf("No es primo");
	}
	
	return 0;
}