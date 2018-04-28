/********************************************
 *  Program for calculating Jacobi symbol   *
 *                                          *
 *  CG - August 2008                        *
 ********************************************/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int jacobi(int, int);

int main()
{
    int a, n;

    cout << "a? ";
    cin >> a;
    cout << "n? ";
    cin >> n;

    cout << "\nThe Jacobi symbol is " << jacobi(a,n) << endl;

    return 0;
}

/* Precondition: a, n >= 0; n is odd */
int jacobi(int a, int n) {
    int ans;

    if (a == 0)
    {
        ans = (n == 1) ? 1 : 0;           // pgcd (2)
        printf("// (2)\n");
	}
    else if (a == 2) {                   // si a = 2  (5)
        switch ( n % 8 ) {
            case 1:
            case 7:
                    ans = 1;
                    printf("// (5)\n");
                    break;
            case 3:
            case 5:
                    ans = -1;
                    printf("// (5)\n");
                    break;
        }
    }
    else if ( a >= n )
    {
        ans = jacobi(a%n, n);				//       (1)
        printf("// (1)\n");
	}
    else if ( a % 2 == 0 )					// Si a est pair  (3)
    {
        ans = jacobi(2,n)*jacobi(a/2, n);
        printf("// (3)\n");
	}
    else                                   //         (6)
    {
        ans = ( a % 4 == 3 && n % 4 == 3 ) ? -jacobi(n,a) : jacobi(n,a);
        printf("// (6)\n");
	}
    return ans;
}
