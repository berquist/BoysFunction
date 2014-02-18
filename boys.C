#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

template <typename T>
inline T Factorial(T n)
{
    assert(n<21); /* overflows 64-bit integer */

    T r=1;
    for (T i=1; i<=n; ++i) r *= i;
    return r;
}

template <typename T>
inline double InverseFactorial(T n)
{
    double r=1.0;
    for (T i=1; i<=n; ++i) r /= i;
    return r;
}

template <typename T>
inline T DoubleFactorial(T n)
{
    //assert(n<21); /* overflows 64-bit integer */

    T r=1;
    for (T i=1; i<=n; ++i) r *= (2*i - 1);
    return r;
}

template <typename T>
inline double InverseDoubleFactorial(T n)
{
    //assert(n<21); /* overflows 64-bit integer */

    double r=1.0;
    for (T i=1; i<=n; ++i) r /= (2*i - 1);
    return r;
}

template <typename T>
inline T Sign(T k)
{
    return ( k>0 ? T(1) : T(-1) );
}

template <typename T>
inline T SignPow(T k)
{
    // return ( k%2==0 ? T(1) : T(-1) );
    // return ( k%2>0 ? T(-1) : T(1) );
    return ( T( 1-2*(k%2) ) );
}

template <typename T, typename F>
inline F BoysAsymp(T n, F x)
{
    return DoubleFactorial(2*n-1)/pow(2,n+1) * sqrt( M_PI / pow(x,2*n+1) );
}

template <typename T, typename F>
inline F BoysTerm(T n, T k, F x)
{
    //return pow(-x,k) / ( Factorial(k) * (2*k + 2*n + 1) );
    return pow(-x,k) * InverseFactorial(k) / (2*k + 2*n + 1);
}

int main(int argc, char* argv[])
{
    int n = ( argc>1 ? atoi(argv[1]) : 0 );
    double x  = ( argc>2 ? atof(argv[2]) : 0.0 );

    printf("n = %d x = %lf \n", n, x );

    int k;
    double xpower;
    double invkfact;
    double knsum;
    double bterm_even;
    double bterm_odd;
    double fast;      

#if DEBUG
    fflush(stdout);
    printf("===================================================\n");

    //printf("%4s %2s %10s %14s %14s %14s %14s %2s \n", "n", "k", "x", "1/k!", "x^k", "k,BoysTerm(n,k,x)", "sum" );
    if (x<10.0) printf("%3s %3s %10s %14s %14s %14s %14s      \n", "n", "k", "x", "k!", "x^k", "BoysTerm(n,k,x)", "sum" );
    else        printf("%3s %3s %10s %14s %14s %14s %14s %14s \n", "n", "k", "x", "k!", "x^k", "BoysTerm(n,k,x)", "sum", "BoysAsymp(n,x)" );
    double b = 0.0, s = 0.0;
    k = 0;
    while( (k<200) && (fabs(b=BoysTerm(n,k,x))>1.0e-16) && (fabs(b)<1.0e10))
    {
        //b = BoysTerm(n,k,x);
        s += b;
        //printf("%4ld %2ld %10.5lf %18ld %14.7e %14.7e %14.7e %2ld %c \n", n, k, x, Factorial(k), pow(-x,k), b, s, k+1, fabs(b)<1.0e-14 ? '*' : ' ' );
        if (x<10.0) printf("%3d %3d %10.5lf %14.7e %14.7e %14.7e %14.7e %3d %c \n", 
                            n, k, x, InverseFactorial(k), pow(x,k), b, s, k+1, 
                            fabs(b)<1.0e-14 ? '*' : ' ' );
        else        printf("%3d %3d %10.5lf %14.7e %14.7e %14.7e %14.7e %14.7e %3d %c %s \n", 
                            n, k, x, InverseFactorial(k), pow(x,k), b, s, BoysAsymp(n,x), k+1, 
                            fabs(b)<1.0e-14 ? '*' : ' ' , fabs(BoysAsymp(n,x)-s)<1.0e-10 ? "**" : "  " );
        k++;
    }

    fflush(stdout);
    printf("===================================================\n");

    printf("%4s %2s %10s %12s %12s %14s %14s %14s %14s %14s %14s %14s %14s \n", "n", "k", "x", "2n+2k+1", "fast 2n+2k+1", "1/k!", "fast 1/k!","x^k", "fast x^k", "BoysTerm", "BoysTermFast", "sum", "fast sum" );

    k = 0;

    double sum;

    xpower      = 1.0;
    invkfact    = 1.0;
    knsum       = 2*n+1;
    bterm_even  = 1/knsum;
    bterm_odd   = 0.0;
    sum         = BoysTerm(n,k,x);
    fast        = 1/knsum;

    printf("%4d %2d %10.5lf %12.6e %12.6e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e \n", 
            n, k, x, 2*n+2*k+1.0, knsum, InverseFactorial(k), invkfact, pow(x,k), xpower, BoysTerm(n,k,x), bterm_even, sum, fast);

    for (k=1; k<20; ++k)
    {
        sum      += BoysTerm(n,k,x);

        xpower    *= x;
        invkfact  /= k;                           /* use table lookup for division by integers on PPC */
        knsum     += 2;                           /* cast should be done by compiler */
        bterm_odd  = -xpower * invkfact / knsum ; /* use table lookup for division by integers on PPC ??? - matrix or compute for every n */
        fast      += bterm_odd;
        printf("%4d %2d %10.5lf %12.6e %12.6e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e \n", 
                n, k, x, 2*n+2*k+1.0, knsum, InverseFactorial(k), invkfact, pow(x,k), xpower, BoysTerm(n,k,x), bterm_odd, sum, fast);

        ++k;

        sum      += BoysTerm(n,k,x);

        xpower    *= x;
        invkfact  /= k;                           /* use table lookup for division by integers on PPC */
        knsum     += 2;                           /* cast should be done by compiler */
        bterm_even =  xpower * invkfact / knsum ; /* use table lookup for division by integers on PPC ??? - matrix or compute for every n */
        fast      += bterm_even;

        printf("%4d %2d %10.5lf %12.6e %12.6e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e \n", 
                n, k, x, 2*n+2*k+1.0, knsum, InverseFactorial(k), invkfact, pow(x,k), xpower, BoysTerm(n,k,x), bterm_even, sum, fast);
    }
    printf("===================================================\n");
    fflush(stdout);
#endif

    printf("%4s %20s %20s \n", "k", "term", "sum" );
    printf("===================================================\n");

    k           = 0;
    xpower      = 1.0;
    invkfact    = 1.0;
    knsum       = 2*n+1;
    bterm_even  = 1.0/knsum;
    bterm_odd   = 0.0;

    double * fast_vec = (double *) malloc( 1001*sizeof(double) );
    fast_vec[0] = bterm_even;

    for (k=1; k<1000; ++k )
    {
        xpower      *= x;
        invkfact    /= k;                           /* use table lookup for division by integers on PPC */
        knsum       += 2;                           /* cast should be done by compiler */
        fast_vec[k]  = -xpower * invkfact / knsum ; /* use table lookup for division by integers on PPC ??? - matrix or compute for every n */
        ++k;
        xpower      *= x;
        invkfact    /= k;                           /* use table lookup for division by integers on PPC */
        knsum       += 2;                           /* cast should be done by compiler */
        fast_vec[k]  =  xpower * invkfact / knsum ; /* use table lookup for division by integers on PPC ??? - matrix or compute for every n */
    }

    fast = 0.0;
    k = 0;
    while ( (k<1000) && (fabs(fast_vec[k])>1e-14) )
    {
        fast += fast_vec[k];
        printf("%4d %24.14e %24.14e \n",k, fast_vec[k], fast );
        ++k;
        if (fast > 1.0e12 ) break;
    }

    if (fast > 1.0/(2*n+1) )
        printf("%lf cannot be greater than %lf \n", fast, 1.0/(2*n+1) ); 

    printf("===================================================\n");
    fflush(stdout);

    return(0);
}
