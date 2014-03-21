#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>

static inline double Factorial(int n)
{
    double r = 1;
    for (int i=1; i<=n; ++i) {
        r *= i;
    }
    return r;
}

static inline double InverseFactorial(int n)
{
    double r=1.0;
    for (int i=1; i<=n; ++i) {
        r /= i;
    }
    return r;
}

static inline double DoubleFactorial(int n)
{
    double r=1;
    for (int i=1; i<=n; ++i) {
        r *= (2*i - 1);
    }
    return r;
}

static inline double PowerOfTwo(int n)
{
    /* This is not efficient and should use bit-shift instead... */
    double r=1;
    for (int i=1; i<=(n+1); ++i) {
        r *= 2;
    }
    return r;
}

/*
 * Helgaker and Taylor, equation 151
 * F_n(x) ~= (2n-1)!! / 2^(n+1) sqrt( PI / x^(2n+1) )
 *           ...for x large
 */
double BoysAsymptotic(double x, int n)
{
#ifdef DEBUG
    printf("=====================================\n");
#endif

    double f = 0.0;
    if (x<1.0e-9) {
         f = 1.0/(2*n+1);
    } else {
        double a = DoubleFactorial(2*n-1); /* (2n-1)!!              */
        double b = PowerOfTwo(n+1);        /* 2^(n+1)               */
        double c = pow(x,(double)(2*n+1)); /* x^(2n+1)              */
        double d = sqrt(M_PI / c);         /* sqrt( PI / x^(2n+1) ) */
        double e = a/b;                    /* (2n-1)!! / 2^(n+1)    */
        f = d*e;                           /* final expression      */
#ifdef DEBUG
        printf("(2n-1)!! = %e 2^(n+1) = %e x^(2n+1) = %e sqrt( PI / x^(2n+1) ) = %e (2n-1)!! / 2^(n+1) = %e \n",
                a, b, c, d, e);
#endif
    }

#ifdef DEBUG
    printf("BoysAsymptotic(%8.5lf,%d) = %20.10e \n", x, n, f);
#endif
    return f;

}

/*
 * Helgaker and Taylor, equation 168
 * F_n(x) = \Sum_{k=0}^{\infty} \frac{(-x)^k}{k!(2n+2k+1)}
 */
double BoysSum(double x, int n, double tol)
{
    /* The sum will diverge if x>4.0 */
    if (x>4.0) return sqrt(-1);

#ifdef DEBUG
        printf("=====================================\n");
#endif

    int    k = 0;
    double f = 0.0;
    double d = 0.0;
    do {
        double a = pow(x,(double)k);     /* x^k                        */
        if (k%2) a = -a;                 /* (-x)^k                     */
        double b = InverseFactorial(k);  /* 1/k!                       */
        double c = (2*n+2*n+1);          /* 2n+2k+1                    */
        d = a*b/c;                       /* (-x)^k * 1/k! / (2n+2k+1)) */
        f += d;                          /* accumulate term to sum     */
#ifdef DEBUG
        printf("k = %d (-x)^k = %e k! = %e 2n+2k+1 = %e term = %e\n",
                k, a, b, c, d);
#endif
        ++k;
    } while (fabs(d)>tol && fabs(d)<1.0e7);

#ifdef DEBUG
    printf("BoysSum(%8.5lf,%d,%8.5e) = %20.10e \n", x, n, f, tol);
#endif
    return f;

}

int main(int argc, char* argv[])
{
    int n = (argc>1) ? atoi(argv[1]) : 0;

    for (int i=0; i<30; ++i) {
        double x = (double)i;
        double a = BoysAsymptotic(x,n);
        double s = BoysSum(x,n, 1.0e-13);
        printf("Boys(%8.5lf,%d) = %20.10e (asymp) %20.10e (sum) \n", x, n, a, s);
    }

    return 0;
}
