
// ---------------------------------------------------------------------------------------
//
// Liear recurrences and primality tests 
// 
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "lnrcr.h"

// ------------------------------------------------------------------------
// Linear recurrences third order
//
// S[n] = p*S[n-1]+q*S[n-2]+r*S[n-3] 
//
// ------------------------------------------------------------------------

// test primality with third-order linear recurrences
int main(int argc, char **argv)
{
	if (argc <= 5) {
		printf("S[n] = p * S[n-1] + q * S[n-2] + r * S[n-3];\n");
		printf("input parameters : p, q, r, n1, n2, n3, n4 ....\n");
		printf("e.g. p = 0 q = 1 r = 1 display the restricted pseudoprimes\n");
		exit(1);
	}

	mpz_t n;
	mpz_inits(n, 0);
	const int64_t p = atoll(argv[1]);
	const int64_t q = atoll(argv[2]);
	const int64_t r = atoll(argv[3]);
	for (int i = 4; i < argc; i++) {
		mpz_set_str(n, argv[i], 10);
		bool b = lnrcr3(p, q, r, n);
		bool s = mpz_probab_prime_p(n, 25) > 0;	// assume it is deterministic
		gmp_printf("%Zd : ", n);
		if (b && s) {
			printf("prime\n");
		}
		if (!b && !s) {
			printf("composite\n");
		}
		if (b && !s) {
			printf("pseudoprime\n");
		}
		if (s && !b) {
			printf("false_negative !!!!! Wrong !\n");
		}
	}
	mpz_clears(n, 0);

	// Branch coverage
	lnrcr3_debug();
	return 0;
}
