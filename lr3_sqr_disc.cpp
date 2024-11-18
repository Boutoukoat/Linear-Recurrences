
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
//
// Just a look at what happens when the discriminant of the cubic polynomial is a square
//
int main(int argc, char **argv)
{
	mpz_t m[3], n, psp;
	mpz_inits(n, m[0], m[1], m[2], psp, 0);
	mpz_set_ui(m[0], 100000);
	mpz_set_ui(m[1], 1000000);
	mpz_set_ui(m[2], 10000000);
	mpz_set_ui(psp, 0);
	int64_t k = atoll(argv[1]);
	int64_t k_end = atoll(argv[2]);
	gmp_printf("k          : a %16s : %8Zd %8Zd %8Zd\n", "", m[0], m[1], m[2]);

	while (k <= k_end) {
		uint64_t pseudoprime_count = 0;
		int64_t a = 7 + k * (k - 1);

		printf("k %8ld : %18ld : ", k, a );
		mpz_set_ui(n, 3);
		mpz_set_ui(psp, 0);
		for (int j = 0; j < 3; j++) {
			for (; mpz_cmp(n, m[j]) <= 0; mpz_add_ui(n, n, 2)) {
				bool b = lnrcr3(0, a, a, n);
				if (b) {
					bool s = mpz_probab_prime_p(n, 25) > 0;	// assume it is deterministic
					if (b && !s) {
						pseudoprime_count++;
						if (mpz_cmp_ui(psp, 0) == 0) {
							mpz_set(psp, n);	// remember the first psp
						}
					}
				}
			}
			printf("%8ld ", pseudoprime_count);
			if (pseudoprime_count != 0)
				break;
		}
		gmp_printf(" (%Zd)\n", psp);
		k++;
	}
	mpz_clears(n, m, 0);

	// Branch coverage
	lnrcr3_debug();
	return 0;
}
