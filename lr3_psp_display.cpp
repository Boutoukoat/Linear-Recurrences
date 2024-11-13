
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

// test one third-order linear recurrences, display the pseudoprimes
int main(int argc, char **argv)
{
	uint64_t prime_count = 0;
	uint64_t composite_count = 0;
	uint64_t pseudoprime_count = 0;
	uint64_t error_count = 0;

	if (argc != 5) {
		printf("S[n] = p * S[n-1] + q * S[n-2] + r * S[n-3], S[2] = 0, S[1] = 1, S[0] = 0;\n");
		printf("input parameters : p, q, r, testing limit\n");
		printf("displays the pseudoprimes for this linear recurrence\n");
		printf("e.g. p = 1 q = 1 r = 0 limit = 170000 displays https://oeis.org/A005845\n");
		exit(1);
	}

	mpz_t n, m;
	mpz_inits(m, n, 0);
	const int64_t p = atoll(argv[1]);
	const int64_t q = atoll(argv[2]);
	const int64_t r = atoll(argv[3]);
	mpz_set_str(m, argv[4], 10);
	for (mpz_set_si(n, 1); mpz_cmp(n, m) <= 0; mpz_add_ui(n, n, 1)) {
		bool b = lnrcr3(p, q, r, n);
		bool s = mpz_probab_prime_p(n, 15) > 0;	// assume it is deterministic
		if (b && s) {
			prime_count++;
		}
		if (!b && !s) {
			composite_count++;
		}
		if (b && !s) {
			pseudoprime_count++;
			gmp_printf("pseudoprime %Zu\n", n);
		}
		if (s && !b) {
			error_count++;
			gmp_printf("false_negative %Zu !!!!! Wrong !\n", n);
		}
	}
	mpz_clears(m, n, 0);

	printf("\n");
	printf("Prime count .......... : %lu\n", prime_count);
	printf("Pseudoprime count .... : %lu\n", pseudoprime_count);
	printf("Composite count ...... : %lu\n", composite_count);
	printf("Error count .......... : %lu (MUST be zero)\n", error_count);
	printf("\n");

	// Branch coverage
	lnrcr3_debug();
	return 0;
}
