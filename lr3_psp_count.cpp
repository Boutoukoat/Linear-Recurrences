
// ---------------------------------------------------------------------------------------
//
// Liear recurrences and primality tests 
// 
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "lnrcr.h"

static uint64_t exit_case[6] = { 0, 0, 0, 0, 0, 0 };

// ------------------------------------------------------------------------
// Linear recurrences third order
//
// S[n] = p*S[n-1]+q*S[n-2]+r*S[n-3] 
//
// ------------------------------------------------------------------------

// test many third-order linear recurrences, display the number of pseudoprimes
int main(int argc, char **argv)
{
	uint64_t pseudoprime_count = 0;
	uint64_t error_count = 0;
	int64_t p, q, r;

	mpz_t n, m;
	mpz_inits(m, n, 0);
	mpz_set_ui(m, 1000000);
	gmp_printf("    p    q    r     #pseudoprimes <= %Zd\n", m);
	fflush(stdout);
	for (p = 1; p <= 10; p += 1) {
		for (q = -10; q <= 10; q += 1) {
			for (r = -10; r <= 10; r += 1) {
				pseudoprime_count = 0;
				for (mpz_set_si(n, 1); mpz_cmp(n, m) <= 0; mpz_add_ui(n, n, 2)) {
					bool b = lnrcr3(p, q, r, n);
					bool s = mpz_probab_prime_p(n, 15) > 0;	// assume gmp primality test is deterministic
					if (b && !s) {
						pseudoprime_count++;
					}
					if (s && !b) {
						error_count++;
						gmp_printf("false_negative %Zu !!!!! Wrong !\n", n);
					}
				}
				printf("%5ld %5ld %5ld %8lu\n", p, q, r, pseudoprime_count);
				fflush(stdout);
			}
		}
	}
	mpz_clears(m, n, 0);

	printf("\n");
	printf("Error count .......... : %lu (MUST be zero)\n", error_count);
	printf("\n");
	// Branch coverage (internal debug)
	lnrcr3_debug();

	return 0;
}
