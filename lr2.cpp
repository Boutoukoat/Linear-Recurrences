
// ---------------------------------------------------------------------------------------
//
// Liear recurrences and primality tests 
// 
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "gmp.h"

static uint64_t exit_case[6] = { 0, 0, 0, 0, 0, 0 };

// ------------------------------------------------------------------------
// Linear recurrences second order
//
// S[n] = p*S[n-1]+q*S[n-2]
//
// return false : composite for sure
// return true : might be prime
// ------------------------------------------------------------------------

bool lnrcr2(int64_t p, int64_t q, mpz_t n)
{
	bool result = true;
	if (mpz_cmp_si(n, 11) < 0) {
		int64_t nn = mpz_get_si(n);
		if (nn == 1)
			return false;	// match the results of primecount
		if (nn == 2)
			return true;
		if (nn == 3)
			return true;
		if (nn == 5)
			return true;
		if (nn == 7)
			return true;
		return false;
	}
	if (mpz_tstbit(n, 0) == 0) {
		return false;
	}

	mpz_t ym, zm, y2, z2, yz, ty, tz, t;
	mpz_inits(ym, zm, y2, z2, yz, ty, tz, t, 0);
	mpz_t z1, zz1, zz2, yq, zq, yq1, zq1;
	mpz_inits(z1, zz1, zz2, yq, zq, yq1, zq1, 0);
	mpz_t ya1, ya2, za1, za2, ymyn;
	mpz_inits(ya1, ya2, za1, za2, ymyn, 0);

	do {
		uint64_t bit = mpz_sizeinbase(n, 2) - 1;
		mpz_set_ui(ym, 1);
		mpz_set_ui(zm, 0);
		while (bit) {
			bit -= 1;
			mpz_mul(y2, ym, ym);
			mpz_mul(z2, zm, zm);
			mpz_mul(yz, ym, zm);
			mpz_add(yz, yz, yz);

			mpz_mul_si(ty, y2, p);
			mpz_add(ty, ty, yz);
			mpz_mul_si(tz, y2, q);
			mpz_add(tz, tz, z2);

			if (mpz_tstbit(n, bit)) {
				mpz_set(t, ty);
				mpz_mul_si(ty, t, p);
				mpz_add(ty, ty, tz);
				mpz_mul_si(tz, t, q);
			}
			mpz_mod(ym, ty, n);
			mpz_mod(zm, tz, n);
		}

		// compute S[n]
		mpz_mul_si(z1, ym, p);
		mpz_add(z1, z1, zm);
		mpz_add(z1, z1, zm);
		mpz_mod(zz1, z1, n);

		// compute S[1]
		mpz_set_si(z2, p);
		mpz_mod(zz2, z2, n);

		// gmp_printf("S %Zd : %Zd %Zd\n", n, zz1, zz2);
		// check S[n] == S[1]
		if (mpz_cmp(zz1, zz2) != 0) {
			result = false;	// composite for sure
			exit_case[0] += 1;
			break;
		}
		if (q != 1) {

			// compute S[n+1]
			mpz_mul_si(yq1, ym, p);
			mpz_add(yq1, yq1, zm);
			mpz_mod(yq1, yq1, n);
			mpz_mul_si(zq1, ym, q);
			mpz_mod(zq1, zq1, n);

			// compute q^n
			mpz_mul(z1, ym, zq1);
			mpz_mul(z2, zm, yq1);
			mpz_sub(z1, z1, z2);
			mpz_mod(zz1, z1, n);

			// compute q^1
			mpz_set_si(z2, q);
			mpz_mod(zz2, z2, n);

			// gmp_printf("q^n %Zd : %Zd %Zd\n", n, zz1, zz2);
			// check q^n == q^1    (Fermat test)
			if (mpz_cmp(zz1, zz2) != 0) {
				result = false;	// composite for sure
				exit_case[1] += 1;
				break;
			}
		}
		// compute S[-n]
		mpz_set(yq, ym);
		mpz_mul_si(zq, ym, p);
		mpz_add(zq, zq, zm);
		mpz_mod(zq, zq, n);
		mpz_sub(zq, n, zq);
		mpz_mul_si(z1, yq, p);
		mpz_add(z1, z1, zq);
		mpz_add(z1, z1, zq);
		mpz_mod(zz1, z1, n);

		// compute S[-1] = (1, -p)
		mpz_set_si(yq, 1);
		mpz_set_si(zq, p);
		mpz_mod(zq, zq, n);
		mpz_sub(zq, n, zq);
		mpz_mul_si(z2, yq, p);
		mpz_add(z2, z2, zq);
		mpz_add(z2, z2, zq);
		mpz_mod(zz2, z2, n);

		// gmp_printf("-n %Zd : %Zd %Zd\n", n, zz1, zz2);
		// check S[-n] == S[-1]
		if (mpz_cmp(zz1, zz2) != 0) {
			result = false;	// composite for sure
			exit_case[2] += 1;
			break;
		}

		mpz_set(ya1, ym);
		mpz_set(za1, zm);
		mpz_set_ui(ya2, 1);
		mpz_set_ui(za2, 0);
		for (uint64_t i = 2; i < 3; i++) {
			if (mpz_cmp(ya1, ya2) == 0 && mpz_cmp(za1, za2) == 0) {
				exit_case[3] += 1;
				break;	// might be prime
			}
			// compute S[i*n] 
			mpz_mul(ymyn, ya1, ym);
			mpz_mul_si(ty, ymyn, p);
			mpz_mul(t, za1, ym);
			mpz_add(ty, t, ty);
			mpz_mul(t, ya1, zm);
			mpz_add(ty, t, ty);
			mpz_mod(ya1, ty, n);
			mpz_mul_si(tz, ymyn, q);
			mpz_mul(t, za1, zm);
			mpz_add(tz, t, tz);
			mpz_mod(za1, tz, n);
			mpz_mul_si(z1, ya1, p);
			mpz_add(z1, z1, za1);
			mpz_add(z1, z1, za1);
			mpz_mod(zz1, z1, n);

			// compute S[i]
			mpz_mul_si(ty, ya2, p);
			mpz_add(ty, ty, za2);
			mpz_mul_si(tz, ya2, q);
			mpz_mod(ya2, ty, n);
			mpz_mod(za2, tz, n);
			mpz_mul_si(z2, ya2, p);
			mpz_add(z2, z2, za2);
			mpz_add(z2, z2, za2);
			mpz_mod(zz2, z2, n);

			// gmp_printf("-i*%lu %Zd : %Zd %Zd\n", i, n, zz1, zz2);
			// check S[i*n] == S[i]
			if (mpz_cmp(zz1, zz2) != 0) {
				result = false;	// composite for sure
				exit_case[4] += 1;
				break;
			}
		}
		if (result == true)
			exit_case[5] += 1;
	} while (0);
	mpz_clears(zm, y2, z2, yz, ty, tz, t, 0);
	mpz_clears(ya1, ya2, za1, za2, ymyn, 0);
	mpz_clears(z1, zz1, zz2, yq, zq, yq1, zq1, 0);
	return result;
}

// test one second-order linear recurrences, display the pseudoprimes
int tmain(int argc, char **argv)
{
	uint64_t prime_count = 0;
	uint64_t composite_count = 0;
	uint64_t pseudoprime_count = 0;
	uint64_t error_count = 0;

	if (argc != 4) {
		printf("S[n] = p * S[n-1] + q * S[n-2], S[1] = 1, S[0] = 0;\n");
		printf("parameters p, q, testing limit\n");
		printf("e.g. p = 1 q = 1 limit = 170000 displays https://oeis.org/A005845\n");
	}

	mpz_t n, m;
	mpz_inits(m, n, 0);
	const int64_t p = atoll(argv[1]);
	const int64_t q = atoll(argv[2]);
	mpz_set_str(m, argv[3], 10);
	for (mpz_set_si(n, 1); mpz_cmp(n, m) <= 0; mpz_add_ui(n, n, 1)) {
		bool r = lnrcr2(p, q, n);
		bool s = mpz_probab_prime_p(n, 15) > 0;	// assume it is deterministic
		if (r && s)
			prime_count++;
		if (!r && !s)
			composite_count++;
		if (r && !s) {
			pseudoprime_count++;
			gmp_printf("pseudoprime %Zu\n", n);
		}
		if (s && !r) {
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
	// Branch coverage
	printf("Exit cases : ");
	for (unsigned i = 0; i < 6; i++)
		printf("%lu, ", exit_case[i]);
	printf("\n");
	return 0;
}

// test many second-order linear recurrences, display the number of pseudoprimes
int main(int argc, char **argv)
{
	uint64_t pseudoprime_count = 0;
	uint64_t error_count = 0;
	int64_t p, q;

	mpz_t n, m;
	mpz_inits(m, n, 0);
	mpz_set_ui(m, 1000000);
	gmp_printf("    p    q     #pseudoprimes <= %Zd\n", m);
	for (p = -10; p <= 10; p += 1) {
		for (q = -10; q <= 10; q += 1) {
			pseudoprime_count = 0;
			for (mpz_set_si(n, 1); mpz_cmp(n, m) <= 0; mpz_add_ui(n, n, 2)) {
				bool r = lnrcr2(p, q, n);
				bool s = mpz_probab_prime_p(n, 15) > 0;	// assume it is deterministic
				if (r && !s) {
					pseudoprime_count++;
				}
				if (s && !r) {
					error_count++;
					gmp_printf("false_negative %Zu !!!!! Wrong !\n", n);
				}
			}
			printf("%5ld %5ld %8lu\n", p, q, pseudoprime_count);
			fflush(stdout);
		}
	}
	mpz_clears(m, n, 0);

	printf("\n");
	printf("Error count .......... : %lu (MUST be zero)\n", error_count);
	// Branch coverage
	printf("Exit cases : ");
	for (unsigned i = 0; i < 6; i++)
		printf("%lu, ", exit_case[i]);
	printf("\n");
	return 0;
}




