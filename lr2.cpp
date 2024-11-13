
// ---------------------------------------------------------------------------------------
//
// Liear recurrences and primality tests 
// 
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "lnrcr.h"

static uint64_t exit_case[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

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
		mpz_set_si(ym, 1);
		mpz_set_si(zm, 0);
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
		mpz_set_si(ya2, 1);
		mpz_set_si(za2, 0);
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

void lnrcr2_debug(void)
{
	// Branch coverage and debug
	printf("Exit cases : ");
	for (unsigned i = 0; i < 9; i++)
		printf("%lu, ", exit_case[i]);
	for (unsigned i = 0; i < 9; i++)
		exit_case[i] = 0;
	printf("\n");
}
