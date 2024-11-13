
// ---------------------------------------------------------------------------------------
//
// Liear recurrences and primality tests 
// 
// ---------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdint.h>

#include "lnrcr.h"

static uint64_t exit_case[6] = { 0, 0, 0, 0, 0, 0 };

// ------------------------------------------------------------------------
// Linear recurrences third order
//
// S[n] = p*S[n-1]+q*S[n-2]+r*S[n-3] 
//
// return false : composite for sure
// return true : might be prime
// ------------------------------------------------------------------------

bool lnrcr3(int64_t p, int64_t q, int64_t r, mpz_t n)
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
		return false;	// even number
	}
	// precomputed constants and initialization sequence
	int64_t bp = p * p + q;
	int64_t bq = p * q + r;
	int64_t br = r * p;
	int64_t S0 = 3;
	int64_t S1 = p;
	int64_t S2 = q + bp;

	mpz_t t, xm, ym, zm, rt0, rt1, rt2, qx, qy, qz;
	mpz_inits(t, xm, ym, zm, rt0, rt1, rt2, qx, qy, qz, 0);
	mpz_t z1, zz1, z2, zz2;
	mpz_inits(z1, zz1, z2, zz2, 0);
	mpz_t xq0, yq0, zq0, xq1, yq1, zq1, xq2, yq2, zq2;
	mpz_inits(xq0, yq0, zq0, xq1, yq1, zq1, xq2, yq2, zq2, 0);
	mpz_t xm1, ym1, zm1, xa1, ya1, za1;
	mpz_inits(xm1, ym1, zm1, xa1, ya1, za1, 0);

	do {

		uint64_t bit = mpz_sizeinbase(n, 2) - 1;

		mpz_set_si(xm, 0);
		mpz_set_si(ym, 1);
		mpz_set_si(zm, 0);

		while (bit) {
			bit -= 1;

			// compute S[2n] from S[n]
			mpz_mul(rt0, xm, xm);
			mpz_mul(rt1, xm, ym);
			mpz_add(rt1, rt1, rt1);
			mpz_mul_si(qx, rt0, bp);
			mpz_mul_si(qy, rt0, bq);
			mpz_mul_si(qz, rt0, br);
			mpz_mul_si(t, rt1, p);
			mpz_add(qx, qx, t);
			mpz_mul_si(t, rt1, q);
			mpz_add(qy, qy, t);
			mpz_mul_si(t, rt1, r);
			mpz_add(qz, qz, t);
			mpz_mul(t, ym, ym);
			mpz_add(qx, qx, t);
			mpz_mul(t, xm, zm);
			mpz_add(qx, qx, t);
			mpz_add(qx, qx, t);
			mpz_mul(t, ym, zm);
			mpz_add(qy, qy, t);
			mpz_add(qy, qy, t);
			mpz_mul(t, zm, zm);
			mpz_add(qz, qz, t);

			if (mpz_tstbit(n, bit)) {
				// compute S[n+1] from S[n]
				mpz_set(rt2, qx);
				mpz_mul_si(qx, rt2, p);
				mpz_add(qx, qx, qy);
				mpz_mul_si(qy, rt2, q);
				mpz_add(qy, qy, qz);
				mpz_mul_si(qz, rt2, r);
			}
			// modular reduction once since intermediate numbers are O(n^2)
			mpz_mod(xm, qx, n);
			mpz_mod(ym, qy, n);
			mpz_mod(zm, qz, n);
		}

		// compute S[n]
		mpz_mul_si(z1, xm, S2);
		mpz_mul_si(t, ym, S1);
		mpz_add(z1, z1, t);
		mpz_mul_si(t, zm, S0);
		mpz_add(z1, z1, t);
		mpz_mod(zz1, z1, n);

		// compute S[1]
		mpz_set_si(z2, S1);
		mpz_mod(zz2, z2, n);

		// check S[n] == S[1]
		if (mpz_cmp(zz1, zz2) != 0) {
			exit_case[0] += 1;
			result = false;	// composite for sure
			break;
		}
		// compute S[n+0]
		mpz_set(xq0, xm);
		mpz_set(yq0, ym);
		mpz_set(zq0, zm);

		// compute S[n+1]
		mpz_mul_si(xq1, xq0, p);
		mpz_add(xq1, xq1, yq0);
		mpz_mul_si(yq1, xq0, q);
		mpz_add(yq1, yq1, zq0);
		mpz_mul_si(zq1, xq0, r);

		// compute S[n+2]
		mpz_mul_si(xq2, xq1, p);
		mpz_add(xq2, xq2, yq1);
		mpz_mul_si(yq2, xq1, q);
		mpz_add(yq2, yq2, zq1);
		mpz_mul_si(zq2, xq1, r);

		if (r != 1) {
			// compute r^n = matrix determinant
			mpz_mul(qx, yq1, zq0);
			mpz_mul(t, yq0, zq1);
			mpz_sub(qx, qx, t);
			mpz_mod(qx, qx, n);

			mpz_mul(qy, xq1, zq0);
			mpz_mul(t, xq0, zq1);
			mpz_sub(qy, qy, t);
			mpz_mod(qy, qy, n);
			mpz_sub(qy, n, qy);

			mpz_mul(qz, xq1, yq0);
			mpz_mul(t, xq0, yq1);
			mpz_sub(qz, qz, t);
			mpz_mod(qz, qz, n);

			mpz_mul(z1, xq2, qx);
			mpz_mul(t, yq2, qy);
			mpz_add(z1, z1, t);
			mpz_mul(t, zq2, qz);
			mpz_add(z1, z1, t);
			mpz_mod(zz1, z1, n);

			// compute r^1
			mpz_set_si(z2, r);
			mpz_mod(zz2, z2, n);

			// check r^n == r^1     (fermat-test)
			if (mpz_cmp(zz1, zz2) != 0) {
				exit_case[1] += 1;
				result = false;	// composite for sure
				break;
			}
		}
		// compute S[-n] from S[n], S[n+1], S[n+2]  (divisor r^n omitted)
		mpz_mul(qx, xq1, yq0);
		mpz_mul(t, xq0, yq1);
		mpz_sub(qx, qx, t);
		mpz_mod(qx, qx, n);

		mpz_mul(qy, xq0, yq2);
		mpz_mul(t, xq2, yq0);
		mpz_sub(qy, qy, t);
		mpz_mod(qy, qy, n);

		mpz_mul(qz, xq2, yq1);
		mpz_mul(t, xq1, yq2);
		mpz_sub(qz, qz, t);
		mpz_mod(qz, qz, n);

		mpz_mul_si(z1, qx, S2);
		mpz_mul_si(t, qy, S1);
		mpz_add(z1, z1, t);
		mpz_mul_si(t, qz, S0);
		mpz_add(z1, z1, t);
		mpz_mod(zz1, z1, n);

		// compute S[-1] (divisor r omitted) as 1,-p,-q
		mpz_set_si(xq1, 1);
		mpz_set_si(yq1, p);
		mpz_mod(yq1, yq1, n);
		mpz_sub(yq1, n, yq1);
		mpz_set_si(zq1, q);
		mpz_mod(zq1, zq1, n);
		mpz_sub(zq1, n, zq1);

		mpz_mul_si(z2, xq1, S2);
		mpz_mul_si(t, yq1, S1);
		mpz_add(z2, z2, t);
		mpz_mul_si(t, zq1, S0);
		mpz_add(z2, z2, t);
		mpz_mod(zz2, z2, n);

		// check S[-n] * r == S[-1] * r^n  with r == r^n already verified
		if (mpz_cmp(zz1, zz2) != 0) {
			exit_case[2] += 1;
			result = false;	// composite for sure
			break;
		}
#define add(rx, ry, rz, ax, ay, az, bx, by, bz, tx, ty, tz) \
	do { \
			mpz_mul(rt0, ax, bx); \
			mpz_mul(rt1, bx, ay); \
			mpz_mul(t, by, ax);   \
			mpz_add(rt1, rt1, t); \
			mpz_mul_si(tx, rt0, bp);  \
			mpz_mul_si(ty, rt0, bq);  \
			mpz_mul_si(tz, rt0, br);  \
			mpz_mul_si(t, rt1, p);    \
			mpz_add(tx, tx, t);       \
			mpz_mul_si(t, rt1, q);    \
			mpz_add(ty, ty, t);       \
			mpz_mul_si(t, rt1, r);    \
			mpz_add(tz, tz, t);       \
			mpz_mul(t, bx, az);       \
			mpz_add(tx, tx, t);       \
			mpz_mul(t, ax, bz);       \
			mpz_add(tx, tx, t);       \
			mpz_mul(t, ay, by);       \
			mpz_add(tx, tx, t);       \
			mpz_mul(t, by, az);       \
			mpz_add(ty, ty, t);       \
			mpz_mul(t, ay, bz);       \
			mpz_add(ty, ty, t);       \
			mpz_mul(t, az, bz);       \
			mpz_add(tz, tz, t);       \
			mpz_mod(rx, tx, n); \
			mpz_mod(ry, ty, n); \
			mpz_mod(rz, tz, n); \
	} while (0)

		// S[-n]
		mpz_set(xm1, qx);
		mpz_set(ym1, qy);
		mpz_set(zm1, qz);

		// S[-1]
		mpz_set(xa1, xq1);
		mpz_set(ya1, yq1);
		mpz_set(za1, zq1);

		for (unsigned i = 2; i < 4; i++) {
			if (mpz_cmp(xa1, xm1) == 0 && mpz_cmp(ya1, ym1) == 0 && mpz_cmp(za1, zm1) == 0) {
				break;
			}
			// compute S[-i*n] = add(S[(i+1)*n], S[-n])
			add(xm1, ym1, zm1, xm1, ym1, zm1, qx, qy, qz, xq2, yq2, zq2);

			mpz_mul_si(z1, xm1, S2);
			mpz_mul_si(t, ym1, S1);
			mpz_add(z1, z1, t);
			mpz_mul_si(t, zm1, S0);
			mpz_add(z1, z1, t);
			mpz_mod(zz1, z1, n);

			// compute S[-i]  from S[-i + 1]
			mpz_set(xq1, za1);
			mpz_mul_si(yq1, xa1, r);
			mpz_mul_si(t, za1, p);
			mpz_mod(t, t, n);
			mpz_sub(t, n, t);
			mpz_add(yq1, yq1, t);
			mpz_mul_si(zq1, ya1, r);
			mpz_mul_si(t, za1, q);
			mpz_mod(t, t, n);
			mpz_sub(t, n, t);
			mpz_add(zq1, zq1, t);
			mpz_set(xa1, xq1);
			mpz_mod(ya1, yq1, n);
			mpz_mod(za1, zq1, n);

			mpz_mul_si(z2, xa1, S2);
			mpz_mul_si(t, ya1, S1);
			mpz_add(z2, z2, t);
			mpz_mul_si(t, za1, S0);
			mpz_add(z2, z2, t);
			mpz_mod(zz2, z2, n);

			// check S[-i*n] == S[-i]
			if (mpz_cmp(zz1, zz2) != 0) {
				exit_case[3] += 1;
				result = false;
				break;
			}
		}
		if (result == false)
			break;

		// S[n]
		mpz_set(xa1, xm);
		mpz_set(ya1, ym);
		mpz_set(za1, zm);

		// S[1]
		mpz_set_si(xq1, 0);
		mpz_set_si(yq1, 1);
		mpz_set_si(zq1, 0);

		for (unsigned i = 2; i < 4; i++) {
			if (mpz_cmp(xa1, xq1) == 0 && mpz_cmp(ya1, yq1) == 0 && mpz_cmp(za1, zq1) == 0) {
				break;
			}
			// compute S[i*n] = add(S[(i-1)*n], S[n])
			add(xa1, ya1, za1, xa1, ya1, za1, xm, ym, zm, xq2, yq2, zq2);

			mpz_mul_si(z1, xa1, S2);
			mpz_mul_si(t, ya1, S1);
			mpz_add(z1, z1, t);
			mpz_mul_si(t, za1, S0);
			mpz_add(z1, z1, t);
			mpz_mod(zz1, z1, n);

			// compute S[i] from S[i-1]
			mpz_mul_si(qx, xq1, p);
			mpz_add(qx, qx, yq1);
			mpz_mul_si(qy, xq1, q);
			mpz_add(qy, qy, zq1);
			mpz_mul_si(qz, xq1, r);
			mpz_mod(xq1, qx, n);
			mpz_mod(yq1, qy, n);
			mpz_mod(zq1, qz, n);

			mpz_mul_si(z2, xq1, S2);
			mpz_mul_si(t, yq1, S1);
			mpz_add(z2, z2, t);
			mpz_mul_si(t, zq1, S0);
			mpz_add(z2, z2, t);
			mpz_mod(zz2, z2, n);

			// check S[i*n] == S[i]
			if (mpz_cmp(zz1, zz2) != 0) {
				exit_case[4] += 1;
				result = false;
				break;
			}
		}
		if (result == false)
			break;
	} while (0);

	mpz_clears(t, xm, ym, zm, rt0, rt1, rt2, qx, qy, qz, 0);
	mpz_clears(z1, zz1, z2, zz2, 0);
	mpz_clears(xq0, yq0, zq0, xq1, yq1, zq1, xq2, yq2, zq2, 0);
	mpz_clears(xm1, ym1, zm1, xa1, ya1, za1, 0);

	return result;
}

void lnrcr3_debug(void)
{
	// Branch coverage and debug
	printf("Exit cases : ");
	for (unsigned i = 0; i < 6; i++)
		printf("%lu, ", exit_case[i]);
	for (unsigned i = 0; i < 6; i++)
		exit_case[i] = 0;
	printf("\n");
}
