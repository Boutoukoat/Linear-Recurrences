
#pragma once

#include <stdint.h>
#include <gmp.h>

bool lnrcr2(int64_t p, int64_t q, mpz_t n);
bool lnrcr3(int64_t p, int64_t q, int64_t r, mpz_t n);
bool lnrcr4(int64_t p, int64_t q, int64_t r, int64_t s, mpz_t n);

void lnrcr2_debug(void);
void lnrcr3_debug(void);
void lnrcr4_debug(void);
