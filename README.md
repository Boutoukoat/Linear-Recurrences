# Linear-Recurrences with Constant Coefficients

# First Order Linear recurrences

```
S[n] = p * S[n-1]
```

# Second Order Linear recurrences

```
S[n] = p * S[n-1] + q * S[n-1]
```

# Third Order Linear recurrences

```
S[n] = p * S[n-1] + q * S[n-1] + r *S[n-3]
```

# Fourth Order Linear recurrences

```
S[n] = p * S[n-1] + q * S[n-1] + r * S[n-3] + s * S[n-4]
```

# Fifth Order Linear recurrences

```
S[n] = p * S[n-1] + q * S[n-1] + r * S[n-3] + s * S[n-4] + t * S[n-5]
```

Linear recurrences are found in many domains : biology, physics, mathematics.

https://en.wikipedia.org/wiki/Linear_recurrence_with_constant_coefficients

https://en.wikipedia.org/wiki/Linear_recurrence_with_constant_coefficients#Order_1

https://en.wikipedia.org/wiki/Linear_recurrence_with_constant_coefficients#Order_2

https://en.wikipedia.org/wiki/Linear_recurrence_with_constant_coefficients#General_solution

The sequences are known under multiple different names for specific values of coefficients p, q, r ...

https://en.wikipedia.org/wiki/Lucas_sequence (second order)

https://en.wikipedia.org/wiki/Lucas_sequence#Specific_names   (second order)

https://en.wikipedia.org/wiki/Frobenius_pseudoprime  (second order)

https://en.wikipedia.org/wiki/Perrin_number  (third order)

https://en.wikipedia.org/wiki/Padovan_sequence (third order)

https://en.wikipedia.org/wiki/Supergolden_ratio#Narayana_sequence  (third order)


The terms S[n] are not very difficult to compute. Similar to exponentiation, there is an efficient "double and shift" addition chain in O(log2(n)) steps, making a general complexity close to O(log2(n) (order^2)) large number multiplications. When FFT techniques are involved, and taking into account the sparse matrix involved, complexity of each step is measured as O(order(3 + order) / 2) direct and inverse FFTs.

Special values for coefficients

# Primality test

https://en.wikipedia.org/wiki/Perrin_number#Perrin_primality_test

Adams and Shanks noted that primes also satisfy the congruence S(−p) ≡ S(−1) ≡ −1 (mod p).
If S(−n) = −1 and S(n) = 0 then n is a probable prime, that is: actually prime or a restricted Perrin pseudoprime.

??? noticed that a number of pseudoprimes can be removed by checking S[-i * n] == S[-i] and S[i * n] == S[i] for small values of i , e.g. 1 <= i <= 3 . This can be done with a low additional cost.

Going further, r^n mod n can be deduced for S[n+1] and S[n+2], this allow an additional Fermat test (ref ????)

Some special values and combinations of coefficients p,q,r,s ... help to speed-up the computations.

The density of pseudoprimes produced by a sequence of any order is not the same across sequences, depending on starting values and coefficients of the sequence.

# Pseudoprimes hunt in primality tests

https://ntheory.org/pseudoprimes.html

https://arxiv.org/pdf/1706.01265

# Deterministic primality test ?

There is an (unproven) assumption than lucas pseudoprime and Fermat pseudoprimes base-2 are exclusive sets of numbers. No counter-example has been found yet. It is conjectured the intersection of these 2 sets is infinitely large.




