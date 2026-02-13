#include "bqspc.h"

namespace bqspc {

/* Computes the greatest common divisor using the Euclidean algorithm. */
long ProductSignature::pairwiseGCD(long value1, long value2)
{
	if (value1 < 0) {
		value1 = -value1;
	}

	if (value2 < 0) {
		value2 = -value2;
	}

	if (value2 == 0) {
		return value1;
	}

	return pairwiseGCD(value2, value1 % value2);
}

/* Returns the factor $n \geq 1$ to which the product signature is dilated, or
 * 1 if the product is not dilated. If the product is dilated, that means the
 * product can be expressed $f(q^n)$ where $f(q)$ is a power series in $q$. */
long ProductSignature::dilation(void)
{
	long gcd = this->powers[0];

	/* This simply computes the GCD of every power in the pattern. */
	for (int index = 1; index < this->period; ++index) {
		if (this->powers[index] == 0) {
			continue;
		}

		gcd = pairwiseGCD(index + 1, gcd);

		if (gcd == 1) {
			break;
		}
	}

	return gcd;
}

/* The entry precomputedDivisorList[index] is an array containing every
 * divisor of the value index. */
long *precomputedDivisorList[MaxSeriesLimit];

/* The entry precomputedDivisorFunction[index] is the number of divisors of
 * the value index. */
int precomputedDivisorFunction[MaxSeriesLimit];

/* Factorizes the $q$-series as a product of geometric series of the form
 * $\prod_{n \geq 1} \frac{1}{(1-q^n)^{a_n}}$ using an algorithm found in
 * George Andrews's book The Theory of Partitions, guaranteeing equality for
 * all coefficients before series.limit. The constant coefficient must equal
 * 1 for this to behave correctly. */
void ProductSignature::factorize(QSeries& series)
{
	long powers[series.limit - 1];

	/* Compute the geometric series powers using dynamic programming. Let
	 * $r(n)$ be the $n$th q-series coefficient. From the relationship
	 * $\sum_{n\geq 0} r(n)q^n = \prod_{n\geq 1}\frac{1}{(1-q^n)^{a_n}}$
	 * we can derive using logarithmic differentiation the recurrence
	 * $a_n = r(n) - \frac{1}n \sum_{k=1}^{n}r(n-k)
	 * \times \sum_{d\mid k, d \neq n}d a_d$. This loop computes the values of
	 * $a_n$ and stores them in powers[n - 1]. */
	for (int nIndex = 1; nIndex < series.limit; ++nIndex) {
		long *divisors;
		int length;
		long power = 0;

		for (int kIndex = 1; kIndex < nIndex; ++kIndex) {
			divisors = precomputedDivisorList[kIndex];
			length = precomputedDivisorFunction[kIndex];

			for (int dIndex = 0; dIndex < length; ++dIndex) {
				power -= series.coefficients[nIndex - kIndex]
				      * divisors[dIndex]
				      * powers[divisors[dIndex] - 1];
			}
		}

		divisors = precomputedDivisorList[nIndex];
		length = precomputedDivisorFunction[nIndex];

		/* This handles the special case excluded from the above loop where
		 * kIndex == nIndex. */
		for (int dIndex = 0; dIndex < length - 1; ++dIndex) {
			power -= divisors[dIndex] * powers[divisors[dIndex] - 1];
		}

		power /= nIndex;
		power += series.coefficients[nIndex];
		powers[nIndex - 1] = power;
	}

	/* The factorization algorithm is finished now. Proceed to look for a
	 * repeating pattern. The period is set to 0 if none is found. */
	for (this->period = 1; this->period <= MaxProductSignatureLength;
	     this->period++) {

		bool success = true;

		for (int index = period; index < series.limit - 1; ++index) {
			if (powers[index] != powers[index % period]) {
				success = false;
				break;
			}
		}

		if (!success) continue;

		/* The first pattern to work is minimal. */
		for (int index = 0; index < period; ++index) {
			this->powers[index] = powers[index];
		}

		return;
	}

	/* No period worked if this is reached. */
	this->period = 0;
	return;
}

};

