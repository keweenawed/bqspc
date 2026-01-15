#include "bqspc.h"

namespace bqspc {

/* Computes the positive greatest common divisor of the arguments using the
 * Euclidean algorithm. */
static long gcdPair(long value1, long value2)
{
	if (value2 < 0) value2 = -value2;

	if (value1 < 0) value1 = -value1;

	if (value2 == 0) return value1;

	return gcdPair(value2, value1 % value2);
}

/* Returns the value to which the product is dilated. If that value is 1,
 * there is no dilation. */
long ProductSignature::dilation(void)
{
	long gcd = this->period;

	for (int n = 0; n < this->period; ++n) {
		if (this->powers[n] == 0) continue;

		gcd = gcdPair(n + 1, gcd);

		if (gcd == 1) break;
	}

	return gcd;
}

};

