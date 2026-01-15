#include "bqspc.h"

namespace bqspc
{

/* Table of divisors that are needed by SeriesUv::factorize. */
int precomputedDivisorValues[MaxSeriesLimit];
long *precomputedDivisorList[MaxSeriesLimit];

/* Computes the truncated univariate q-Pochhammer symbol with form
 * $((prefix)q^{d1}; q^{d2})_{subscript}^{power}$. See SeriesBv.cpp for more
 * documentation. */
void SeriesUv::qPochhammer(int prefix, int d1, int d2,
			   int power, int subscript)
{
	this->zero();
	this->coefficients[0] = 1;

	if (power > 0) {
		for (int k = 0; k < subscript; ++k) {
			class SeriesUv factor;
			int qp = d1 + k * d2;

			if (qp >= this->qLimit) break;

			factor.qLimit = this->qLimit;
			factor.zero();
			factor.coefficients[0] = 1;
			factor.coefficients[qp] += -prefix;
			*this *= factor;
		}
	} else if (power < 0) {
		for (int k = 0; k < subscript; ++k) {
			class SeriesUv factor;
			int qp = d1 + k * d2;

			if (qp >= this->qLimit) break;

			factor.qLimit = this->qLimit;
			factor.zero();

			for (int m = 0, pp = 1; qp * m < this->qLimit; ++m,
			     pp *= prefix) {
				factor.coefficients[qp * m] += pp;
			}

			*this *= factor;
		}
	}

	if (power != 0) {
		class SeriesUv copy;

		copy = *this;
		power = (power < 0) ? -power : power;

		while (power != 1) {
			*this *= copy;
			--power;
		}
	}
}

/* Finds the truncated univariate q-series coefficients determined by the
 * given parameters. */
void SeriesUv::qSeries(class Parameters& parameters)
{
	this->zero();

	for (int n1 = 0;; ++n1) {
		class SeriesUv term;
		int qPower;

		qPower = parameters.qScalarDeg1 * n1 +
			 parameters.qScalarDeg2 * n1 * n1;

		if (parameters.dividePowerBy2) qPower /= 2;

		if (qPower >= this->qLimit) break;

		term.qLimit = this->qLimit - qPower;
		term.zero();
		term.coefficients[0] = 1;

		for (int n2 = 0; n2 < parameters.qPSLength; ++n2) {
			class SeriesUv factor;
			int prefix;

			if (parameters.qPS[6 * n2 + 3] > 0) {
				prefix = -1;
			} else {
				prefix = 1;
			}

			factor.qLimit = term.qLimit;
			factor.qPochhammer(prefix,
				parameters.qPS[6 * n2 + 1],
				parameters.qPS[6 * n2 + 2],
				parameters.qPS[6 * n2 + 3],
				parameters.qPS[6 * n2 + 4] * n1 +
				parameters.qPS[6 * n2 + 5]);
			term *= factor;
		}

		term.qLimit = this->qLimit;
		term.translate(qPower);

		if (parameters.alternatingSign && (n1 % 2) == 1) {
			term = -term;
		}
	
		*this += term;
	}
}

/* Given truncated univariate q-series coefficients with a constant term equal
 * to 1, this method uniquely factorizes the series as a product of geometric 
 * series using an algorithm derived from George Andrews's book The Theory of
 * Partitions, so that equality holds up to the largest coefficient not
 * truncated. If this follows a pattern that suggests the series may be
 * written as a finite product of infinite q-Pochhammers symbols, the
 * so-called product signature holds this information. */
void SeriesUv::factorize(ProductSignature& signature)
{
	long powers[this->qLimit - 1];

	/* Recursively compute the geometric series powers. */
	for (int n = 1; n < this->qLimit; ++n) {
		long *divisors;
		long power = 0;
		int length;

		for (int k = 1; k < n; ++k) {
			divisors = precomputedDivisorList[k];
			length = precomputedDivisorValues[k];

			for (int m = 0; m < length; ++m) {
				power -= this->coefficients[n - k]
				      * divisors[m]
				      * powers[divisors[m] - 1];
			}
		}

		divisors = precomputedDivisorList[n];
		length = precomputedDivisorValues[n];

		for (int k = 0; k < length - 1; ++k) {
			power -= divisors[k] * powers[divisors[k] - 1];
		}

		power /= n;
		power += this->coefficients[n];
		powers[n - 1] = power;
	}

	signature.period = 0;

	/* Look for a repeating period in the powers. If none is found, the
	 * period is kept 0. */
	for (int period = 1; period <= MaxProductSignatureLength; ++period) {
		bool success = true;

		for (int n = period; n < this->qLimit - 1; ++n) {
			if (powers[n] != powers[n % period]) {
				success = false;
				break;
			}
		}

		if (!success) continue;

		/* If a pattern of minimal length is found, store the powers
		 * and finish. */
		for (int n = 0; n < period; ++n) {
			signature.powers[n] = powers[n];
		}

		signature.period = period;
		break;
	}

	return;
}

};

