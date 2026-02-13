#include "bqspc.h"

namespace bqspc {

/* Computes the truncated reciprocal of the $q$-series. This method assumes
 * that the constant coefficient equals 1 since fractional coefficients are
 * not supported. */
void QSeries::reciprocal(void)
{
	QSeries copy = *this;

	/* If $1/(\sum_{n \geq 0} a_nq^n) = \sum_{n \geq 0} b_nq^n$, then we have
	 * $b_n = -\sum_{k=1}^n a_{k} b_{n-k}$ for $n \geq 1$. This algorithm
	 * computes $a_n$ using dynamic programming and this relation. The time
	 * complexity is quadratic in the number of coefficients. */
	for (int nIndex = 1; nIndex < this->limit; ++nIndex) {
		this->coefficients[nIndex] = 0;

		for (int kIndex = 1; kIndex <= nIndex; ++kIndex) {
			this->coefficients[nIndex] -= this->coefficients[nIndex - kIndex]
										* copy.coefficients[kIndex];
		}
	}
}

/* The truncated $q$-series raised to the given power. For negative powers,
 * the constant coefficient must equal 1 as the above function is called. */
void QSeries::raiseToPower(int power)
{
	if (power == 0) {
		this->zero();
		this->coefficients[0] = 1;
		return;
	}

	if (power < 0) {
		this->reciprocal();
		power = -power;
	}

	if (power == 1) {
		return;
	}

	if (power == 2) {
		QSeries copy = *this;
		*this *= copy;
		return;
	}

	/* For larger powers than 2 in magnitude, we use a classic algorithm with
	 * logarithmic time complexity in the power. This breaks the problem into
	 * recursively computing half the power, then squaring that, separated
	 * into an even and odd case. */
	if (power % 2 == 0) {
		this->raiseToPower(power / 2);
		this->raiseToPower(2);
	} else {
		QSeries copy = *this;

		this->raiseToPower(power / 2);
		this->raiseToPower(2);
		*this *= copy;
	}
}

/* The truncated finite $q$-Pochhammer symbol of the form
 * $(\pm q^{dilation1}; q^{dilation2})_{subscript}$ where $\pm$ is
 * negative when the value of negativePrefix is set to true. The definition
 * of the finite $q$-Pochhammer symbol in general is the polynomial
 * $(z;q)_n = \prod_{k=0}^{n-1}(1-zq^k)$. */
void QSeries::qPochhammer(int dilation1, int dilation2, bool negativePrefix,
						  int subscript)
{
	QSeries shiftedCopy(this->limit);

	/* We start by setting the result to one, and then repeatedly multiply by
	 * the polynomial $1 \pm q^{shift}$. This is sped up by avoiding the full
	 * multiplication algorithm, and so uses quadratic time complexity in the
	 * subscript value times number of coefficients. */
	this->zero();
	this->coefficients[0] = 1;

	if (subscript == 0) return;

	/* Loop over every factor in the $q$-Pochhammer symbol. */
	for (int nIndex = 0; nIndex < subscript; ++nIndex) {
		int shift = dilation1 + nIndex * dilation2;

		/* As a function of nIndex, shift is weakly increasing and so we are
		 * done the first time shift grows this large. */
		if (shift >= this->limit) break;

		/* Copy the coefficients, but shifted to the right to mimic a
		 * multiplication by $q^{shift}$. */
		if (negativePrefix) {
			for (int kIndex = 0; kIndex < limit - shift; ++kIndex) {
				shiftedCopy.coefficients[kIndex + shift] =
					this->coefficients[kIndex];
			}
		} else {
			for (int kIndex = 0; kIndex < limit - shift; ++kIndex) {
				shiftedCopy.coefficients[kIndex + shift]
					= -this->coefficients[kIndex];
			}
		}

		/* Now, add the shifted coefficients to the original. */
		for (int kIndex = shift; kIndex < this->limit; ++kIndex) {
			this->coefficients[kIndex] += shiftedCopy.coefficients[kIndex];
		}
	}
}

/* The truncated $q$-binomial coefficient. Both top and bottom must be
 * non-negative. If top is smaller than, bottom the result is zero, and is
 * otherwise computed as $(q;q_{top}/((q;q)_{bottom}(q;q)_{top - bottom})$.*/
void QSeries::qBinomial(int top, int bottom)
{
	if (bottom > top) {
		this->zero();
	} else if (bottom == top || bottom == 0) {
		this->zero();
		this->coefficients[0] = 1;
	} else {
		int bigSubscript;
		QSeries divisor(this->limit);

		/* Find out which subscript on the denominator is the biggest. The
		 * associated $q$-Pochhammer symbol will be cancelled from the
		 * numerator to start with. */
		bigSubscript = (bottom > top - bottom) ? bottom : top - bottom;
		this->qPochhammer(bigSubscript + 1, 1, true, top - bigSubscript);

		/* Perform division by multiplying by the reciprocal. */
		divisor.qPochhammer(1, 1, true, top - bigSubscript);
		divisor.reciprocal();
		*this *= divisor;
	}
}

/* Computes the truncated coefficients of a particular term in a $q$-series,
 * determined by the combination of indices. The calling function handles the
 * shift by $q^{c(n_0, \dots, n_\ell)}$. */
void QSeries::qSeriesTerm(Parameters& parameters, int (&indices)[MaxIndices])
{	
	this->zero();
	this->coefficients[0] = 1;

	/* Multiply the term by all the $q$-Pochhammer symbols. */
	for (int qPSIndex = 0; qPSIndex < parameters.qPSInUse; ++qPSIndex) {
		QSeries factor(this->limit);
		int subscript = 0;

		/* Compute $s_i(n_0, \dots, n_\ell)$. */
		for (int index = 0; index < parameters.indicesInUse; ++index) {
			subscript += parameters.qPS[qPSIndex].subScalars[index]
					   * indices[index];
		}

		factor.qPochhammer(parameters.qPS[qPSIndex].dilation1,
						   parameters.qPS[qPSIndex].dilation2,
						   parameters.qPS[qPSIndex].negativePrefix,
						   subscript);
		factor.raiseToPower(parameters.qPS[qPSIndex].power);
		*this *= factor;
	}

	if (parameters.alternatingSign) {
		int indexSum = 0;

		for (int index = 0; index < parameters.indicesInUse; ++index) {
			indexSum += indices[index];
		}

		/* If the $q$-series indices sum to an odd value, negate the term. */
		if (indexSum % 2 == 1) {
			*this = -(*this);
		}
	}
}

/* Computes the value of of the power of $q$ $c(n_0, \dots, n_\ell)$ for a
 * $q$-series. */
int QSeries::qSeriesPower(Parameters& parameters, int (&indices)[MaxIndices])
{
	int power = 0;

	for (int nIndex = 0; nIndex < parameters.indicesInUse; ++nIndex) {

		/* The pure degree 1 and 2 contributions. */
		power += parameters.qScalarsDegree1[nIndex] * indices[nIndex]
			   + parameters.qScalarsDegree2Pure[nIndex]
			   * indices[nIndex] * indices[nIndex];

		/* The mixed degree 2 contributions. */
		for (int kIndex = nIndex + 1; kIndex <
			 parameters.indicesInUse; ++kIndex) {

			/* See bqspc.h for an explanation of this indexing. */
			int mIndex = nIndex * (parameters.indicesInUse - 1)
					   - nIndex * (nIndex + 1) / 2 + kIndex - 1;

			power += parameters.qScalarsDegree2Mixed[mIndex]
				   * indices[nIndex] * indices[kIndex];
		}
	}

	if (parameters.dividePowerBy2) {
		power /= 2;
	}

	return power;
}

/* Computes the truncated $q$-series coefficients determined by the given
 * parameters. */
void QSeries::qSeries(Parameters& parameters)
{
	int indices[MaxIndices];
	int index;

	/* Start with the $q$-series summation indices $n_0, \dots, n_\ell$
	 * all equal to zero. */
	for (index = 0; index < parameters.indicesInUse; ++index) {
		indices[index] = 0;
	}

	/* The term given by the indices all equaling zero is identically 1 for
	 * any combination of parameters. */
	this->zero();
	this->coefficients[0] = 1;

	/* Iterate through the other possible combinations of indices. */
	for (;;) {
		int power;

		for (index = 0; index < parameters.indicesInUse; ++index) {
			++indices[index];

			if (indices[index] < MaxSeriesLimit) {
				break;
			}

			indices[index] = 0;
		}

		/* This condition is met only when we are finished. */
		if (index == parameters.indicesInUse) {
			return;
		}

		power = this->qSeriesPower(parameters, indices);

		/* Compute the term and add the contribution if there are any
		 * coefficients that will not be truncated. */
		if (power < this->limit) {
			QSeries term(this->limit - power);

			term.qSeriesTerm(parameters, indices);
			term.limit = this->limit;
			term.translate(power);
			*this += term;
		}
	}
}

};

