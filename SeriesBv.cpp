#include "bqspc.h"

namespace bqspc
{

/* Computes the truncated bivariate q-Pochhammer symbol with form
 * $((prefix)z^{d1}q^{d2}; q^{d3})_{subscript}^{power}$. */
void SeriesBv::qPochhammer(int prefix, int d1, int d2, int d3,
			   int power, int subscript)
{
	/* There is no reason for this not to be on the stack other than
	 * the memory footprint possibly causing a stack overflow. */
	class SeriesBv *factor;

	this->zero();
	this->coefficients[0][0] = 1;

	/* The result is simply 1 in this case due to the truncation. */
	if (d1 >= this->zLimit) return;

	/* Start by ignoring the power. */
	if (power > 0) {
		factor = new SeriesBv();

		for (int k = 0; k < subscript; ++k) {
			int qp = d2 + k * d3;

			if (qp >= this->qLimit) break;

			/* Multiply by $1-(prefix)z^{d1}q^{qp}$. */
			factor->qLimit = this->qLimit;
			factor->zLimit = this->zLimit;
			factor->zero();
			factor->coefficients[0][0] = 1;
			factor->coefficients[qp][d1] += -prefix;
			*this *= *factor;
		}
	} else if (power < 0) {
		factor = new SeriesBv();

		for (int k = 0; k < subscript; ++k) {
			int qp = d2 + k * d3;

			if (qp >= this->qLimit) break;

			factor->qLimit = this->qLimit;
			factor->zLimit = this->zLimit;
			factor->zero();

			/* Multiply by $1/(1-(prefix)z^{d1}q^{qp})$ expressed
			 * $\sum_{m=0}^\infty ((prefix)z^{d1}q^{qp})^m$. */
			for (int m = 0, pp = 1; d1 * m < this->zLimit &&
			     qp * m < this->qLimit; ++m, pp *= prefix) {
				factor->coefficients[qp * m][d1 * m] += pp;
			}

			*this *= *factor;
		}
	}

	/* Now handle the power by repeatedly multiplying the result by
	 * itself. */
	if (power != 0) {
		*factor = *this;
		power = (power < 0) ? -power : power;

		while (power != 1) {
			*this *= *factor;
			--power;
		}

		delete factor;
	}
}

/* Finds the truncated bivariate q-series coefficients determined by the given
 * parameters. */
void SeriesBv::qSeries(class Parameters& parameters)
{

	class SeriesBv *term = new SeriesBv();
	class SeriesBv *factor = new SeriesBv();

	this->zero();

	/* Each instance of this loop builds one of the terms. The index
	 * n1 here is the summation index of the q-series. */
	for (int n1 = 0;; ++n1) {
		int qPower;
		int zPower;

		qPower = parameters.qScalarDeg1 * n1 +
			 parameters.qScalarDeg2 * n1 * n1;

		if (parameters.dividePowerBy2) qPower /= 2;

		/* We are finished when either the power of q or z exceeds
		 * their respective truncation limits. This assumes powers can
		 * only strictly increase. */
		if (qPower >= this->qLimit) break;

		zPower = parameters.zScalar * n1;

		if (zPower >= this->zLimit) break;

		/* Only compute coefficients that will not be thrown out. */
		term->qLimit = this->qLimit - qPower;
		term->zLimit = this->zLimit - zPower;
		term->zero();
		term->coefficients[0][0] = 1;

		/* Handle the q-Pochhammer symbols. */
		for (int n2 = 0; n2 < parameters.qPSLength; ++n2) {
			int prefix;

			if (parameters.qPS[6 * n2 + 3] > 0) {
				prefix = -1;
			} else {
				prefix = 1;
			}

			factor->qLimit = term->qLimit;
			factor->qPochhammer(prefix,
				parameters.qPS[6 * n2 + 0],
				parameters.qPS[6 * n2 + 1],
				parameters.qPS[6 * n2 + 2],
				parameters.qPS[6 * n2 + 3],
				parameters.qPS[6 * n2 + 4] * n1 +
				parameters.qPS[6 * n2 + 5]);
			*term *= *factor;
		}

		term->qLimit = this->qLimit;
		term->zLimit = this->zLimit;

		/* Multiply the term by the power of q and z. */
		term->translate(qPower, zPower);

		if (parameters.alternatingSign && (n1 % 2) == 1) {
			*term = -(*term);
		}
	
		*this += *term;
	}

	delete factor;
	delete term;
}

};
