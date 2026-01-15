#include <iostream>
#include <string>
#include <sstream>
#include "bqspc.h"

namespace bqspc
{

/* Prints out a conjectured identity formatted for LaTeX. */
void WorkerThread::reportIdentity(class Parameters& parameters,
				  class ProductSignature& signature)
{
	bool bivariate = false;
	std::string sum;
	std::string prod;
	std::string sumNum;
	std::string sumDen;
	std::string prodNum;
	std::string prodDen;
	std::stringstream output;

	sum += "\\sum_{n=0}^\\infty ";

	auto prettyPrint = [&](int value, std::string variable) {
		if (value == 0) return std::string("1");

		if (value == 1) return variable;

		return variable
			+ "^{"
			+ std::to_string(value)
			+ "}";
	};

	if (parameters.alternatingSign) {
		sumNum += "(-1)^n";
	}

	if (parameters.zScalar != 0) {
		sumNum += "z^{";
		bivariate = true;

		if (parameters.zScalar != 1) {
			sumNum += std::to_string(parameters.zScalar);
		}

		sumNum += "n}";
	}

	sumNum += "q^{";

	if (parameters.dividePowerBy2) {
		sumNum += "\\frac{";
	}

	if (parameters.qScalarDeg2 != 0) {
		if (parameters.qScalarDeg2 != 1) {
			sumNum += std::to_string(parameters.qScalarDeg2);
		}

		sumNum += "n^2";
	}
	
	if (parameters.qScalarDeg1 != 0) {
		if (parameters.qScalarDeg2 != 0) {
			sumNum += " + ";
		}

		if (parameters.qScalarDeg1 != 1) {
			sumNum += std::to_string(parameters.qScalarDeg1);
		}

		sumNum += "n";
	}

	if (parameters.dividePowerBy2) {
		sumNum += "}{2}";
	}

	sumNum += "}";

	for (int n = 0; n < parameters.qPSLength; ++n) {
		int power;
		std::string qPS;

		qPS += "(";

		if (parameters.qPS[6 * n + 3] > 0) {
			qPS += "-";
		}

		if (parameters.qPS[6 * n + 0] != 0) {
			bivariate = true;
			qPS += prettyPrint(parameters.qPS[6 * n + 0], "z");
		}

		if (parameters.qPS[6 * n + 1] != 0) {
			qPS += prettyPrint(parameters.qPS[6 * n + 1], "q");
		}

		if (parameters.qPS[6 * n + 0] == 0 &&
		    parameters.qPS[6 * n + 1] == 0) {
			qPS += "1";
		}
		qPS += "; " + prettyPrint(parameters.qPS[6 * n + 2], "q")
		    + ")_{";

		if (parameters.qPS[6 * n + 4] != 1) {
			qPS += std::to_string(parameters.qPS[6 * n + 4]);
		}

		qPS += "n";

		if (parameters.qPS[6 * n + 5] != 0) {
		    qPS += " + "
		        + std::to_string(parameters.qPS[6 * n + 5]);
		}

		qPS += "}";
		power = parameters.qPS[6 * n + 3];

		if (power < 0) {
			power = -power;
		}

		if (power > 1) {
			qPS += "^{" + std::to_string(power) + "}";
		}

		if (parameters.qPS[6 * n + 3] > 0) {
			sumNum += qPS;
		} else {
			sumDen += qPS;
		}
	}

	if (sumDen != "") {
		sum += "\\frac{" + sumNum
		       + "}{" + sumDen + "}";
	} else {
		sum += sumNum;
	}

	for (int n = 0; n < signature.period; ++n) {
		int power;
		std::string qPS;

		if (signature.powers[n] == 0) continue;

		qPS += "(";

		if (bivariate) {
			qPS += "z";
		}

		qPS += prettyPrint(n + 1, "q") + "; "
		    + prettyPrint(signature.period, "q") + ")_{\\infty}";

		power = signature.powers[n];

		if (power != 1) {
			if (power < 0) power = -power;

			qPS += "^{" + std::to_string(power) + "}";
		}

		if (signature.powers[n] > 0) {
			prodDen += qPS;
		} else {
			prodNum += qPS;
		}
	}

	if (signature.period == 1 && signature.powers[0] == 0) {
		prod = "1";
	} else {
		if (prodDen == "") {
			prod = prodDen;
		} else {
			if (prodNum == "") prodNum = "1";

			prod = "\\frac{" + prodNum + "}{" + prodDen + "}";
		}
	}

	output << "$$" + sum + " = " + prod + "$$\n";
	std::cout << output.str();
}

/* Determines if the given parameters may lead to a q-series identity. */
void WorkerThread::tryCombinationUv(Parameters& parameters)
{
	ProductSignature signature;
	SeriesUv candidate;

	/* Generate the univariate q-series coefficients and factor them. */
	candidate.qSeries(parameters);
	candidate.factorize(signature);

	if (signature.period == 0) return;

	/* Compute the GCD of all entries in the signature. If this is not
	 * equal to 1, we have a dilated result which should be thrown out. */
	if (signature.dilation() > 1) return;

	/* If a pattern was detected, we have a conjectured univariate
	 * identity. Report this. */
	this->reportIdentity(parameters, signature);
}

/* Acquires and executes jobs from the parameter generator on loop. */
void WorkerThread::jobLoop(class ParameterGenerator& generator)
{
	for (;;) {

		/* Get some work. */
		generator.populateJobQueue(*this);

		/* The parameter generator will notify the worker threads that
		 * the work is finished by not providing any jobs here. */
		if (this->jobQueueLength == 0) return;

		for (int n = 0; n < this->jobQueueLength; ++n) {
			this->tryCombinationUv(this->jobQueue[n]);

			/* If the power of q has odd coefficients on both the
			 * degree 1 and 2 terms, we can divide the total power
			 * by 2, which may lead an identity. */
			if ((this->jobQueue[n].qScalarDeg1 % 2) == 1 &&
			    (this->jobQueue[n].qScalarDeg2 % 2) == 1) {
				this->jobQueue[n].dividePowerBy2 = true;
				this->tryCombinationUv(this->jobQueue[n]);
			}
		}
	}
}

};

