#include <iostream>
#include <mutex>
#include <string>
#include <sstream>
#include "bqspc.h"

namespace bqspc
{

/* Prints out a conjectured sum-product identity fully formatted for LaTeX. */
void WorkerThread::reportIdentity(Parameters& parameters,
								  ProductSignature &signature)
{
	std::string sum;
	std::string prod;
	std::string sumNum;
	std::string sumDen;
	std::string prodNum;
	std::string prodDen;
	std::stringstream output;

	bool powerUsePlus = false;

	/* Helper function for printing powers of $q$. */
	auto prettyPrint = [&](int value) {
		if (value == 0) return std::string("1");

		if (value == 1) return std::string("q");

		return "q^{" + std::to_string(value) + "}";
	};

	/* Sigma notation. */
	if (parameters.indicesInUse == 1) {
		sum += "\\sum_{n_0\\geq 0} ";
	} else if (parameters.indicesInUse == 2) {
		sum += "\\sum_{n_0,n_1\\geq 0} ";
	} else {
		sum += "\\sum_{n_0,\\dots,n_{" + std::to_string(
			   parameters.indicesInUse - 1) + "}\\geq 0} ";
	}

	sumNum += "q^{";

	/* The function $c(n_0, \dots, n_\ell)$. */
	for (int index = 0; index < parameters.indicesInUse; ++index) {
		std::string term;

		if (parameters.qScalarsDegree2Pure[index] == 0) {
			continue;
		} else if (parameters.qScalarsDegree2Pure[index] != 1) {
			term += std::to_string(parameters.qScalarsDegree2Pure[index]);
		}

		term += "n_{" + std::to_string(index) + "}^2";

		if (powerUsePlus) {
			sumNum += "+" + term;
		} else {
			sumNum += term;
			powerUsePlus = true;
		}
	}

	for (int nIndex = 0; nIndex < parameters.indicesInUse; ++nIndex) {
		for (int kIndex = nIndex + 1; kIndex <
			 parameters.indicesInUse; ++kIndex) {

			std::string term;

			int mIndex = nIndex * (parameters.indicesInUse - 1)
					   - nIndex * (nIndex - 1) / 2 + kIndex - 1;

			if (parameters.qScalarsDegree2Mixed[mIndex] == 0) {
				continue;
			} else if (parameters.qScalarsDegree2Mixed[mIndex] != 1) {
				term += std::to_string(
						parameters.qScalarsDegree2Mixed[mIndex]);
			}

			term += "n_{" + std::to_string(nIndex) + "}"
			      + "n_{" + std::to_string(kIndex) + "}";

			if (powerUsePlus) {
				sumNum += "+" + term;
			} else {
				sumNum += term;
				powerUsePlus = true;
			}
		}
	}

	for (int index = 0; index < parameters.indicesInUse; ++index) {
		std::string term;

		if (parameters.qScalarsDegree1[index] == 0) {
			continue;
		} else if (parameters.qScalarsDegree1[index] != 1) {
			term += std::to_string(parameters.qScalarsDegree1[index]);
		}

		term += "n_{" + std::to_string(index) + "}";

		if (powerUsePlus) {
			sumNum += "+" + term;
		} else {
			sumNum += term;
			powerUsePlus = true;
		}
	}

	sumNum += "}";

	/* The $q$-Pochhammer symbols. */
	for (int nIndex = 0; nIndex < parameters.qPSInUse; ++nIndex) {
		int powerAbs = parameters.qPS[nIndex].power;
		std::string qPS;

		qPS += "(";

		if (parameters.qPS[nIndex].negativePrefix) {
			qPS += "-";
		}

		if (parameters.qPS[nIndex].dilation1 == 0) {
			qPS += "1";
		} else {
			qPS += prettyPrint(parameters.qPS[nIndex].dilation1);
		}

		qPS += "; " + prettyPrint(parameters.qPS[nIndex].dilation1) + ")_{";

		powerUsePlus = false;

		for (int kIndex = 0; kIndex < parameters.indicesInUse; ++kIndex) {
			std::string term;

			if (parameters.qPS[nIndex].subScalars[kIndex] == 0) {
				continue;
			} else if (parameters.qPS[nIndex].subScalars[kIndex] != 1) {
				term = std::to_string(
					   parameters.qPS[nIndex].subScalars[kIndex]);
			}

			term += "n_{" + std::to_string(kIndex) + "}";

			if (powerUsePlus) {
				qPS += "+" + term;
			} else {
				qPS += term;
				powerUsePlus = true;
			}
		}

		qPS += "}";

		if (powerAbs < 0) {
			powerAbs = -powerAbs;
		}

		if (powerAbs > 1) {
			qPS += "^{" + std::to_string(powerAbs) + "}";
		}

		if (parameters.qPS[nIndex].power > 0) {
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

	/* Now for the product side. */
	for (int index = 0; index < signature.period; ++index) {
		int powerAbs;
		std::string qPS;

		if (signature.powers[index] == 0) continue;

		qPS += "(" + prettyPrint(index + 1) + "; "
		    + prettyPrint(signature.period) + ")_{\\infty}";
		powerAbs = signature.powers[index];

		if (powerAbs < 0) {
			powerAbs = -powerAbs;
		}

		if (powerAbs != 1) {
			qPS += "^{" + std::to_string(powerAbs) + "}";
		}

		if (signature.powers[index] > 0) {
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

	/* Write the result to stdout. This is threadsafe and does not require
	 * holding a mutex. */
	output << "\\begin{equation}\n" + sum + " = "
			  + prod + "\n\\end{equation}\n";
	std::cout << output.str();
}

/* Attempts to find a sum-product identity from the given parameters. */
void WorkerThread::tryCombination(Parameters& parameters)
{
	ProductSignature signature;
	QSeries candidate;

	/* Generate the $q$-series coefficients and factor them. */
	candidate.qSeries(parameters);
	signature.factorize(candidate);

	/* If there is no sum-product identity found or if the identity is dilated
	 * then this parameter combination is considered a failure. */
	if (signature.period == 0 || signature.dilation() > 1) return;

	/* Otherwise, report the identity and move on. */
	this->reportIdentity(parameters, signature);
}

/* Acquires and executes jobs from the generator on loop. */
void WorkerThread::jobLoop(void)
{
	for (;;) {

		/* Get some work. */
		this->generator->populate(*this);

		/* The generator will notify the worker threads that
		 * the work is finished by not providing any jobs here. */
		if (this->jobQueueLength == 0) return;

		for (int index = 0; index < this->jobQueueLength; ++index) {
			this->tryCombination(*this->jobQueue[index]);
		}
	}
}

};

