#include <mutex>
#include "bqspc.h"

namespace bqspc {

/* These constants define the largest values the generator will allow for
 * the respective parameters. */
const static int Max_qScalarsDegree1 = 2;
const static int Max_qScalarsDegree2Pure = 2;
const static int Max_qScalarsDegree2Mixed = 2;
const static int Max_qPS_power = 2;
const static int Max_qPS_dilation1 = 2;
const static int Max_qPS_dilation2 = 2;
const static int Max_qPS_subScalars = 2;

/* Advances the state of the generator. */
void ParameterGenerator::advance(void)
{
	/* First the function $c(n_0, \dots, n_\ell)$ is generated. */
	for (int index = 0; index < this->indicesInUse; ++index) {
		this->qScalarsDegree1[index]++;

		if (this->qScalarsDegree1[index] <= Max_qScalarsDegree1) {
			return;
		}
		
		this->qScalarsDegree1[index] = 0;
	}

	for (int index = 0; index < this->indicesInUse; ++index) {
		this->qScalarsDegree2Pure[index]++;

		if (this->qScalarsDegree2Pure[index] <= Max_qScalarsDegree2Pure) {
			return;
		}
		
		this->qScalarsDegree2Pure[index] = 0;
	}

	for (int index = 0; index < this->indicesInUse
		 * (this->indicesInUse - 1) / 2; ++index) {

		this->qScalarsDegree2Mixed[index]++;

		if (this->qScalarsDegree2Mixed[index] <= Max_qScalarsDegree2Mixed) {
			return;
		}
		
		this->qScalarsDegree2Mixed[index] = 0;
	}

	/* Reset to the default state when all combinations are exhausted. */
	this->qScalarsDegree1[0] = 1;

	/* The parameters of the $q$-Pochhammer symbols. */
	for (int nIndex = 0; nIndex < this->qPSInUse; ++nIndex) {
		this->qPS[nIndex].dilation1++;

		if (this->qPS[nIndex].dilation1 <= Max_qPS_dilation1) {
			return;
		}

		this->qPS[nIndex].dilation1 = 1;
		this->qPS[nIndex].dilation2++;

		if (this->qPS[nIndex].dilation2 <= Max_qPS_dilation2) {
			return;
		}

		this->qPS[nIndex].dilation2 = 1;
		this->qPS[nIndex].power++;

		/* Do not allow a power of 0. */
		if (this->qPS[nIndex].power == 0) this->qPS[nIndex].power++;

		if (this->qPS[nIndex].power < 0) {
			this->qPS[nIndex].negativePrefix = false;
		} else {
			this->qPS[nIndex].negativePrefix = true;
		}

		if (this->qPS[nIndex].power <= Max_qPS_power) {
			return;
		}

		this->qPS[nIndex].power = -Max_qPS_power;

		/* Generate the function $s_i(n_0, \dots n_\ell)$. */
		for (int kIndex = 0; kIndex < this->indicesInUse; ++kIndex) {
			this->qPS[nIndex].subScalars[kIndex]++;

			if (this->qPS[nIndex].subScalars[kIndex] <= Max_qPS_subScalars) {
				return;
			}

			this->qPS[nIndex].subScalars[kIndex] = 0;
		}

		this->qPS[nIndex].subScalars[0] = 1;
	}

	/* Number of $q$-Pochhammer symbols. */
	this->qPSInUse++;

	if (this->qPSInUse <= MaxQPS) return;

	this->qPSInUse = 0;

	/* Number of indices. */
	this->indicesInUse++;

	if (this->indicesInUse <= MaxIndices) return;

	/* When this is reached, we have exhausted every parameter combination. */
	this->continueWorking = false;
}

/* Must be held whenever the state of the generator is accessed in any way. */
static std::mutex generatorLock;

/* Gives the worker thread some jobs. */
void ParameterGenerator::populate(WorkerThread& worker)
{
	/* This is the only place where the mutex needs to be held. */
	std::scoped_lock<std::mutex> lock(generatorLock);

	worker.jobQueueLength = 0;

	while (this->continueWorking) {
		*worker.jobQueue[worker.jobQueueLength++] = *this;
		if (worker.jobQueueLength == JobQueueLimit) break;

		this->advance();
	}
}

/* Sets all the parameters to their initial state. */
ParameterGenerator::ParameterGenerator(void)
{
	this->continueWorking = true;
	this->alternatingSign = false;
	this->dividePowerBy2 = false;
	this->indicesInUse = 2;
	this->qPSInUse = 0;

	for (int index = 0; index < MaxIndices; ++index) {
		this->qScalarsDegree1[index] = 0;
		this->qScalarsDegree2Pure[index] = 0;
	}

	for (int index = 0; index < MaxIndices * (MaxIndices - 1) / 2; ++index) {
		this->qScalarsDegree2Mixed[index] = 0;
	}

	this->qScalarsDegree1[0] = 1;

	for (int nIndex = 0; nIndex < MaxQPS; ++nIndex) {
		this->qPS[nIndex].dilation1 = 1;
		this->qPS[nIndex].dilation2 = 1;
		this->qPS[nIndex].negativePrefix = false;
		this->qPS[nIndex].power = -Max_qPS_power;

		for (int kIndex = 1; kIndex < MaxIndices; ++kIndex) {
			this->qPS[nIndex].subScalars[kIndex] = 0;
		}

		this->qPS[nIndex].subScalars[0] = 1;
	}
}

};

