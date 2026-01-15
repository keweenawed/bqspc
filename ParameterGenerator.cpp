#include <mutex>
#include <thread>
#include "bqspc.h"

namespace bqspc
{

/* Must be held by any thread while acquiring jobs. */
static std::mutex generatorLock;

/* Some checks to avoid redundancy. Could be made completely effective with
 * much more work, but there are very few duplicated answers as it is now. */
void ParameterGenerator::checkParameterRedundancy(void)
{
	if (this->qPSLength <= 1) return;

	for (int n = 1; n < this->qPSLength; ++n) {

		/* Insist that the powers on the q-Pochhammer symbols are
		 * weakly decreasing from left to right. */
		if (this->qPS[n * 6 + 3] > this->qPS[n * 6 - 3]) {
			this->advanceState();
			return;
		}

		/* If the powers on two neighboring q-Pochhammer symbols are 
		 * equal, order the d3 in the same way next. */
		if (this->qPS[n * 6 + 3] == this->qPS[n * 6 - 3]) {
			if (this->qPS[n * 6 + 2] > this->qPS[n * 6 - 4]) {
				this->advanceState();
				return;
			}
		}
	}
}

/* Called to get the next combination of parameters. */
void ParameterGenerator::advanceState(void)
{
	++this->qScalarDeg1;

	if (this->qScalarDeg1 <= MaxQPowerDegree1) return;
		
	this->qScalarDeg1 = 0;
	++this->qScalarDeg2;

	if (this->qScalarDeg2 <= MaxQPowerDegree2) return;

	this->qScalarDeg2 = 1;

	/* Go through all the q-Pochhammer symbol parameters. */
	for (int n = 0; n < this->qPSLength * 6; ++n) {

		/* Ignore anything to do with z here. */
		if (n % 6 == 0) continue;

		this->qPS[n]++;

		switch (n % 6) {
		case 1:
		case 2:
		case 4:
			if (qPS[n] <= MaxQPSParameters[n % 6]) return;

			qPS[n] = 1;
			continue;

		case 3:
			/* Disallow powers of 0 to avoid duplication. */
			if (qPS[n] == 0) qPS[n] = 1;

			if (qPS[n] <= MaxQPSParameters[n % 6]) return;

			qPS[n] = -MaxQPSParameters[n % 5];
			continue;

		case 5:
			if (qPS[n] <= MaxQPSParameters[n % 6]) return;

			qPS[n] = 0;
			continue;
		}
	}

	++this->qPSLength;

	/* We are finished when more q-Pochhammer symbols than are permitted
	 * would be required to advance the state. */
	if (this->qPSLength > MaxNumberQPS) {
		this->keepWorking = false;
		return;
	}

	/* Initialize the newest q-Pochhammer symbol. */
	this->qPS[6 * (qPSLength - 1) + 0] = 0;
	this->qPS[6 * (qPSLength - 1) + 1] = 1;
	this->qPS[6 * (qPSLength - 1) + 2] = 1;
	this->qPS[6 * (qPSLength - 1) + 3] = -MaxQPSParameters[3];
	this->qPS[6 * (qPSLength - 1) + 4] = 1;
	this->qPS[6 * (qPSLength - 1) + 5] = 0;
}

/* Fills the job cache of the given worker thread. */
void ParameterGenerator::populateJobQueue(WorkerThread& worker)
{
	/* This is the only place where the mutex needs to be held. */
	std::scoped_lock<std::mutex> lock(generatorLock);
	worker.jobQueueLength = 0;

	while (this->keepWorking) {
		this->checkParameterRedundancy();
		worker.jobQueue[worker.jobQueueLength++] = *this;

		if (worker.jobQueueLength == MaxJobQueueLimit) break;

		this->advanceState();
	}
}

};

