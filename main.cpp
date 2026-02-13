#include <iostream>
#include <thread>
#include "bqspc.h"

namespace bqspc {

/* The number of worker threads to use. */
const static int WorkerThreadsToUse = 10;

extern long *precomputedDivisorList[MaxSeriesLimit];
extern int precomputedDivisorFunction[MaxSeriesLimit];

};

using namespace bqspc;

static void workerThreadEntry(ParameterGenerator *generator)
{
	WorkerThread worker(generator);

	worker.jobLoop();
}

/* No arguments are parsed. The range of parameters used must be specified at
 * compile time for now. The result will be written in the format of a LaTeX
 * file that can immediately be built into a pdf without extra work. */
int main(void)
{
	ParameterGenerator generator;
	std::thread threads[WorkerThreadsToUse];

	/* Populate the precomputed divisors using a brute force algorithm since
	 * this is only done once, and does not benefit from extra efficiency. */
	precomputedDivisorFunction[0] = 0;
	precomputedDivisorList[0] = nullptr;

	for (long nIndex = 1; nIndex < MaxSeriesLimit; ++nIndex) {
		long buffer[MaxSeriesLimit];
		int length = 0;

		/* Fill the buffer with every value that divides nIndex. */
		for (long kIndex = 1; kIndex <= nIndex / 2; ++kIndex) {
			if (nIndex % kIndex == 0) {
				buffer[length++] = kIndex;
			}
		}

		buffer[length++] = nIndex;

		/* Create the next entry in the divisor list and copy the divisors
		 * over, and save the value of the divisor function. */
		precomputedDivisorList[nIndex] = new long[length];

		for (int kIndex = 0; kIndex < length; ++kIndex) {
			precomputedDivisorList[nIndex][kIndex] = buffer[kIndex];
		}

		precomputedDivisorFunction[nIndex] = length;
	}

	/* Header for the LaTeX output. */
	std::cout << "\\documentclass{article}\n"\
				 "\\usepackage[margin=1in]{geometry}\n"\
				 "\\begin{document}\n\n";

	/* Create the worker threads. */
	for (int index = 0; index < WorkerThreadsToUse; ++index) {
		threads[index] = std::thread(workerThreadEntry, &generator);
	}

	/* Cleanup. */
	for (int index = 0; index < WorkerThreadsToUse; ++index) {
		threads[index].join();
	}

	for (int index = 1; index < MaxSeriesLimit; ++index) {
		delete precomputedDivisorList[index];
	}

	/* Footer for the LaTeX output. */
	std::cout << "\\end{document}\n";

	return 0;
}

