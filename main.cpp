#include <iostream>
#include <thread>
#include "bqspc.h"

namespace bqspc {

extern int precomputedDivisorValues[MaxSeriesLimit];
extern long *precomputedDivisorList[MaxSeriesLimit];

};

using namespace bqspc;

/* Entry point for the worker threads. */
static void worker_thread_entry(ParameterGenerator& generator)
{
	/* This seems silly, but this is absolutely necessary to prevent
	 * stack overflow on some platforms. */
	WorkerThread *worker = new WorkerThread();

	worker->jobLoop(generator);
	delete worker;
}

int main(void)
{
	ParameterGenerator generator;
	std::thread threads[WorkerThreadNumber];

	/* Populate the divisor table. This uses a brute force algorithm, and
	 * does not really benefit from a more efficient method  since this is
	 * only performed once. */
	for (long n = 1; n < MaxSeriesLimit; ++n) {
		long buffer[MaxSeriesLimit];
		int length = 0;

		for (long k = 1; k <= n/ 2; ++k) {
			if (n % k == 0) {
				buffer[length++] = k;
			}
		}

		buffer[length++] = n;
		precomputedDivisorList[n] = new long[length];

		for (int k = 0; k < length; ++k) {
			precomputedDivisorList[n][k] = buffer[k];
		}

		precomputedDivisorValues[n] = length;
	}

	/* Header for the LaTeX output. */
	std::cout << "\\documentclass{article}\n"\
		     "\\usepackage[margin=1in]{geometry}\n"\
		     "\\begin{document}\n\n";

	/* Create worker threads to hunt for identities in parallel. */
	for (int n = 0; n < WorkerThreadNumber; ++n) {
		threads[n] = std::thread(worker_thread_entry,
					 std::ref(generator));
	}

	for (int n = 0; n < WorkerThreadNumber; ++n) {
		threads[n].join();
	}

	for (int n = 1; n < MaxSeriesLimit; ++n) {
		delete precomputedDivisorList[n];
	}

	/* Footer for the LaTeX output. */
	std::cout << "\\end{document}\n";

	return 0;
}

