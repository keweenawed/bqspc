
namespace bqspc {

/* The largest allowed coefficient to truncate $q$-series computations at. */
const static int MaxSeriesLimit = 100;

/* Longest pattern of powers in the truncated product to search for. */
const static int MaxProductSignatureLength = 50;

/* The largest number of $q$-series summation indices allowed. */
const static int MaxIndices = 2;

/* The largest number of distinct $q$-Pochhammer symbols allowed in a
 * single $q$-series, ignoring multiplicity. */
const static int MaxQPS = 2;

/* Number of $q$-series parameters to cache per worker thread. */
const static int JobQueueLimit = 100;

/* The parameters that fully determine a particular $q$-series of the form
 * $\sum_{n_0, \dots, n_\ell \geq 0} (-1)^{d \times (n_0 + \cdots + n_\ell))}
 * \times q^{c(n_0 \dots, n_\ell)} \prod_{i=0}^k
 * (\pm q^{a_i};q^{b_i})_{s_i(n_0, \dots, n_\ell)}^{p_i}$, where $\ell \geq 0$
 * determines the number of summation indices to use, $d \in \{0, 1\}$, $c$ is
 * a nonnegative degree 1 or 2 polynomial in $n_0, dots, n_\ell$, $k \geq 0$
 * determines the number of $q$-Pochhammer symbols, $p_i$ is any nonzero power
 * on them, $s_i$ is a linear function in $n_0, \dots, n_\ell$ that is not
 * identically zero, and the dilations satisfy $a_i, b_i \geq 1$ with
 * $a_i = 0$ additionally allowed if $p_i > 0$. */
class Parameters
{
public:

	/* Set to true to for $d=1$, and false for $d=0$. */
	bool alternatingSign;

	/* Set to true to divide the coefficients of $c$ all by 2. This is only
	 * permitted if doing so preserves $c$ being integer valued. */
	bool dividePowerBy2;

	/* The number of summation indices $\ell + 1$ to use. */
	int indicesInUse;

	/* The number of distinct q-Pochhammer symbols $k + 1$ to use. */
	int qPSInUse;
	
	/* Coefficients of the pure degree 1 and 2 powers in $c$, respectively,
	 * stored so the coefficient of $n_i$ and $n_i$^2 is at index $i$. */
	int qScalarsDegree1[MaxIndices];
	int qScalarsDegree2Pure[MaxIndices];

	/* Here, for $0 \leq i < j \leq \ell$, the coefficient of $n_in_j$ is
	 * stored respecting the ordering $n_0n_1, \dots, n_0n_\ell, n_1n_2,
	 * \dots, n_1n_\ell, \dots, n_{\ell-1}n_\ell$. */
	int qScalarsDegree2Mixed[MaxIndices * (MaxIndices - 1) / 2];

	/* Parameters of $(\pm q^{a_i};q^{b_i})_{s(n_0, \dots, n_\ell)}^{p_i}$ are
	 * stored at index $i$, where $\pm$ is given by negativePrefix, dilation1
	 * is $a_i$, dilation2 is $b_i$, power is $p_i$, and the coefficients of
	 * $s_i(n_0, \dots, n_\ell) = f_0 n_0 + \cdots + f_\ell n_\ell$ are stored
	 * in subScalars so $f_j$ is at index $j$. */
	class {
	public:
		int dilation1;
		int dilation2;
		bool negativePrefix;
		int power;
		int subScalars[MaxIndices];
	} qPS[MaxQPS];
};

class WorkerThread;

/* Iterates through all parameter combinations to provide the threads work. */
class ParameterGenerator : private Parameters
{
	friend class WorkerThread;

	bool continueWorking;

	void advance(void);
	void populate(WorkerThread&);

public:
	ParameterGenerator(void);
};

class QSeries;

/* Encodes the product signature of a truncated $q$-series, which is the
 * minimal length pattern $a_1, \dots, a_\ell$ of repeating powers appearing
 * when the $q$-series is factored as the infinite product of geometric series 
 * $\prod_{n \geq 1} \frac{1}{(1-q^n)^{a_n}}$, if such a pattern exists. */
class ProductSignature
{
	friend class QSeries;
	friend class WorkerThread;

	/* The length of the pattern. */
	int period;

	/* The sequence of powers that forms the pattern. */
	long powers[MaxProductSignatureLength];

	long pairwiseGCD(long, long);

public:
	long dilation(void);
	void factorize(QSeries&);
};

/* Stores the truncated coefficients of a $q$-series, and provides all
 * functionality for truncated $q$-series arithmetic and manipulations. */
class QSeries
{
	friend class ProductSignature;

	/* The $q$-series coefficients, stored so the coefficient of $q^i$ is at
	 * the index $i$. */
	long coefficients[MaxSeriesLimit];

	/* The coefficient to truncate all computations at. */
	int limit;

	void reciprocal(void);
	void raiseToPower(int);
	void qPochhammer(int, int, bool, int);
	void qBinomial(int, int);
	int qSeriesPower(Parameters&, int (&)[MaxIndices]);
	void qSeriesTerm(Parameters&, int (&)[MaxIndices]);

public:

	void qSeries(Parameters&);

	QSeries(int limit = MaxSeriesLimit) {this->limit = limit;}

	/* Sets all coefficients to zero. */
	inline void zero(void)
	{
		for (int index = 0; index < MaxSeriesLimit; ++index) {
			this->coefficients[index] = 0;
		}
	}

	/* Shifts all coefficients over by the value of power, or equivalently,
	 * computes the truncated multiplication by $q^{power}$. */
	inline void translate(int power)
	{
		QSeries copy = *this;

		for (int index = 0; index < power; ++index) {
			this->coefficients[index] = 0;
		}

		for (int index = 0; index < this->limit - power; ++index) {
			this->coefficients[index + power] = copy.coefficients[index];
		}
	}


	/* All of the following methods that act on two $q$-series assume both
	 * of them have an equal value of limit. The calling functions are always
	 * responsible for ensuring this is true. */

	/* Two $q$-series are deemed equal when they have equal coefficients. */
	inline bool operator==(const QSeries& series)
	{
		for (int index = 0; index < this->limit; ++index) {
			if (this->coefficients[index] != series.coefficients[index]) {
				return false;
			}
		}

		return true;
	}

	/* Adds the coefficients of two $q$-series. */
	inline QSeries operator+(const QSeries& series)
	{
		QSeries result(this->limit);

		for (int index = 0; index < this->limit; ++index) {
			result.coefficients[index] = this->coefficients[index]
									   + series.coefficients[index];
		}

		return result;
	}

	inline QSeries& operator+=(const QSeries& series)
	{
		QSeries copy;

		copy = *this;
		*this = copy + series;
		return *this;
	}

	/* Computes the product of two $q$-series using the quadratic time Cauchy
	 * product. Experimentally, this is faster than the asymptotically
	 * superior Karatsuba's algorithm for this use case. */
	inline QSeries operator*(const QSeries& series)
	{
		QSeries result(this->limit);

		for (int nIndex = 0; nIndex < this->limit; ++nIndex) {
			result.coefficients[nIndex] = 0;

			for (int kIndex = 0; kIndex <= nIndex; ++kIndex) {
				result.coefficients[nIndex]
					+= this->coefficients[kIndex]
					* series.coefficients[nIndex - kIndex];
			}
		}

		return result;
	}

	inline QSeries& operator*=(const QSeries& series)
	{
		QSeries copy = *this;

		*this = copy * series;
		return *this;

	}

	/* Negates the coefficients of the $q$-series. */
	inline QSeries operator-(void)
	{
		QSeries result(this->limit);

		for (int index = 0; index < this->limit; ++index) {
			result.coefficients[index] = -this->coefficients[index];
		}

		return result;
	}
};

/* Data and methods for each worker thread. */
class WorkerThread
{
	friend class ParameterGenerator;

	/* Points to the universal generator. */
	ParameterGenerator *generator;

	/* Cache of parameters to try. */
	Parameters *jobQueue[JobQueueLimit];

	/* Number of parameters in the cache. */
	int jobQueueLength;

	void reportIdentity(Parameters&, ProductSignature&);
	void tryCombination(Parameters&);

public:
	void jobLoop(void);

	WorkerThread(ParameterGenerator *generator)
	{
		this->generator = generator;

		/* The amount of memory this class can take up if new is not used here
		 * may cause a stack overflow. */
		for (int index = 0; index < JobQueueLimit; ++index) {
			jobQueue[index] = new Parameters();
		}
	}

	~WorkerThread(void)
	{
		for (int index = 0; index < JobQueueLimit; ++index) {
			delete jobQueue[index];
		}
	}
};

};

