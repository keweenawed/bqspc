
namespace bqspc
{

/* Largest coefficient to truncate all series computations at. */
const static int MaxSeriesLimit = 100;

/* Largest modulus q-series product to consider in the search. */
const static int MaxProductSignatureLength = 50;

/* Largest number of q-Pochhammer symbols to allow on the sum side. */
const static int MaxNumberQPS = 2;

/* Largest degree 1 and 2 coefficients for the sum side power of q.*/
const static int MaxQPowerDegree1 = 5;
const static int MaxQPowerDegree2 = 5;

/* Largest degree 1 coefficient for the sum side power of z.*/
const static int MaxZPower = 2;

/* Largest values the parameters in the q-Pochhammer symbols on the sum side
 * can take. Written $(z^a q^b; q^c)^d_{en+f}^f)$, these are in alphabetical
 * order. */
const static int MaxQPSParameters[6] = {0, 6, 6, 4, 3, 3};

/* Number of worker threads to use. */
const static int WorkerThreadNumber = 8;

/* Size of the job queue each worker thread uses. */
const static int MaxJobQueueLimit = 20;

class Parameters;

/* Parent class for both kinds of series. */
class Series
{
protected:

	/* The power of the coefficient to truncated computations at. This can
	 * be at most the value of MaxSeriesLimit. */
	int qLimit;
};

/* Stores the parameters of infinite q-Pochhammer symbols that determine the
 * product side of a Rogers-Ramanujan type identity. */
class ProductSignature
{
public:
	int period;
	long powers[MaxProductSignatureLength];

	long dilation(void);
};

/* Encodes a truncated univariate series. */
class SeriesUv : public Series
{
private:
	/* The coefficient of q^n is given by coefficients[n]. */
	long coefficients[MaxSeriesLimit];

	void qPochhammer(int, int, int, int, int);

public:
	SeriesUv(void) {this->qLimit = MaxSeriesLimit;}

	void qSeries(class Parameters&);
	void factorize(ProductSignature&);

	inline bool operator==(const class SeriesUv& series)
	{
		for (int n = 0; n < this->qLimit; ++n) {
			if (this->coefficients[n] != series.coefficients[n]) {
				return false;
			}
		}

		return true;
	}

	inline class SeriesUv operator+(const class SeriesUv& series)
	{
		class SeriesUv result;

		result.qLimit = this->qLimit;

		for (int n = 0; n < this->qLimit; ++n) {
			result.coefficients[n] = this->coefficients[n]
					       + series.coefficients[n];
		}

		return result;
	}

	inline class SeriesUv& operator+=(const class SeriesUv& series)
	{
		class SeriesUv copy;

		copy = *this;
		*this = copy + series;
		return *this;
	}

	inline class SeriesUv operator*(const class SeriesUv& series)
	{
		class SeriesUv result;

		result.qLimit = this->qLimit;

		/* Compute the Cauchy product. */
		for (int n = 0; n <= this->qLimit; ++n) {
			result.coefficients[n] = 0;

			for (int k = 0; k <= n; ++k) {
				result.coefficients[n]
					+= this->coefficients[k]
					* series.coefficients[n - k];
			}
		}

		return result;
	}

	inline class SeriesUv& operator*=(const class SeriesUv& series)
	{
		class SeriesUv copy;

		copy = *this;
		*this = copy * series;
		return *this;

	}

	inline class SeriesUv operator-(void)
	{
		class SeriesUv result;

		result.qLimit = this->qLimit;

		for (int n = 0; n < this->qLimit; ++n) {
			result.coefficients[n] = -this->coefficients[n];
		}

		return result;
	}

	inline void zero(void)
	{
		for (int n = 0; n < this->qLimit; ++n) {
			this->coefficients[n] = 0;
		}
	}

	/* Computes the truncated q-binomial Coefficient, which is defined
	 * $(q;q)_{top}/((q;q)_{bottom}(q;q){top - bottom})$. */
	inline void qBinomial(int top, int bottom)
	{
		class SeriesUv factor;

		this->zero();
		this->coefficients[0] = 1;
		factor.qLimit = this->qLimit;
		factor.qPochhammer(1, 1, 1, 1, top);
		*this *= factor;
		factor.qPochhammer(1, 1, 1, -1, bottom);
		*this *= factor;
		factor.qPochhammer(1, 1, 1, -1, top - bottom);
		*this *= factor;
	}

	/* Multiplies the series by $q^{power}$. */
	inline void translate(int qPower)
	{
		class SeriesUv copy = *this;

		for (int n = 0; n < qPower; ++n) {
			this->coefficients[n] = 0;
		}

		for (int n = 0; n < this->qLimit - qPower; ++n) {
			this->coefficients[n + qPower] = copy.coefficients[n];
		}
	}
};

/* Encodes a truncated bivariate series. This is currently not used in the
 * program, but will be in future versions. */
class SeriesBv : public Series
{
private:

	/* All coefficients at or past $q^{qLimit}z^{zLimit}$ are truncated.
	 * This value can be at most MaxSeriesLimit. */
	int zLimit;

	/* The coefficient of $q^n z^k$ is coefficients[n][k]$. */
	long coefficients[MaxSeriesLimit][MaxSeriesLimit];

	void qPochhammer(int, int, int , int, int, int);

public:
	SeriesBv(void)
	{
		this->qLimit = MaxSeriesLimit;
		this->zLimit = MaxSeriesLimit;
	}

	void qSeries(class Parameters&);
	void generalizeProduct(ProductSignature&);

	inline bool operator==(const class SeriesBv& series)
	{
		for (int n = 0; n < this->qLimit; ++n) {
			for (int k = 0; k < this->zLimit; ++k) {
				if (this->coefficients[n][k]
				    != series.coefficients[n][k]) {
					return false;
				}
			}
		}

		return true;
	}

	inline bool operator!=(const class SeriesBv& series)
	{
		return !(*this == series);
	}

	inline class SeriesBv operator+(const class SeriesBv& series)
	{
		class SeriesBv result;

		result.qLimit = this->qLimit;
		result.zLimit = this->zLimit;

		for (int n = 0; n < this->qLimit; ++n) {
			for (int k = 0; k < this->zLimit; ++k) {
				result.coefficients[n][k]
					= this->coefficients[n][k]
					+ series.coefficients[n][k];
			}
		}

		return result;
	}

	inline class SeriesBv& operator+=(const class SeriesBv& series)
	{
		class SeriesBv copy;

		copy = *this;
		*this = copy + series;
		return *this;
	}

	inline class SeriesBv operator*(const class SeriesBv& series)
	{
		class SeriesBv result;

		result.qLimit = this->qLimit;
		result.zLimit = this->zLimit;
		result.zero();

		/* Double Cauchy product. */
		for (int n1 = 0; n1 < this->qLimit; ++n1) {
			for (int n2 = 0; n2 <= n1; ++n2) {
				long buffer[MaxSeriesLimit];
				const long *c1 = this->coefficients[n2];
				const long *c2 = series.coefficients[n1 - n2];

				for (int k1 = 0; k1 < this->zLimit; ++k1) {
					buffer[k1] = 0;

					for (int k2 = 0; k2 <= k1; ++k2) {
						buffer[k1] += c1[k1 - k2]
							    * c2[k2];
					}
				}

				for (int k = 0; k < this->zLimit; ++k) {
					result.coefficients[n1][k]
						+= buffer[k];
				}
			}
		}

		return result;
	}

	inline class SeriesBv& operator*=(const class SeriesBv& series)
	{
		class SeriesBv copy;

		copy = *this;
		*this = copy * series;
		return *this;

	}

	inline class SeriesBv operator-(void)
	{
		class SeriesBv result;

		result.qLimit = this->qLimit;
		result.zLimit = this->zLimit;

		for (int n = 0; n < this->qLimit; ++n) {
			for (int k = 0; k < this->zLimit; ++k) {
				result.coefficients[n][k]
					= -this->coefficients[n][k];
			}
		}

		return result;
	}

	inline void zero(void)
	{
		for (int n = 0; n < this->qLimit; ++n) {
			for (int k = 0; k < this->zLimit; ++k) {
				this->coefficients[n][k] = 0;
			}
		}
	}

	/* Multiplies the series by $q^{qPower}z^{zPower}$. $*/
	inline void translate(int qPower, int zPower)
	{
		class SeriesBv copy = *this;

		this->zero();

		for (int n = 0; n < this->qLimit - qPower; ++n) {
			for (int k = 0; k < this->zLimit - zPower; ++k) {
				this->coefficients[n + qPower][k + zPower]
					= copy.coefficients[n][k];
			}
		}
	}
};

/* The parameters that determine the sum side of a q-series identity. */
class Parameters
{
public:
	/* Coefficients for the degree 1 and 2 terms on the power of q. */
	int qScalarDeg1;
	int qScalarDeg2;

	/* Set true if the power of q is to be divided by 2. */
	bool dividePowerBy2;

	/* Set true to put an alternating sign on the series. */
	bool alternatingSign;

	/* Number of q-Pochhammer symbols. */
	int qPSLength;

	/* Coefficient for the degree 1 term on the power of z. */
	int zScalar;

	/* The q-Pochhammer symbols are stored so that the nth from the left
	 * begins at the address of qPS[n * 6]. */
	int qPS[MaxNumberQPS * 6];
};

class WorkerThread;

/* Universal interface that threads use to get more work. */
class ParameterGenerator: private Parameters
{
private:

	/* Set true until all parameter combinations are exhausted. */
	bool keepWorking;

	void advanceState(void);
	void checkParameterRedundancy(void);

public:
	void populateJobQueue(WorkerThread&);

	ParameterGenerator(void)
	{
		this->qScalarDeg1 = 0;
		this->qScalarDeg2 = 1;
		this->dividePowerBy2 = false;
		this->alternatingSign = false;
		this->qPSLength = 0;
		this->zScalar = 0;

		for (int n = 0; n < MaxNumberQPS * 6; ++n) {
			this->qPS[n] = 0;
		}

		this->keepWorking = true;
	}
};

/* Data and methods for each worker thread. */
class WorkerThread
{
private:

	/* Buffer of parameters to try that belong only to the worker. */
	int jobQueueLength;
	Parameters jobQueue[MaxJobQueueLimit];

	void tryCombinationUv(Parameters&);
	void reportIdentity(Parameters&, ProductSignature&);

public:
	friend class ParameterGenerator;

	void jobLoop(ParameterGenerator&);
};

};

