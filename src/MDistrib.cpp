#include "MDistrib.h"
#include "distrib_header.h"

MDistrib::MDistrib() : mseed(0), weight(1) {
	std::cout << "=== MDistrib::MDistrib()" << std::endl;
	size_t size = sizeof(mseed);
	std::ifstream urandom("/dev/urandom", std::ios::in|std::ios::binary);
	
	// check if stream is open
	if(urandom) {
		urandom.read(reinterpret_cast<char*>(&mseed), size); // read from urandom
		// check if stream is ok, read succeeded
		if(!urandom) {
			std::cerr << "*** Could not read from /dev/urandom - abort !" << std::endl;
			exit(-1);
		}
		urandom.close(); //close stream
	}
	else {
		std::cerr << "*** Could not open /dev/urandom - abort !" << std::endl;
		exit(-1);
	}
	
	// set seed
	random_engine.seed(mseed);
	real_distribution = std::uniform_real_distribution<double>(0.0, 1.0);
}
MDistrib::~MDistrib() {
	std::cout << "=== MDistrib::~MDistrib()" << std::endl;
}

void MDistrib::SetSeed(unsigned int seed) {
	mseed = seed;
	random_engine.seed(mseed);
}
double MDistrib::GetUniform(const double& a, const double& b) {
	if(a == b) {
		weight = 1;
		return a;
	}
	
	const double u = real_distribution(random_engine);
	double x = a + (b - a) * u;
	weight = b - a;
	return x;
}

double MDistrib::GetNormal(const double& mu, const double& sigma) {
	if(sigma == 0) {
		weight = 1;
		return mu;
	}
	
	const double u = real_distribution(random_engine);
	const double v = real_distribution(random_engine);
	
	const double two_pi = 2 * M_PI;
	const double mag = std::sqrt(-2 * std::log(u));
	
	double x = sigma * mag * std::cos(two_pi * v) + mu;
	
	// compute weight = inv(pdf)
	const double xi = (x - mu) / sigma;
	const double nume = std::exp(-0.5 * xi * xi);
	const double deno = sigma * std::sqrt(two_pi);
	weight = deno / nume;
	return x;
}

double *MDistrib::GetSolidAngle(double angles[2], const double& azMin, const double& azMax, 
		const double& elMin, const double& elMax) {
	const double deg = M_PI / 180;
	const double u = real_distribution(random_engine);
	const double v = real_distribution(random_engine);
	
	double azWeight, elWeight;
	if(azMin == azMax) {
		azWeight = 1;
		angles[0] = azMin;
	}
	else {
		azWeight = deg * (azMax - azMin); 	// in rad
		angles[0] = azMin + u * (azMax - azMin);
	}
	
	const double sinMin = std::sin(elMin * deg);
	const double sinMax = std::sin(elMax * deg);
	if(sinMin == sinMax) {
		elWeight = 1;
		angles[1] = elMin;		//FIXME sinus equality does not mean angles equality unless same quadrant
	}
	else {
		elWeight = sinMax - sinMin;			// in rad
		angles[1] = std::asin(sinMin + v * (sinMax - sinMin));
		angles[1] /= deg;					// in deg
	}
	
	weight = azWeight * elWeight;
	
	return angles;
}

double MDistrib::GetLogarithm(const double& min, const double& max) {
	if(min == max) {
		weight = 1;
		return min;
	}
	
	double u = real_distribution(random_engine);
	double r = std::log(max / min);
	double x = min * std::exp(r * u);
	weight = r * x;							// Let's assume that min > 0
	return x;
}

// algorithm https://fr.wikipedia.org/wiki/Loi_de_Poisson
double MDistrib::GetPoissonKnuth(const double& lambda) {
	double p = 1;
	long long k = 0;
	
	while(p > std::exp(-lambda)) {
		const double u = real_distribution(random_engine);
		if ((u <= 0.) || (u >= 1.)) continue;
		
		p *= u;
		k += 1;
	}
	
	double x = (double) (k - 1);
	
	long long factK = 1;
	for(long long i = 2; i < k - 1; i++) {
		factK *= i;
	}
	
	const double pdf = std::exp(-lambda) * std::pow(lambda, x) / (double)factK;
	weight = 1 / pdf;
	
	return x;
}
// reference: https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
double MDistrib::GetPoissonCook(const double& lambda) {
	const double c = 0.767 - 3.36 / lambda;
	const double beta = M_PI / std::sqrt(3.0 * lambda);
	const double alpha = beta * lambda;
	const double k = std::log(c) - lambda - std::log(beta);
	
	while(true) {
		const double u = real_distribution(random_engine);
		const double x = (alpha - std::log((1.0 - u)/u)) / beta;
		const long n = std::floor(x + 0.5);
		if(n < 0) continue;
		
		const double v = real_distribution(random_engine);
		const double y = alpha - beta * x;
		const double t = 1.0 + std::exp(y);
		const double lhs = y + std::log(v / (t*t));
		// lgamma c++ function https://en.cppreference.com/w/cpp/numeric/math/lgamma
		// N.B log(n!) ~= lgamma(n + 1)
		const double rhs = k + n*std::log(lambda) - std::lgamma(n + 1);
		
		if(lhs <= rhs) return n;
	}
}
// https://www.johndcook.com/blog/2010/06/14/generating-poisson-random-values/
double MDistrib::GetPoisson(const double& lambda) {
	if(lambda < 30) return this->GetPoissonKnuth(lambda);
	else return this->GetPoissonCook(lambda);
}

double MDistrib::GetUniformBinCenter(const double& min, const double& max, const double& delta) {
	if(min == max) {
		weight = 1;
		return min;
	}
	
	double value = MAX_DOUBLE;
	
	const long int n = std::floor((max - min) / delta);
	const double pdf = 1.0 / n;
	weight = 1 / pdf;
	// get random value
	double u = real_distribution(random_engine);
	
	double cdf = 0;
	for(long int i = 1; i < n + 1; i++) {
		cdf += pdf;
		
		if(u < cdf) {
			value = min + (i - 0.5) * delta;
			break;
		}
	}
	return value;
}
