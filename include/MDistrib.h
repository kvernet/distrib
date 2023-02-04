#ifndef M_DISTRIB_H
#define M_DISTRIB_H

#include <random>

class MDistrib {
	public:
		MDistrib();
		~MDistrib();
		
		unsigned int GetSeed() const { return mseed; }
		void SetSeed(unsigned int seed);
		double GetWeight() const { return weight; }		
		
		double GetU();
		double GetUniform(const double& a = 0, const double& b = 1);
		double GetNormal(const double& mu = 0, const double& sigma = 1);
		double *GetSolidAngle(double angles[2], const double& azMin = 0, const double& azMax = 1, 
		const double& elMin = 0, const double& elMax = 1);
		double GetLogarithm(const double& min = 0, const double& max = 1);
		double GetPoissonKnuth(const double& lambda = 0);
		double GetPoissonCook(const double& lambda = 0);
		double GetPoisson(const double& lambda = 0);
		double GetUniformBinCenter(const double& min = 0, const double& max = 1, const double& delta = 1);
	
	private:
		unsigned int mseed;
		std::default_random_engine random_engine;
		std::uniform_real_distribution<double> real_distribution;
		double weight;
};

#endif	//M_DISTRIB_H
