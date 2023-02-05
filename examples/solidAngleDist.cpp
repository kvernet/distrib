#include "MDistrib.h"

#include <iostream>

int main(int argc, char **argv) {
	if(argc < 6) {
		std::cerr << "Usage:" << std::endl;
		std::cerr << "\t" << argv[0] << " azMin azMax elMin elMax n_events" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	const double azMin = std::stod(argv[1]);
	const double azMax = std::stod(argv[2]);
	const double elMin = std::stod(argv[3]);
	const double elMax = std::stod(argv[4]);
	const long n_events = std::stol(argv[5]);
	
	MDistrib *mDistrib = new MDistrib();
	const long seed = mDistrib->GetSeed();
	printf("=== seed: %ld\nx    weight\n==========================\n", seed);
	
	double angles[2];
	for(long i = 0; i < n_events; i++) {
		mDistrib->GetSolidAngle(angles, azMin, azMax, elMin, elMax);
		double weight = mDistrib->GetWeight();
		printf("(%.3f, %.3f)\t%.3f\n", angles[0], angles[1], weight);
	}
	
	delete mDistrib;
	
	exit(EXIT_SUCCESS);
}
