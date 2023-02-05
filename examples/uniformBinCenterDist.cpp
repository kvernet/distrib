#include "MDistrib.h"

#include <iostream>

int main(int argc, char **argv) {
	if(argc < 5) {
		std::cerr << "Usage:" << std::endl;
		std::cerr << "\t" << argv[0] << " min mMax delta n_events" << std::endl;
		exit(EXIT_FAILURE);
	}
	
	const double min = std::stod(argv[1]);
	const double max = std::stod(argv[2]);
	const double delta = std::stod(argv[3]);
	const long n_events = std::stol(argv[4]);
	
	MDistrib *mDistrib = new MDistrib();
	const long seed = mDistrib->GetSeed();
	printf("=== seed: %ld\nx    weight\n==========================\n", seed);
	
	for(long i = 0; i < n_events; i++) {
		double x = mDistrib->GetUniformBinCenter(min, max, delta);
		double weight = mDistrib->GetWeight();
		printf("%.3f\t%.3f\n", x, weight);
	}
	
	delete mDistrib;
	
	exit(EXIT_SUCCESS);
}
