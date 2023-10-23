#include <iostream>
#include "vmd.h"
#include <Eigen/core>
#include <math.h>
#define _USE_MATH_DEFINES

using namespace std;

int main() {
	int T = 1000;
	float fs = 1. / T;
	Eigen::VectorXf t = Eigen::VectorXf::LinSpaced(T, 0, 1);
	Eigen::VectorXf input(T);

	float f1 = 2.0f, f2 = 24.0f, f3 = 288.0f;
	for (int i = 0; i < T; i++) {
		input(i) = cos(2 * M_PI * f1 * t(i)) 
				 + cos(2 * M_PI * f2 * t(i)) 
				 + cos(2 * M_PI * f3 * t(i));
	}
	input += 0.1 * Eigen::VectorXf::Random(T);

	vmd result(input, 3, 1e-7, 2000);

	cout << result.omega;
	//cout << result.u;

	
	return 0;
}