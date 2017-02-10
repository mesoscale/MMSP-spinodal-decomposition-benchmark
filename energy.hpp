// energy.hpp
// Energy functions for 2D and 3D Cahn-Hilliard model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef CAHNHILLIARD_ENERGY
#define CAHNHILLIARD_ENERGY
#include<cmath>

const double deltaX = 1.0;
const double Ca = 0.3;
const double Cb = 0.7;
const double rhoS = 5.0;
const double kappa = 2.0;
const double M = 5.0;
const double C0 = 0.5;
const double epsilon = 0.01;
const double CFL = 0.25;
const double dt = std::pow(deltaX, 4)*CFL/(24.0*M*kappa);

double chemenergy(const double& C)
{
	// Equation 2
	const double A = C-Ca;
	const double B = Cb-C;
	return rhoS * A*A * B*B;
}

double dfdc(const double& C)
{
	// d(chemenergy)/dc
	const double A = C-Ca;
	const double B = Cb-C;
	return 2.0 * rhoS * A * B * (Ca + Cb - 2.0 * C);
}

double cheminit(const double& x, const double& y)
{
	// Equation 12
	return C0 + epsilon * ( std::cos(0.105*x)          * std::cos(0.11*y)
	                      + std::pow(std::cos(0.13*x)  * std::cos(0.087*y), 2.0)
	                      + std::cos(0.025*x - 0.15*y) * std::cos(0.07*x - 0.02*y)
                          );
}

#endif
