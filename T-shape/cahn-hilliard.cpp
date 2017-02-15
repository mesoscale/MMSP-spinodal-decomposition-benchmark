// cahn-hilliard.hpp
// Algorithms for 2D and 3D Cahn-Hilliard model
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#ifndef CAHNHILLIARD_UPDATE
#define CAHNHILLIARD_UPDATE
#include"MMSP.hpp"
#include<cmath>
#include"cahn-hilliard.hpp"
#include"../energy.hpp"

namespace MMSP {

bool isOutside(const vector<int>& x)
{
	if ((x[1]<100) && ((x[0]<0) || (x[0]>19)))
		return true;
	return false;
}

// custom vector calc expressions for boundary points
template <int dim, typename T>
T zfLaplacian(const grid<dim,T>& GRID, const vector<int>& x)
{
	T laplacian = 0.0;
	vector<int> s = x;
	const T& y = GRID(x);

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = isOutside(s) ? y : GRID(s);
		s[i] -= 2;
		const T& yl = isOutside(s) ? y : GRID(s);
		s[i] += 1;

		double weight = 1.0 / (dx(GRID, i) * dx(GRID, i));
		laplacian += weight * (yh - 2.0 * y + yl);
	}
	return laplacian;
}

template <int dim, typename T>
vector<T> zfGradient(const grid<dim,T>& GRID, const vector<int>& x)
{
	vector<T> gradient(dim, 0.0);
	vector<int> s = x;
	const T& y = GRID(x);

	for (int i=0; i<dim; i++) {
		s[i] += 1;
		const T& yh = isOutside(s) ? y : GRID(s);
		s[i] -= 2;
		const T& yl = isOutside(s) ? y : GRID(s);
		s[i] += 1;

		double weight = 1.0 / (2.0 * dx(GRID, i));
		gradient[i] = weight * (yh - yl);
	}
	return gradient;
}

template <int dim,typename T>
double Helmholtz(const grid<dim,T>& GRID)
{
	double dV = 1.0;
	for (int d=0; d<dim; d++)
		dV *= dx(GRID, d);

	double f = 0.0;
	double g = 0.0;

	for (int n=0; n<nodes(GRID); n++) {
		vector<int> x = position(GRID, n);
		if (!isOutside(x)) {
			vector<double> gradc = zfGradient(GRID, x);
			f += chemenergy(GRID(x));
			g += gradc*gradc;
		}
	}

	double F = dV*(f + 0.5*kappa*g);
	#ifdef MPI_VERSION
	double myF(F);
        MPI::COMM_WORLD.Barrier();
	MPI::COMM_WORLD.Allreduce(&myF, &F, 1, MPI_DOUBLE, MPI_SUM);
	#endif
	return F;
}

void generate(int dim, const char* filename)
{
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	if (dim!=2 && rank==0) {
		std::cerr<<"ERROR: CHiMaD problems are 2-D, only!"<<std::endl;
		std::exit(-1);
	}

	if (dim==2) {
		GRID2D initGrid(1,-40,60,0,120);

		for (int d=0; d<dim; d++){
			dx(initGrid,d) = deltaX;
			if (x0(initGrid,d)==g0(initGrid,d))
				b0(initGrid,d) = Neumann; // enumerated in MMSP.utility.hpp
			else if (x1(initGrid,d)==g1(initGrid,d))
				b1(initGrid,d) = Neumann; // enumerated in MMSP.utility.hpp
		}

		for (int n=0; n<nodes(initGrid); n++) {
			vector<int> x = position(initGrid, n);
			if (isOutside(x))
				initGrid(x) = C0;
			else
				initGrid(x) = cheminit(dx(initGrid,0)*x[0], dx(initGrid,1)*x[1]);
		}

		ghostswap(initGrid);
		output(initGrid,filename);

	}
}

template <int dim, typename T>
void update(grid<dim,T>& oldGrid, int steps)
{
	ghostswap(oldGrid);

	grid<dim,T> newGrid(oldGrid);
	grid<dim,T> lapGrid(oldGrid);

	for (int step=0; step<steps; step++) {
		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid,n);
			if (isOutside(x)) {
				lapGrid(x) = 0.0;
			} else {
				lapGrid(x) = dfdc(oldGrid(x)) - kappa*zfLaplacian(oldGrid,x);
			}
		}

		ghostswap(lapGrid);

		for (int n=0; n<nodes(oldGrid); n++) {
			vector<int> x = position(oldGrid,n);
			if (isOutside(x)) {
				newGrid(x) = C0;
			} else {
				newGrid(x) = oldGrid(x) + dt*M*zfLaplacian(lapGrid,x);
			}
		}

		swap(oldGrid,newGrid);
		ghostswap(oldGrid);
	}
}

} // MMSP
#endif

#include"../main.cpp"
