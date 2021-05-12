#include "LiGu_AlgorithmLib/Mat.h"
/******************************************************************************
*                    Schrödinger 方程
*	[定义]: iℏ·∂ψ(r,t)/∂t = [-ℏ/2m▽² + U(r,t)]·ψ(r,t)
		or	iℏ·d/dt·|ψ(t)> = H |ψ(t)>
		[符号]:
			ℏ: 约化Planck常量 = h / 2π    U(r,t): 势
			H: Hamiltonian算子 = -ℏ/2m▽² + U(r)
			ψ: the state vector of the quantum system, letter psi.
		[* Time-independent Schrödinger]:
			[-ℏ/2m▽² + U(r)]·ψ(r)  = Eψ(r)
		or	H|ψ> = E|ψ>
		[符号]: E: 系统的能级, 常量.
		* In the language of linear algebra, this equation is an eigenvalue equation.
		  Therefore, the wave function is an eigenfunction of the Hamiltonian operator
		  with corresponding eigenvalue(s) E.
*	[算法]: 有限差分法
		* ∂ψ/∂x   = [ψ(x+1,...) - ψ(x-1,...)] / 2Δt
		* ∂²ψ/∂x² = [ψ(x+1,...) - 2·ψ(x,...) + ψ(x-1,...)] / Δt²
******************************************************************************/
void Schrodinger_1D(double m, double dx, double xs, double xe, double(*U)(double x)
	, Mat<double>& E, Mat<double>& Psi) {
	const double hbar = 1.05e-34;
	// init
	int N = (xe - xs) / dx;
	Mat<double> U0(N, N);
	for (int i = 0; i < N; i++)U0(i, i) = U(xs + dx * i);
	// 二阶导算子
	Mat<double> D2(N, N);
	for (int i = 0; i < N; i++) {
		D2(i, i) = -2;
		if (i > 0)D2(i, i - 1) = 1;
		if (i < N - 1)D2(i, i + 1) = 1;
	}
	D2.mult(1.0 / (dx * dx), D2);
	D2(0, 0) = D2(0, 1) = D2(1, 0) = 0;								//So that f(xs) = 0
	D2(N - 1, N - 2) = D2(N - 2, N - 1) = D2(N - 1, N - 1) = 0;		//So that f(xe) = 0
	// Hamiltonian算子
	Mat<double> H;
	H.add(H.mult(-hbar * hbar / (2 * m), D2), U0);
	// Time-independent Schrödinger
	H.eig(1e-30, Psi, E);
	E.diag(E);
}

/*
double U_potential(double x) { if (x > 7 * 5e-10 || x < 3 * 5e-10)return 10000; else return 0; }
int main() {
	Mat<double> E, Psi;
	Schrodinger_1D(1e-29, 1e-10, 0, 5e-9, U_potential, E, Psi);
	for (int i = 0; i < E.rows; i++)printf("%e ", E[i]); printf("\n\n");
	for (int i = 0; i < Psi.rows; i++) {
		for (int j = 0; j < Psi.cols; j++) { printf("%e ", Psi(i, j)); }printf("\n");
	}
}*/