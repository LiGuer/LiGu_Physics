#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/Tensor.h"
#include "Calculus.h"
/******************************************************************************
*                    计算电磁学 -  FDTD
*	[公式]: 麦克斯韦方程组 (旋度那两个)
	[1] ▽ × H = ∂D/∂y + J
		∂Hz/∂y - ∂Hy/∂z = ε∂Ex/∂t + σEx
		∂Hx/∂z - ∂Hz/∂x = ε∂Ey/∂t + σEy
		∂Hy/∂x - ∂Hx/∂y = ε∂Ez/∂t + σEz
	[2] ▽ × E = - ∂B/∂y - Jm
		∂Ez/∂y - ∂Ey/∂z =-μ∂Hx/∂t - σHx
		∂Ex/∂z - ∂Ez/∂x =-μ∂Hy/∂t - σHy
		∂Ey/∂x - ∂Ex/∂y =-μ∂Hz/∂t - σHz
-------------------------------------------------------------------------------
*	[Example]:
	Tensor<> H[3], E[3], Htmp[3], Etmp[3]; Mat<> x(3), dx(3), Et, Ht;  dx.fill(2);
	for (int i = 0; i < 3; i++) {
		H[i].zero(N, N, N); Htmp[i].zero(N, N, N);
		E[i].zero(N, N, N); Etmp[i].zero(N, N, N);
	}
	for (int i = 0; i < E[0].size(); i++) {
		x.getData(E[0].i2x(i) * dx[0], E[0].i2y(i) * dx[1], E[0].i2z(i) * dx[2]);
		Electromagnetics(
			x, dx, 1,
			[&E, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : E[0](x, y, z); },
			[&E, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : E[1](x, y, z); }, 
			[&E, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : E[2](x, y, z); },
			[&H, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : H[0](x, y, z); },
			[&H, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : H[1](x, y, z); },
			[&H, &dx](Mat<>& pos) { int x = pos[0] / dx[0], y = pos[1] / dx[1], z = pos[2] / dx[2]; return (x < 0 || y < 0 || z < 0 || x >= N || y >= N || z >= N) ? 0.0 : H[2](x, y, z); },
		Et, Ht); for (int j = 0; j < Et.size(); j++) Etmp[j][i] = Et[j], Htmp[j][i] = Ht[j];
	}
	for (int j = 0; j < 3; j++) E[j] = Etmp[j], H[j] = Htmp[j]; E[2](N / 2, N / 2, N / 2) = sin(t / (2 * PI));
******************************************************************************/
template<typename FEx, typename FEy, typename FEz, typename FHx, typename FHy, typename FHz>
void Electromagnetics(Mat<>& x, Mat<>& dx, double dt, 
	FEx&& Ex, FEy&& Ey, FEz&& Ez, FHx&& Hx, FHy&& Hy, FHz&& Hz,
	Mat<>& E_ans, Mat<>& H_ans, double mu = 1, double ep = 1, double sigmaE = 0, double sigmaM = 0
) {
	double t1 = sigmaE * dt / 2 / ep,
		   t2 = sigmaM * dt / 2 / mu;
	Mat<> tMat;
	H_ans.add(
		H_ans.mul((1 - t2) / (1 + t2), H_ans.alloc(x.rows).getData(Ex(x), Ey(x), Ez(x))),
		tMat. mul(-dt / mu / (1 + t2), Calculus::Curl(x, dx, Ex, Ey, Ez, tMat))
	);
	E_ans.add(
		E_ans.mul((1 - t1) / (1 + t1), E_ans.alloc(x.rows).getData(Ex(x), Ey(x), Ez(x))), 
		tMat. mul( dt / ep / (1 + t1), Calculus::Curl(x, dx, Hx, Hy, Hz, tMat))
	);
}
/******************************************************************************
						静电场  /  恒磁场
*	[公式]: 
		▽²A  = -4π/c·J		磁矢势 (恒磁场)
		▽²φ = -ρ/ε0			电  势 (静电场)
		电场强度: E = -▽ φ
		磁场强度: H = -▽×A
		▽·E = ρ/ε0			(静电场)
		▽×E = 0
		▽·H = 0				(恒磁场)
		▽×H = 4π/c·J
*	[算法]:	Poisson's方程		
		当ρ=0时, ▽²φ = 0		Laplace's方程
		解Poisson's方程，Green's函数，得 φ(r) = - 4π/ε0 ∫∫∫ f(rt) / |r-rt| d³rt
*	[静电场唯一性定理]:
		对于各种边界条件，Poisson's方程有许多种解，但每个解梯度相同.
		静电场下, 意味边界条件下满足Poisson's方程的势函数，所解得电场唯一确定.
-------------------------------------------------------------------------------
*	[Example]://静电场
	Mat<> Phi(N, N), x(2), dx(2), St(2), Ed(2); dx.fill(1); Ed.fill(N);
		for (int i = 0; i < Phi.size(); i++)
			Phi[i] = Calculus::PoissonEquation( x.getData(Phi.i2x(i), Phi.i2y(i)), dx, St, Ed, [](Mat<>& x) { return sin(x.norm() / (10 * 2 * PI)); });
		Mat<> Ex(N, N), Ey(N, N),Et(2);
		for (int i = 0; i < Phi.size(); i++) {
			Calculus::Grad(x.getData(Phi.i2x(i), Phi.i2y(i)), dx, [&Phi](Mat<>& _x) {
				int x = _x[0], y = _x[1];
				x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
				y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
				return Phi(x, y);
			}, Et); Ex[i] = Et[0]; Ey[i] = Et[1];
		}
******************************************************************************/


/******************************************************************************
*                   Navier Stokes 流体方程
*	[公式]: ∂\vec u/∂t + (\vec v·▽)\vec v =  + η/ρ▽²\vec v - 1/ρ·▽p + \vec g
******************************************************************************/
template<typename F1, typename F2, typename F3, typename F4>
Mat<>& NavierStokesEquations(Mat<>& x, Mat<>& dx, double dt, double density, F1&& vx, F2&& vy, F3&& vz, F4&& p, Mat<>& g, Mat<>& ans) {
	(Calculus::Grad(x, dx, p, ans) *= -1 / density) += g;
	Mat<> vt(x.rows), tmp; vt.getData(vx(x), vy(x), vz(x));
	ans[0] += -vt.dot(Calculus::Grad(x, dx, vx, tmp)) + Calculus::LaplaceOperator(x, dx, vx);
	ans[1] += -vt.dot(Calculus::Grad(x, dx, vy, tmp)) + Calculus::LaplaceOperator(x, dx, vy);
	ans[2] += -vt.dot(Calculus::Grad(x, dx, vz, tmp)) + Calculus::LaplaceOperator(x, dx, vz);
	//if (ans[0] != 0)printf("%f \n", ans[0]);
	ans.add(vt, (ans *= dt));
	return ans;
}

/******************************************************************************
*                   Eular 理想流体方程
*	[公式]: ∂\vec u/∂t + (\vec v·▽)\vec v = - 1/ρ·▽p + \vec g
		分量式:
			∂u_x/∂t + v_x·∂v_x/∂x+ v_y·∂v_x/∂y + v_z·∂v_x/∂z = - 1/ρ·∂p/∂x + g_x
******************************************************************************/
template<typename F1, typename F2, typename F3, typename F4>
Mat<>& Eular(Mat<>& x, Mat<>& dx, double dt, double density, F1&& vx, F2&& vy, F3&& vz, F4&& p, Mat<>& g, Mat<>& ans) {
	(Calculus::Grad(x, dx, p, ans) *= -1 / density) += g;
	Mat<> vt(x.rows), tmp; vt.getData(vx(x), vy(x), vz(x));
	ans[0] += -vt.dot(Calculus::Grad(x, dx, vx, tmp));
	ans[1] += -vt.dot(Calculus::Grad(x, dx, vy, tmp));
	ans[2] += -vt.dot(Calculus::Grad(x, dx, vz, tmp));
	//if (ans[0] != 0)printf("%f \n", ans[0]);
	ans.add(vt, (ans *= dt));
	return ans;
}


#endif
