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
******************************************************************************/
template<typename FEx, typename FEy, typename FEz, typename FHx, typename FHy, typename FHz>
void Electromagnetics(Mat<>& x, Mat<>& dx, double dt, 
	FEx&& Ex, FEy&& Ey, FEz&& Ez, FHx&& Hx, FHy&& Hy, FHz&& Hz,
	Mat<>& E_ans, Mat<>& H_ans, double mu = 1, double ep = 1, double sigmaE = 0, double sigmaM = 0
) {
	double
		t1 = sigmaE * dt / 2 / ep, CA = (1 - t1) / (1 + t1), CB = dt / ep / (1 + t1),
		t2 = sigmaM * dt / 2 / mu, CP = (1 - t2) / (1 + t2), CQ =-dt / mu / (1 + t2);
	Mat<> tMat, H(x.rows), E(x.rows);
	H_ans.add(H_ans.mult(CP, H.getData(Ex(x), Ey(x), Ez(x))), tMat.mult(CQ, Calculus::Curl(x, dx, Ex, Ey, Ez, tMat)));
	E_ans.add(E_ans.mult(CA, E.getData(Ex(x), Ey(x), Ez(x))), tMat.mult(CB, Calculus::Curl(x, dx, Hx, Hy, Hz, tMat)));
}
/******************************************************************************
*                   Eular 理想流体方程
*	[公式]: ∂\vec u/∂t + (\vec v·▽)\vec v = - 1/ρ·▽p + \vec g
		分量式:
			∂u_x/∂t + v_x·∂v_x/∂x+ v_y·∂v_x/∂y + v_z·∂v_x/∂z = - 1/ρ·∂p/∂x + g_x
*****************************************************************************
template<typename F1, typename F2, typename F3, typename F4>
Mat<>& Eular(Mat<>& x, Mat<>& dx, double density, F1&& vx, F2&& vy, F3&& vz, F4&& p, Mat<>& g, Mat<>& ans) {
	(Calculus::Grad(x, dx, p, ans) *= -1 / density) += g;
	Mat<> GradVx, GradVy, GradVz, vt;
	vt.zero(3).getData(vx(x), vy(y), vz(z));
	ans[0] -= vx.dot(Grad(x, dx, vx, GradVx));
	ans[1] -= vy.dot(Grad(x, dx, vy, GradVy));
	ans[2] -= vz.dot(Grad(x, dx, vz, GradVz));
	return ans;
}*/
/*----------------[ Electrostatic Field 静电场 ]----------------
*	[公式]: ▽·E = ρ/ε0    ▽×E = 0
	[电势]: E = -▽φ
		=>	▽²φ = -ρ/ε0		Poisson's方程
		当ρ=0时, ▽²φ = 0	 Laplace's方程
		[三维直角坐标系]: ∂²φ/∂x² + ∂²φ/∂y² + ∂²φ/∂z² = -1/ε0·ρ(x,y,z)
	[算法]: 即 解Poisson's方程，格林函数，得 φ(r) = - 4π/ε0 ∫∫∫ f(rt) / |r-rt| d³rt
		//Tensor<double>* PoissonEquation(Mat<double>st, Mat<double>ed, Mat<double> delta, double (*f) (Mat<double>& x));
	[静电场唯一性定理]:
		对于各种边界条件，Poisson's方程可能有许多种解，但每个解的梯度相同. 
		静电场情况下, 意味着在边界条件下的满足泊松方程的势函数，所解得的电场确定且唯一.
**----------------------------------------------------------------------*/

/*----------------[ Magnetostatic Field 恒磁场 ]----------------
*	[公式]: ▽·H = 0    ▽×H = 4π/c·J
	[磁矢势]: H = -▽×A		(▽·A = 0)
		=>	▽²A = -4π/c·J		Poisson's方程
	[算法]: 即 解Poisson's方程，格林函数, 得 A = 1/c ∫∫∫ J/r dV
		//Tensor<double>* PoissonEquation(Mat<double>st, Mat<double>ed, Mat<double> delta, double (*f) (Mat<double>& x));
**----------------------------------------------------------------------*/


#endif
