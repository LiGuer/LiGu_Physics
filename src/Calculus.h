#ifndef CALCULUS_H
#define CALCULUS_H
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/Tensor.h"
#include "../../LiGu_AlgorithmLib/NumberTheory.h"
#include <functional>
#include <complex>
/******************************************************************************
*                    微积分 / 微分方程
-------------------------------------------------------------------------------
---------------- 微分 ----------------
double diff		(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);		\\导数 (N阶)
double diff_	(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);
double Curvature(double x0, double(*y)(double x),			 double dx = 1E-3);		\\曲率
double PartiDeriv	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数
double PartiDeriv2	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数 (2阶)
double LaplaceOperator	(Mat<>& x, Mat<>& dx, F&& f);								\\Laplace算子
Mat<>& Grad		(Mat<>& x, Mat<>& dx, F&& f, Mat<>& ans);							\\梯度
double Div		(Mat<>& x, Mat<>& dx, FX&& fx, FY&& fy, FZ&& fz);					\\散度
Mat<>& Curl		(Mat<>& x, Mat<>& dx, FX&& f0, FY&& f1, FZ&& f2, Mat<>& ans);		\\旋度
---------------- 级数展开 ------------
Mat<>& TaylorFormula(double x0, double(*f)(double x), Mat<>& Coeff, int N = 3);		\\Taylor展开
double Exp		(double x, int N = 18);												\\常用函数
double Sin		(double x, int N = 18);
double Cos		(double x, int N = 18);
double lnOneAdd	(double x, int N = 18);
double Arctan	(double x, int N = 18);
double PowOneAdd(double x, double p, int N = 18);
---------------- 积分 ----------------
double integral(double xSt, double xEd, F&& f, int n);								\\积分
---------------- 微分方程 ------------
void   RungeKutta			(Mat<>& y, double dx, double x0, F&& f, int enpoch = 1);		\\解常微分方程组: RungeKutta法
double PoissonEquation		(Mat<>& x, Mat<>& dx, F&& u);									\\Poisson's方程
double WaveEquation			(Mat<>& x, Mat<>& dx, double t, double dt, double A, F&& u);	\\波动方程
double DiffusionEquation	(Mat<>& x, Mat<>& dx, double t, double dt, double A, F&& u);	\\扩散方程
******************************************************************************/
namespace Calculus {
#define PI 3.141592653589
/*#############################################################################

								微分

#############################################################################*/
/******************************************************************************
*                    导数  N阶
*	[定义]:
		导数: df/dx = lim_(x->0)  [ f(x+Δx) - f(x-Δx) ] / (2 Δx)
		N阶偏导:
			  d^n f/dx ^n = lim_(x->0)  [ f^(n - 1)(x+Δx) - f^(n - 1)(x-Δx) ] / (2 Δx)
*	[算法]: 中心差分公式
		f'(x) = ( f(x+Δx) -  f(x-Δx) ) / (2 Δx) + Err(f,Δx)
		·截断误差: Err(f,Δx) = h² f^(3)(c) / 6 = O(h²)
		·精度: O(h²)
		f'(x) = ( -f(x+2Δx) + 8·f(x+Δx) - 8·f(x-Δx) + f(x-2Δx) ) / (12 Δx) + Err(f,Δx)
		·截断误差: Err(f,Δx) = h^4 f^(5)(c) / 6 = O(h^4)
		·精度: O(h^4)
******************************************************************************/
template<typename F>
double diff(double x0, F&& f, int N = 1,double dx = 1E-3) {
	return N == 0 ? f(x0) :
		(diff(x0 + dx, f, N - 1) - diff(x0 - dx, f, N - 1)) / (2 * dx);
}
template<typename F>
double diff_(double x0, F&& f, int N = 1, double dx = 1E-3) {
	return N == 0 ? f(x0) : 
		(
		- diff_(x0 + 2 * dx, f, N - 1) 
		+ diff_(x0 +	 dx, f, N - 1) * 8
		- diff_(x0 -	 dx, f, N - 1) * 8
		+ diff_(x0 - 2 * dx, f, N - 1)
		) / (12 * dx);
}
/******************************************************************************
*                    曲率
*	[定义]: 单位弧段弯曲的角度. 也即等效曲率圆的半径的倒数.
*	[公式]: K = |Δα/Δs| = 1 / R = |y''| / (1 + y'²)^(3/2)
******************************************************************************/
template<typename F>
double Curvature(double x0, F&& y, double dx = 1E-3) {
	return fabs(diff(x0, y, 2)) / pow(1 + pow(diff(x0, y), 2), 1.5);
}
/******************************************************************************
*                    偏导数
*	[定义]:
		偏导: ∂f/∂x_i = [ f(..,xi+Δxi,..) -  f(..,xi-Δxi,..) ] / (2 Δxi)
		二阶偏导:
			∂²f/∂x_i² = [ f'(..,xi+Δxi,..) -  f'(..,xi-Δ,..) ] / Δxi
					  = f(..,xi+Δxi) - 2·f'(..,x) + f'(..,xi-Δxi) / Δxi²
******************************************************************************/
template<typename F>
double PartiDeriv(Mat<>& x, int dim, double dx, F&& f) {
	x[dim] += dx;	double t1 = f(x);	x[dim] -= dx;
	x[dim] -= dx;	double t2 = f(x);	x[dim] += dx;
	return (t1 - t2) / (2 * dx);
}
template<typename F>
double PartiDeriv2(Mat<>& x, int dim, double dx, F&& f) {
	double ans = 0;
	x[dim] += dx; ans += f(x); x[dim] -= dx;
	x[dim] -= dx; ans += f(x); x[dim] += dx;
				  ans -= f(x) * 2;
	return ans / (dx * dx);
}
/******************************************************************************
*                    Hamilton 算子、Laplace 算子
*	[定义]: 
		Hamilton: ▽  ≡ ∂/∂x \vec x + ∂/∂y \vec y + ∂/∂z \vec z + ... (看成一个矢量)
		Laplace:  ▽² ≡ ∂²/∂x² + ∂²/∂y² + ∂²/∂z² + ...
		▽² = ▽·▽
*	[算法]: 有限差分法
		[u(x+1,...) - 2·u(x,...) + u(x-1,...)] / Δx² + ...
*	[注]: Hamilton算子程序，见<梯度、散度、旋度>
******************************************************************************/
template<typename F>
inline double LaplaceOperator(Mat<>& x, Mat<>& dx, F&& f) {
	double ans = 0;
	for (int dim = 0; dim < dx.size(); dim++) ans += PartiDeriv2(x, dim, dx[dim], f);
	return ans;
}
/******************************************************************************
*                    梯度、散度、旋度
*	[定义]:
		梯度: ▽f		: 矢量, 函数在该点处变化率最大的方向.
		散度: ▽·\vec f: 标量, 矢量场在该点发散的程度, 表征场的有源性(>0源,<0汇,=0无源)
		旋度: ▽×\vec f: 矢量, 矢量场在该点旋转的程度, 
							方向是旋转度最大的环量的旋转轴, 旋转的方向满足右手定则,
							大小是绕该旋转轴旋转的环量与旋转路径围成的面元面积之比.
*	[公式]: (直角坐标系)
		▽f        = ∂f/∂x \vec x + ∂f/∂y \vec y + ∂f/∂z \vec z + ...
		▽·\vec f =  ∂fx/∂x + ∂fy/∂y + ∂fz/∂z + ...
		▽×\vec f = (∂fz/∂y - ∂fy/∂z) \vec x
				   + (∂fx/∂z - ∂fz/∂x) \vec y
				   + (∂fy/∂x - ∂fx/∂y) \vec z
******************************************************************************/
template<typename F>
Mat<>& Grad(Mat<>& x, Mat<>& dx, F&& f, Mat<>& ans) {
	ans.alloc(x.rows);
	for (int dim = 0; dim < x.size(); dim++) ans[dim] = PartiDeriv(x, dim, dx[dim], f);
	return ans;
}
template<typename FX, typename FY, typename FZ>
double Div(Mat<>& x, Mat<>& dx, FX&& fx, FY&& fy, FZ&& fz) {
	double ans = 0;
	ans += PartiDeriv(x, 0, dx[0], fx);
	ans += PartiDeriv(x, 1, dx[1], fy);
	ans += PartiDeriv(x, 2, dx[2], fz);
	return ans;
}
template<typename FX, typename FY, typename FZ>
Mat<>& Curl(Mat<>& x, Mat<>& dx, FX&& f0, FY&& f1, FZ&& f2, Mat<>& ans) {
	ans.alloc(x.rows);
	ans[0] = PartiDeriv(x, 1, dx[1], f2) - PartiDeriv(x, 2, dx[2], f1);
	ans[1] = PartiDeriv(x, 2, dx[2], f0) - PartiDeriv(x, 0, dx[0], f2);
	ans[2] = PartiDeriv(x, 0, dx[0], f1) - PartiDeriv(x, 1, dx[1], f0);
	return ans;
}
/*#############################################################################

								级数展开

#############################################################################*/
/******************************************************************************
*                    Taylor 展开
*	[定义]: 
		f(x) = f(x0)/0! + f'(x0)/1!·(x-x0) + ... + f^(n)(x0)/n!·(x-x0)
*	[Example]:
		Calculus::TaylorFormula(0, [](double x) { return sin(x); }, Coeff, 10);
******************************************************************************/
template<typename F>
Mat<>& TaylorFormula(double x0, F&& f, Mat<>& Coeff, int N = 3) {
	Coeff.alloc(N + 1);
	for (int i = 0; i <= N; i++) 
		Coeff[i] = diff(x0, f, i) / NumberTheory::Factorial(i);
	return Coeff;
}
/******************************************************************************
*                    常用函数
*	[公式]: (Taylor展开)
		exp(x) = 1 + x + x^2/2! + x^3/3! + ...
		sin(x) = x - x^3/3! + x^5/5! - x^7/7! + ...
		cos(x) = 1 - x^2/2! + x^4/4! - x^6/6! + ...
		ln(1+x)= x - x^2/2  + x^3/3  - x^4/4  + ...	x∈[-1,1]
	  arctan(x)= x - x^3/3  + x^5/5  - x^7/7  + ...	x∈[-1,1]
		(1+x)^p= 1 + px + p(p-1)/2!·x^2 + p(p-1)(p-2)/3!·x^2 + ...	|x| < 1
******************************************************************************/
double Exp(double x, int N = 18) {
	double ans = 0;
	for (int i = 0; i <= N; i++) ans += pow(x, i) / NumberTheory::Factorial(i);
	return ans;
}
double Sin(double x, int N = 18) {
	double ans = 0;
	for (int i = 1; i <= N; i += 2) ans += (i % 4 == 1 ? 1 : -1) * pow(x, i) / NumberTheory::Factorial(i);
	return ans;
}
double Cos(double x, int N = 18) {
	double ans = 0;
	for (int i = 0; i <= N; i += 2) ans += (i % 4 == 0 ? 1 : -1) * pow(x, i) / NumberTheory::Factorial(i);
	return ans;
}
double lnOneAdd(double x, int N = 18) {
	double ans = 0;
	for (int i = 1; i <= N; i++) ans += (i % 2 == 1 ? 1 : -1) * pow(x, i) / i;
	return ans;
}
double Arctan(double x, int N = 18) {
	double ans = 0;
	for (int i = 1; i <= N; i += 2) ans += (i % 4 == 1 ? 1 : -1) * pow(x, i) / i;
	return ans;
}
double PowOneAdd(double x, double p, int N = 18) {
	double ans = 0, pTmp = 1;
	for (int i = 0; i <= N; i++) { ans += pTmp * pow(x, i) / NumberTheory::Factorial(i); pTmp *= (p - i); }
	return ans;
}
/******************************************************************************
*					Fast Fourier Transform 快速Fourier变换
*	[定义]: 离散Fourier变换的高效算法
*	[公式]:
		离散Fourier变换: X[k] = Σ_(n=0)^(N-1)  e^(-j2πnk/N)·x[n]
*	[时间复杂度]: O(N·lgN)
*	[Reference]:
		Thanks for https://www.math.wustl.edu/~victor/mfmm/fourier/fft.c
******************************************************************************/
/*
void FFT(std::complex<double>* v, int n, std::complex<double>* tmp) {
	if (n > 1) return;
	std::complex<double> z, w, * vo = tmp + n / 2, * ve = tmp;
	for (int k = 0; k < n / 2; k++) {
		ve[k] = v[2 * k];
		vo[k] = v[2 * k + 1];
	}
	FFT(ve, n / 2, v);							//FFT 偶数序列 v[] 
	FFT(vo, n / 2, v);							//FFT 奇数序列 v[] 
	for (int m = 0; m < n / 2; m++) {
		w.real = cos(2 * PI * m / (double)n);
		w.imag =-sin(2 * PI * m / (double)n);
		z			 = w * vo[m];
		v[m]		 = ve[m] + z;
		v[m + n / 2] = ve[m] - z;
	}
}
//逆Fourier变换
void iFFT(std::complex<double>* v, int n, std::complex<double>* tmp){
	if (n > 1) return;
	std::complex<double> z, w, * vo = tmp + n / 2, * ve = tmp;
	for (int k = 0; k < n / 2; k++) {
		ve[k] = v[2 * k];
		vo[k] = v[2 * k + 1];
	}
	iFFT(ve, n / 2, v);							//FFT 偶数序列 v[]
	iFFT(vo, n / 2, v);							//FFT 奇数序列 v[]
	for (int m = 0; m < n / 2; m++) {
		w.real = cos(2 * PI * m / (double)n);
		w.imag = sin(2 * PI * m / (double)n);
		z			 = w * vo[m];
		v[m]		 = ve[m] + z;
		v[m + n / 2] = ve[m] - z;
	}
}*/
/*#############################################################################

								积分

#############################################################################*/
/******************************************************************************
*                    积分
*	[定义]:
*	[算法]: NewtonCotes 公式
		∫_a^b f(x) = (b - a) Σ_(k=0)^n  C_k^(n) f(xi)
		C_k^(n) = (-1)^(n-k) / (n·k!(n-k)!) ∫_0^n Π_(k≠j) (t-j)dt 
		n = 1: C = {1/2, 1/2}
		n = 2: C = {1/6, 4/6, 1/6}
		n = 4: C = {7/90, 32/90, 12/90, 32/90, 7/90}
		* NewtonCotes 公式在 n > 8 时不具有稳定性
		复合求积法: 将积分区间分成若干个子区间, 再在每个子区间使用低阶求积公式.
******************************************************************************/
template<typename F>
double integral_NewtonCotes(double xSt, double xEd, F&& f, int n = 4) {
	double ans = 0, dx = (xSt - xEd) / n, xi = xSt,
		C[] = { 7 / 90.0, 32 / 90.0, 12 / 90.0, 32 / 90.0, 7 / 90.0 };
	for (int i = 0; i <= n; i++, xi += dx) ans += C[i] * f(xi);
	return ans *= (xEd - xSt);
}
template<typename F>
double integral(double xSt, double xEd, F&& f, int n) {
	double ans = 0, dx = (xSt - xEd) / n, xi = xSt;
	for (int i = 0; i < n; i++, xi += dx) ans += integral_NewtonCotes(xi, xi + dx, f);
	return ans;
}
/******************************************************************************
*                    重积分
*	[定义]: ∫∫∫ f(r) dr³
******************************************************************************/
template<typename F>
double multIntegral(Mat<>& dx, Mat<>&St, Mat<>& Ed, F&& f) {
	Mat<> x = St;
	double ans = 0, dx_n = dx.product();
	while (true) {
		int dim = 0; x[dim] += dx[dim];
		while (dim <  dx.size() - 1 && x[dim] > Ed[dim]) { x[dim] = St[dim]; dim++; x[dim] += dx[dim]; }
		if    (x[dx.size() - 1] > Ed[dx.size() - 1]) break;
		ans += f(x) * dx_n;
	}
	return ans;
}
/*#############################################################################

								微分方程

#############################################################################*/
/******************************************************************************
*                    解常微分方程组: Runge Kutta 方法
*	[公式]:           ->   ->       ->      ->
		常微分方程组: y' = f(x, y)	y(x0) = y0
		迭代求解 \vec y(x) 的一点/一区间的数值解.
		y(x + dx) = y(x) + dx/6·(k1 + 2·k2 + 2·k3 + k4)
		k1 = f(xn , yn)						//区间开始斜率
		k2 = f(xn + dx/2, yn + dx/2·k1)	//区间中点斜率,通过欧拉法采用k1决定y在xn+dx/2值
		k3 = f(xn + dx/2, yn + dx/2·k2)	//区间中点斜率,采用k2决定y值
		k4 = f(xn + dx	, yn + dx  ·k3)	//区间终点斜率
*	[目的]: 解常微分方程组
		[ y1'(x) = f1(y1 , ... , yn , x)          ->   ->
		| y2'(x) = f2(y1 , ... , yn , x)    =>    y' = f(x , y)
		| ...
		[ yn'(x) = fn(y1 , ... , yn , x)
*	[性质]:
		* RK4法是四阶方法，每步误差是h⁵阶，总积累误差为h⁴阶
******************************************************************************/
template<typename F>
void RungeKutta(Mat<>& y, double dx, double x0, F&& f, int enpoch = 1) {
	Mat<> tmp, k, k1, k2, k3, k4;
	double x = x0;
	while (enpoch--) {
		// k1, k2, k3 ,k4
		k1 = f(x, y);
		k2 = f(x + dx / 2, tmp.add(y, tmp.mul(dx / 2, k1)));
		k3 = f(x + dx / 2, tmp.add(y, tmp.mul(dx / 2, k2)));
		k4 = f(x + dx,     tmp.add(y, tmp.mul(dx,	   k3)));
		// y[n+1] = y[n] + h/6·(k1 + 2·k2 + 2·k3 + k4)
		k.add(k.add(k1, k4), tmp.mul(2, tmp.add(k2, k3)));
		y.add(y, k.mul(dx / 6, k));
	};
}
/******************************************************************************
*                    Poisson's方程
*	[定义]: Δφ = f
			▽²φ = f  (Euclidean空间)
		三维直角坐标系中 (∂²/∂x² + ∂²/∂y² + ∂²/∂z²) φ(x,y,z) = f(x,y,z)
		当f ≡ 0, 得到 Laplace's方程
	[解法]:  Green's函数  φ(r) = - ∫∫∫ f(rt) / 4π|r-rt| d³rt    ,r rt为矢量
	[唯一性定理]:
		对于各种边界条件，Poisson's方程可能有许多种解，但每个解的梯度相同.
		静电场情况下, 意味着在边界条件下的满足Poisson's方程的势函数，所解得的电场唯一确定.
******************************************************************************/
template<typename F>
double PoissonEquation(Mat<>& x, Mat<>& dx, Mat<>& St, Mat<>& Ed, F&& f) {
	Mat<> tmp; 
	return multIntegral(dx, St, Ed, [&x, &f, &tmp](Mat<>& xt) {
		return x == xt ? 0.0 : f(xt) / tmp.sub(x, xt).norm();
	}) / (4 * PI);
}
/******************************************************************************
*                    波动方程
*	[定义]: a ▽²u = ∂²u/∂t²
*	[算法]: 有限差分法
		u(t+1,...) = 2·u(t,...) - u(t-1,...)
				+ Δt²·a{[u(x+1,...) - 2·u(x,...) + u(x-1,...)]/Δx² + ...}
*                    扩散方程
*	[定义]: a ▽²u = ∂u/∂t
	[算法]: 有限差分法
		u(t+1,...) = u(t,r)
				+ Δt·a{[u(x+1,...) - 2·u(x,...) + u(x-1,...)]/Δx² + ...}
-------------------------------------------------------------------------------
*	[Example]: //波动方程
		Mat<> u(N, N), u_pre(N, N), x(2), dx(2);  dx.fill(2);
		for (int i = 0; i < u.size(); i++) {
			x.getData(u.i2x(i) * dx[0], u.i2y(i) * dx[1]);
			u_pre[i] = Calculus::WaveEquation( x, dx, 0, 1, 1,
				[&u, &u_pre, &dx](Mat<>& pos, double t = 0) {
					int x = pos[0] / dx[0],
						y = pos[1] / dx[1];
					if (x < 0 || y < 0 || x >= N || y >= N) return 0.0;
					return t == 0 ? u(x, y) : u_pre(x, y);
				}
			);
		}
		u.swap(u_pre); u(N / 2, N / 2) = sin(t / (2 * PI));
******************************************************************************/
template<typename F>
inline double WaveEquation		(Mat<>& x, Mat<>& dx, double t, double dt, double A, F&& u) {
	return 2 * u(x, t) - u(x, t - 1) + A * dt * dt * LaplaceOperator(x, dx, u);
}
template<typename F>
inline double DiffusionEquation	(Mat<>& x, Mat<>& dx, double t, double dt, double A, F&& u) {
	return u(x, t) + A * dt * LaplaceOperator(x, dx, u);
}
/*#############################################################################

								插值拟合

#############################################################################*/
/******************************************************************************
*					样条插值
*	[算法]: 过求解三弯矩方程组得出曲线函数组的过程
******************************************************************************/
void CubicSpline(double x, double y) {

}
}
#endif