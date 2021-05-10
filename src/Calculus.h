#ifndef CALCULUS_H
#define CALCULUS_H
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/Tensor.h"
#include "../../LiGu_AlgorithmLib/NumberTheory.h"
/******************************************************************************
*                    微积分 / 微分方程
-------------------------------------------------------------------------------
double diff		(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);		\\导数 (N阶)
double diff_	(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);
double Curvature(double x0, double(*y)(double x),			 double dx = 1E-3);		\\曲率
double PartiDeriv	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数
double PartiDeriv2	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数 (2阶)
Mat<>& TaylorFormula(double x0, double(*f)(double x), Mat<>& Coeff, int N = 3);		\\Taylor展开
double Exp		(double x, int N = 18);												\\常用函数
double Sin		(double x, int N = 18);
double Cos		(double x, int N = 18);
double lnOneAdd	(double x, int N = 18);
double Arctan	(double x, int N = 18);
double PowOneAdd(double x, double p, int N = 18);
void   RungeKutta(Mat<>& y, double dt, double t0, int enpoch, Mat<>& (*derivY)(double t, Mat<>& y));
																					\\解常微分方程组: RungeKutta法
	   PoissonEquation();															\\Poisson's方程
	   WaveEquation();																\\波动方程
******************************************************************************/
namespace Calculus {
#define PI 3.141592653589
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
double diff(double x0, double(*f)(double x), int N = 1,double dx = 1E-3) {
	return N == 0 ? f(x0) :
		(diff(x0 + dx, f, N - 1) - diff(x0 - dx, f, N - 1)) / (2 * dx);
}
double diff_(double x0, double(*f)(double x), int N = 1, double dx = 1E-3) {
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
double Curvature(double x0, double(*y)(double x), double dx = 1E-3) {
	return fabs(diff(x0, y, 2)) / pow(1 + pow(diff(x0, y), 2), 1.5);
}
/******************************************************************************
*                    偏导数
*	[定义]:
		偏导: ∂f/∂x_i = [ f(..,xi+Δxi,..) -  f(..,xi-Δxi,..) ] / (2 Δxi)
		二阶偏导:
			∂²f/∂x_i² = [ f'(..,xi+Δxi,..) -  f'(..,xi,..) ] / Δxi
					  = f(..,xi+2Δxi) - f'(..,x) - f'(..,xi+Δxi) + f'(..,xi-Δxi)
******************************************************************************/
double PartiDeriv(Mat<>& x, int index, double dx, double(*f)(Mat<>& x)) {
	Mat<> xt = x;
	xt[index] += dx;		double t1 = f(xt);
	xt[index] -= 2 * dx;	double t2 = f(xt);
	return (t1 - t2) / (2 * dx);
}
double PartiDeriv2(Mat<>& x, int index, double dx, double(*f)(Mat<>& x)) {
	Mat<> xt = x;
	xt[index] += 2 * dx;	double t1 = f(xt);
	xt[index] -= 2 * dx;	double t2 = f(xt);
	xt[index] += dx;		double t3 = f(xt);
	xt[index] -= 2 * dx;	double t4 = f(xt);
	return (t1 - t2 - t3 + t4) / (2 * dx * dx);
}
/******************************************************************************
*                    Taylor 展开
*	[定义]: 
		f(x) = f(x0)/0! + f'(x0)/1!·(x-x0) + ... + f^(n)(x0)/n!·(x-x0)
*	[Example]:
		Calculus::TaylorFormula(0, [](double x) { return sin(x); }, Coeff, 10);
******************************************************************************/
Mat<>& TaylorFormula(double x0, double(*f)(double x), Mat<>& Coeff, int N = 3) {
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
		ln(1+x)= x - x^2/2 + x^3/3 - x^4/4 + ...	x∈[-1,1]
	  arctan(x)= x - x^3/3 + x^5/5 - x^7/7 + ...	x∈[-1,1]
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
double integral_NewtonCotes(double xSt, double xEd, double(*f)(double x), int n = 4) {
	double ans = 0, dx = (xSt - xEd) / n, xi = xSt,
		C[] = { 7 / 90.0, 32 / 90.0, 12 / 90.0, 32 / 90.0, 7 / 90.0 };
	for (int i = 0; i <= n; i++, xi += dx) ans += C[i] * f(xi);
	return ans *= (xEd - xSt);
}
double integral(double xSt, double xEd, double(*f)(double x), int n) {
	double ans = 0, dx = (xSt - xEd) / n, xi = xSt;
	for (int i = 0; i < n; i++, xi += dx) ans += integral_NewtonCotes(xi, xi + dx, f);
	return ans;
}
/******************************************************************************
*                    解常微分方程组: Runge Kutta 方法
*	[公式]:           ->   ->       ->      ->
		对于初值问题: y' = f(t, y)	y(t0) = y0
		y[n+1] = y[n] + h/6·(k1 + 2·k2 + 2·k3 + k4)
		k1 = f(tn , yn)					//区间开始斜率
		k2 = f(tn + h/2 , yn + h/2·k1)	//区间中点斜率,通过欧拉法采用k1决定y在tn+h/2值
		k3 = f(tn + h/2 , yn + h/2·k2)	//区间中点斜率,采用k2决定y值
		k4 = f(tn + h , yn + h·k3)		//区间终点斜率
*	[目的]: 解常微分方程组
		[ y1'(t) = f1(y1 , ... , yn , t)          ->   ->
		| y2'(t) = f2(y1 , ... , yn , t)    =>    y' = f(t , y)
		| ...
		[ yn'(t) = fn(y1 , ... , yn , t)
*	[性质]:
		* RK4法是四阶方法，每步误差是h⁵阶，总积累误差为h⁴阶
******************************************************************************/
void RungeKutta(Mat<>& y, double dt, double t0, int enpoch, Mat<>& (*derivY)(double t, Mat<>& y)) {
	Mat<> tmp, k, k1, k2, k3, k4;
	double t = t0;
	while (enpoch--) {
		// k1, k2, k3 ,k4
		k1 = derivY(t, y);
		k2 = derivY(t + dt / 2, tmp.add(y, tmp.mult(dt / 2, k1)));
		k3 = derivY(t + dt / 2, tmp.add(y, tmp.mult(dt / 2, k2)));
		k4 = derivY(t + dt,     tmp.add(y, tmp.mult(dt,		k3)));
		// y[n+1] = y[n] + h/6·(k1 + 2·k2 + 2·k3 + k4)
		k.add(k.add(k1, k4), tmp.mult(2, tmp.add(k2, k3)));
		y.add(y, k.mult(dt / 6, k));
	};
}
/******************************************************************************
*                    Poisson's方程
*	[定义]: Δφ = f
		流形属于Euclidean空间时, 写成 ▽²φ = f
		三维直角坐标系中 (∂²/∂x² + ∂²/∂y² + ∂²/∂z²) φ(x,y,z) = f(x,y,z)
		当f ≡ 0, 得到 Laplace's方程
	[解法]:  Green's函数  φ(r) = - ∫∫∫ f(rt) / 4π|r-rt| d³rt    ,r rt为矢量
	[唯一性定理]:
		对于各种边界条件，Poisson's方程可能有许多种解，但每个解的梯度相同.
		静电场情况下, 意味着在边界条件下的满足Poisson's方程的势函数，所解得的电场唯一确定.
******************************************************************************/
Tensor<double>* PoissonEquation(Mat<>st, Mat<>ed, Mat<> delta, double (*f) (Mat<>& x)) {
	// init
	Mat<> tmp;
	tmp.sub(ed, st);
	Mat<int> dim; dim.E(st.rows);
	for (int i = 0; i < st.rows; i++) dim[i] = (int)(tmp[i] / delta[i]);
	Tensor<double>* Map = new Tensor<double>(dim.rows, dim.data);
	// compute Green's function
	Mat<> r = st;
	for (int i = 0; i < Map->dim.product(); i++) {
		double t = 0;
		Mat<> rt = st;
		for (int j = 0; j < Map->dim.product(); j++) {
			for (int k = 0; k < rt.rows; k++) {
				rt[k] += delta[k]; 
				if(rt[k] >= ed[k]) rt[k] = st[k]; 
				else break;
			}
			t += f(rt) / tmp.sub(r, rt).norm();
		}
		(*Map)[i] = 1 / (4 * PI) * delta.product() * t;
		// update r
		for (int k = 0; k < r.rows; k++) {
			r[k] += delta[k];
			if (r[k] >= ed[k]) r[k] = st[k]; else break;
		}
	}
	return Map;
}
/******************************************************************************
*                    波动方程
*	[定义]: a ▽²u = ∂²u/∂t²
			// ∂²u/∂t², 即.加速度
			其中 ▽² ≡ ∂²/∂x² + ∂²/∂y² + ∂²/∂z² + ...
		边界条件:
			[1] u(r,t = 0) = f(r)		//初位置
			[2] ∂u(r,t = 0)/∂t = g(r)	//初速度
			[3]	u(a,t) = b				//边界值
	[算法]: 有限差分法
		u(t+1,...) = 2·u(t,...) - u(t-1,...)
				+ Δt²·a{[u(x+1,...) - 2·u(x,...) + u(x-1,...)]/Δx² + ...}
	[推导]:
		∂²u/∂t² = { [u(t+1,rt) - u(t,rt)]/Δt - [u(t,rt) - u(t-1,rt)]/Δt } / Δt
				= [u(t+1,rt) - 2·u(t,rt) + u(t-1,rt)] / Δt²
		▽²u同理, 代入波动方程,
		a { [u(x+1,...) - 2·u(x,...) + u(x-1,...)] / Δx² + ...}
		  = [u(t+1,...) - 2·u(t,...) + u(t-1,...)] / Δt²
		所以, 下一时刻有
		u(t+1,...) = 2·u(t,...) - u(t-1,...)
				+ Δt²·a{[u(x+1,...) - 2·u(x,...) + u(x-1,...)]/Δx² + ...}
		特别的, 对于t = 1 * Δt时刻, 有
		u(r,t=Δt) = u(r,t=0) + ∂u(r,t=0)/∂t·Δt + a▽²u·Δt²
	[*] 暂时只二维
******************************************************************************/
void WaveEquation(Mat<>& Map, Mat<>& veloc, void (*setBoundaryEquations) (Mat<>& x, int time),
	double alpha, double deltaTime, double deltaX, double deltaY, int EndTimes) {
	Mat<> MapNow (Map), 
		  MapPrev(Map);
	//对于t = 1 * Δt时刻, 有 u(r,t=Δt) = u(r,t=0) + ∂u(r,t=0)/∂t·Δt + a▽²u·Δt²
	for (int x = 1; x < Map.rows - 1; x++) {
		for (int y = 1; y < Map.cols - 1; y++) {
			MapNow(x, y) += veloc(x, y) * deltaTime + alpha * deltaTime * deltaTime * (
				(Map(x + 1, y) - 2 * Map(x, y) + Map(x - 1, y)) / (deltaX * deltaX)
			  + (Map(x, y + 1) - 2 * Map(x, y) + Map(x, y - 1)) / (deltaY * deltaY)
			);
		}
	}setBoundaryEquations(MapNow, 0);
	for (int time = 1; time < EndTimes; time++) {
		//u(t+1,...) = 2·u(t,...) - u(t-1,...)+Δt²·a{ [u(x + 1,...) - 2·u(x,...) + u(x - 1,...)] / Δx² + ... }
		for (int x = 1; x < Map.rows - 1; x++) {
			for (int y = 1; y < Map.cols - 1; y++) {
				MapPrev(x, y) = 2 * MapNow(x, y) - MapPrev(x, y) + alpha * deltaTime * deltaTime * (
					(MapNow(x + 1, y) - 2 * MapNow(x, y) + MapNow(x - 1, y)) / (deltaX * deltaX)
				  + (MapNow(x, y + 1) - 2 * MapNow(x, y) + MapNow(x, y - 1)) / (deltaY * deltaY)
				);
			}
		}
		MapPrev.swap(MapNow);
		setBoundaryEquations(MapNow, time);
	}
	Map = MapNow;
}
/******************************************************************************
*                    扩散方程
*	[定义]: a ▽²u = ∂u/∂t
	[算法]: 有限差分法
		u(t+1,...) = u(t,r)
				+ Δt·a{[u(x+1,...) - 2·u(x,...) + u(x-1,...)]/Δx² + ...}
******************************************************************************/
}
#endif