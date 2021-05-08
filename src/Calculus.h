#ifndef CALCULUS_H
#define CALCULUS_H
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "../../LiGu_AlgorithmLib/Tensor.h"
/******************************************************************************
*                    微积分 / 微分方程
******************************************************************************/
namespace Calculus {
#define PI 3.141592653589
/******************************************************************************
*                    Runge Kutta 方法
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