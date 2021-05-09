/*
Copyright 2020 LiGuer. All Rights Reserved.
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
	http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/
#include "../../LiGu_AlgorithmLib/Mat.h"
#include "Calculus.h"
#include <functional>
namespace Dynamics {
/*************************************************************************************************
*							运动
*	[算法]: Runge Kutta
*	[公式]:           ->   ->       ->      ->
		对于初值问题: y' = f(t, y)	y(t0) = y0
		y[n+1] = y[n] + dt/6·(k1 + 2·k2 + 2·k3 + k4)
		k1 = f(tn , yn)						//区间开始斜率
		k2 = f(tn + dt/2 , yn + dt/2·k1)	//区间中点斜率,通过欧拉法采用k1决定y在tn+dt/2值
		k3 = f(tn + dt/2 , yn + dt/2·k2)	//区间中点斜率,采用k2决定y值
		k4 = f(tn + dt ,   yn + dt  ·k3)	//区间终点斜率
		y[n+1] = y[n] + dt/6·(k1 + 2·k2 + 2·k3 + k4)
-------------------------------------------------------------------------------------------------
*	[Example]: //N连杆摆动
		Mat<> angle(10), angleVeloc(10); angle.fill(PI / 2);
		Dynamics::run(
			angle, angleVeloc, 0.01, 1,
			[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
				Dynamics::Lagrange(
					x, dx, ddx,
					[](Mat<>& dx) {return 1.0 / 8 * (dx.dot(dx)); },
					[](Mat<>&  x) {
						double sum = 0;
						for (int i = 0; i < x.size(); i++) sum += (x.size() - i) * cos(x[i]);
						return 10.0 / 2 * (sum);
					}
				);
			}
		);
*	[Example] //重力 + 抛体 + 弹簧
		Mat<> x(2), v(2); v.fill(40);
		Dynamics::run(
			x, v, 0.1, 1,
			[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
				Dynamics::Lagrange(
					x, dx, ddx,
					[](Mat<>& dx) { return 1.0 / 2 * dx.dot(dx); },
					[](Mat<>& x)  { return x[1] * 10 
					+ 0.05 * (pow(x[0], 2) + pow(x[1], 2))
					; }
				);
			}
		);
*	[Example] //弹性碰撞
		#define N 20
		Mat<> x(2 * N), v(2 * N);
		x.rands(2 * N, 1, -200, 400);
		v.rands(2 * N, 1, -50, 50);
		Dynamics::run(
			x, v, 0.1, 1,
			[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
				Dynamics::Lagrange(
					x, dx, ddx,
					[](Mat<>& dx) { return 1.0 / 2 * dx.dot(dx); },
					[](Mat<>& x) {
						double sum = 0;
						for (int i = 0; i < N; i++) {
							for (int j = i + 1; j < N; j++) {
								double d = sqrt(pow(x[2 * i] - x[2 * j], 2) + pow(x[2 * i + 1] - x[2 * j + 1], 2));
								sum += d > 40 ? 0 : 1000 * (40 - d);
							}
						}return sum;
					}
				);
			}
		);
*************************************************************************************************/
/*----------------[ 运动 ]----------------*/
void run(Mat<>& x, Mat<>&dx, double dt, int enpoch, void(*AcceleFun)(Mat<>& x, Mat<>& dx, Mat<>& ddx)) {
	Mat<>
		 xTmp	(x),
		dxTmp	(dx),
		ddx		(x.rows),
		kix, kidx, kx, kdx,tmp;
	while (enpoch--) {
		for (int k = 0; k < 4; k++) {
			AcceleFun(xTmp, dxTmp, ddx);
			kidx = ddx;
			kix  = dxTmp;
			if (k == 0) {
				kx  = kix;					 xTmp.add( x, kix  *= dt / 2);
				kdx = kidx;					dxTmp.add(dx, kidx *= dt / 2);
			}
			if (k == 1) {
				kx  += tmp.mult(2, kix );	 xTmp.add( x, kix  *= dt / 2);
				kdx += tmp.mult(2, kidx);	dxTmp.add(dx, kidx *= dt / 2);	
			}
			if (k == 2) {
				kx  += tmp.mult(2, kix);	 xTmp.add( x, kix  *= dt);
				kdx += tmp.mult(2, kidx);	dxTmp.add(dx, kidx *= dt);
			}
			if (k == 3) {
				kx  += kix ;
				kdx += kidx;
			}
		}
		 x += (kx  *= dt / 6);
		dx += (kdx *= dt / 6);
	}
}
void run(Mat<>* r, Mat<>* v, double* mass, int N, double dt, int enpoch, 
	void(*AcceleFun)(Mat<>*, double*, Mat<>*, int) // r,mess,accele,N
) {
	static Mat<>* rt =		(Mat<>*)calloc(N, sizeof(Mat<>)),
				* vt =		(Mat<>*)calloc(N, sizeof(Mat<>)),
				* kr =		(Mat<>*)calloc(N, sizeof(Mat<>)),
				* kv =		(Mat<>*)calloc(N, sizeof(Mat<>)),
				* accele =	(Mat<>*)calloc(N, sizeof(Mat<>)),
				* kiv, * kir, tmp;
	for (int i = 0; i < N; i++) accele[i].zero(r[0]);
	while (enpoch--) {
		// k1,k2,k3,k4		//r_diff = v;  v_diff = a
		for (int k = 0; k < 4; k++) {
			if (k == 0) {
				kir = v;  AcceleFun(r,  mass, accele, N);
			}
			else {
				kir = vt; AcceleFun(rt, mass, accele, N);
			}
			kiv = accele;
			for (int i = 0; i < N; i++) {
				if (k == 0) {
					kr[i] = kir[i];
					kv[i] = kiv[i];
				}
				else if (k < 3) {
					kr[i] += tmp.mult(2, kir[i]);
					kv[i] += tmp.mult(2, kiv[i]);
				}
				else {
					kr[i] += kir[i];
					kv[i] += kiv[i];
				}
				if (k < 3) {
					rt[i].add(r[i], tmp.mult(dt / 2, kir[i]));
					vt[i].add(v[i], tmp.mult(dt / 2, kiv[i]));
				}
			}
		}
		for (int i = 0; i < N; i++) {
			r[i] += (kr[i] *= dt / 6);
			v[i] += (kv[i] *= dt / 6);
		}
	}
}
/*************************************************************************************************
*							Lagrange 方程
*	[公式]: \frac{d}{dt}(\frac{∂L}{∂\dot x_i}) = \frac{∂L}{∂x_i}
*			L = T - U
*	x: 广义坐标
*************************************************************************************************/
Mat<>& Lagrange(Mat<>& x, Mat<>& dx, Mat<>& ddx, double(*T)(Mat<>& dx), double(*U)(Mat<>& x)) {
	for (int i = 0; i < x.size(); i++)
		ddx[i] = - Calculus::PartiDeriv(x, i, 1E-9, U) / Calculus::PartiDeriv2(dx, i, 1E-4, T);
	return ddx;
}
/*************************************************************************************************
*							MassCentre	质心
*	[定义]:
		[1] Mc = ΣMi
		[2] Mc * rc = Σ Mi * ri
		[3] Mc * vc = Σ Mi * vi
*************************************************************************************************/
double MassCentre(Mat<>* r, Mat<>* v, double* m, int N, Mat<>& rc, Mat<>& vc) {
	rc.zero(r[0]); 
	vc.zero(v[0]);
	Mat<> tmp;
	double mc = 0;
	for (int i = 0; i < N; i++) {
		rc += tmp.mult(m[i], r[i]);
		vc += tmp.mult(m[i], v[i]);
		mc += m[i];
	}
	rc *= 1 / mc;
	vc *= 1 / mc;
	return mc;
}
/*************************************************************************************************
*							Centrifugal Potential 离心势
*	      0   2             1    2  2
*	U = ∫  ω  * r dr = - --- ω  R
*	      R                 2
*************************************************************************************************/
double CentrifugalPotential(double angularVeloc, Mat<>& center, Mat<>& position) {
	Mat<> tmp;
	return  -pow(angularVeloc * tmp.sub(center, position).norm(), 2) / 2;
}
/*************************************************************************************************
*							Gravitation	万有引力
*	[定义]: -> 	       m1·m2  ^
			F = - G·-------- r12
						r12²
*	[微分方程]: 
		方程数 = Point Num × (3 + 3)			//令 Ai = - G·Σj mj / rij²// F/m
		xi' = (xi')  
		yi' = (yi')  
		zi' = (zi')
		(xi')' = Ai_x  
		(yi')' = Ai_y  
		(zi')' = Ai_z
*	[引力势]: U = - G * m1 / r12      
*************************************************************************************************/
static const double G = 6.67408E-11;
/*----------------[ 万有引力 ]----------------*/
void GravitationAcceleration(Mat<>* r, double* mass, Mat<>* accele, int N) {
	Mat<> dr;
	for (int i = 0; i < N; i++) {
		accele[i].zero(r[0]);
		for (int j = 0; j < N; j++) {
			if (j == i || mass[j] == 0) continue;
			accele[i] += (dr *= (-G * mass[j] / pow(dr.sub(r[i], r[j]).norm(), 3)));
		}
	}
}
/*----------------[ 引力势 ]----------------*/
double* GravitatePotential(Mat<>* r, double* mass, int N) {
	static const double G = 6.67408E-11;
	double* U = (double*)calloc(N, sizeof(double));
	Mat<> dr;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j == i || mass[j] == 0) continue;
			U[i] = -G * mass[j] / dr.sub(r[i], r[j]).norm();
		}
	}return U;
}
/*---------------- 二体引力各轴长度 ----------------
*	A: 半长轴		C: 焦距
*	角动量 L = m r vθ
*	能量  E = Ek + Ep = 1/2 m v²  + G M m / r
** ------------------------------------------------*/
void TwobodyGraviAxis(Mat<>& Mr, Mat<>& mr, Mat<>& Mv, Mat<>& mv, double M, double m, double& A, double& C) {
	double r2 = 0, v2 = 0;									//r: PM Pm距离//v: PM Pm相对速度
	for (int dim = 0; dim < Mr.size(); dim++) {
		r2 += pow(Mr[dim] - mr[dim], 2);
		v2 += pow(Mv[dim] - mv[dim], 2);
	}
	double 
		r		= sqrt(r2), 
		vtheta	= sqrt(v2),
		L		= m * r * vtheta,							//L:角动量
		E		= 1.0 / 2 * m * v2 - G * M * m / r,			//E:能量
		p		= (r2 * vtheta * vtheta) / (G * M),
		epsi2	= 1 + (2 * E * p) / (G * M * m);
	A = p / (1 - epsi2);
	C = sqrt(epsi2) * A;
}
/*---------------- 二体引力旋转周期 ----------------
考虑m对M的引力影响，有
约化质量: μm = (M m) / (M +m)
Fm = μm am
即
	 (M + m) m   ^     ->
- G ------------ r = m am
	    r^2
			   3
T = 2π sqrt( A  / G (M + m) )
			 3    2                              2
开普勒定律: A  / T  = k ,  其中 k = (G M) / (4π ) * (1 + m / M)
** ------------------------------------------------*/
double TwobodyGraviPeriod(Mat<>& Mr, Mat<>& mr, Mat<>& Mv, Mat<>& mv, double M, double m) {
	double A, C;
	TwobodyGraviAxis(Mr, mr, Mv, mv, M, m, A, C);
	return sqrt(G * (M + m) / (A * A * A));		//旋转周期的平均角速度
} 
}