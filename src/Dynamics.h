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
namespace Dynamics {
/*************************************************************************************************
*							运动
*	[算法]: Runge Kutta
		y[n+1] = y[n] + h/6·(k1 + 2·k2 + 2·k3 + k4)
*************************************************************************************************/
/*----------------[ 运动 ]----------------*/
void run(Mat<>* r, Mat<>* v, double* mass, int N, double dt, int enpoch, 
	void(*AcceleFun)(Mat<>*, Mat<>*, double*, Mat<>*, int) // r,v,mess,accele,N
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
				kir = v;  AcceleFun(r,  v,  mass, accele, N);
			}
			else {
				kir = vt; AcceleFun(rt, vt, mass, accele, N);
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
double CentrifugalPotential(double Omega, Mat<>& center, Mat<>& position) {
	Mat<> tmp;
	double r = tmp.sub(center, position).norm();
	return  -pow(Omega * r, 2) / 2;
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
void GravitationAcceleration(Mat<>* r, Mat<>* v, double* mass, Mat<>* accele, int N) {
	Mat<> dr;
	for (int i = 0; i < N; i++) {
		accele[i].zero();
		for (int j = 0; j < N; j++) {
			if (j == i || mass[j] == 0) continue;
			double s = dr.sub(r[i], r[j]).norm();
			accele[i] += (dr *= (-G * mass[j] / pow(s, 3)));
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
** ----------------------------------------*/
void TwobodyGravitationAxis(Mat<>& Mr, Mat<>& mr, Mat<>& Mv, Mat<>& mv, double M, double m, double& A, double& C) {
	double r2 = 0, v2 = 0;									//r: PM Pm距离//v: PM Pm相对速度
	for (int dim = 0; dim < 3; dim++) {
		r2 += pow(Mr[dim] - mr[dim], 2);
		v2 += pow(Mv[dim] - mv[dim], 2);
	}
	double r = sqrt(r2), vtheta = sqrt(v2);					//bug:vtheta
	double L = m * r * vtheta;								//L:角动量
	double E = 1.0 / 2 * m * v2 - G * M * m / r;			//E:能量
	double p = (r2 * vtheta * vtheta) / (G * M);
	double epsi2 = 1 + (2 * E * p) / (G * M * m);
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
** ----------------------------------------*/
double TwobodyGravitationPeriod(Mat<>& Mr, Mat<>& mr, Mat<>& Mv, Mat<>& mv, double M, double m) {
	double A, C;
	TwobodyGravitationAxis(Mr, mr, Mv, mv, M, m, A, C);
	return sqrt(G * (M + m) / (A * A * A));		//旋转周期的平均角速度
} 
}