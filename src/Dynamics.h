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
#include "Mat.h"
/*************************************************************************************************
*							MassCentre	质心
*	[定义]:
		[1] Mc = ΣMi
		[2] Mc * rc = Σ Mi * ri
		[3] Mc * vc = Σ Mi * vi
*************************************************************************************************/
double MassCentre(Mat<double>* r, Mat<double>* v, double* m, int N, Mat<double>& rc, Mat<double>& vc) {
	rc.zero(r[0].rows, 1); vc.zero(v[0].rows, 1);
	Mat<double> tmp;
	double MassC = 0;
	for (int i = 0; i < N; i++) {
		rc.add(rc, tmp.mult(m[i], r[i]));
		vc.add(vc, tmp.mult(m[i], v[i]));
		MassC += m[i];
	}
	rc.mult(1 / MassC, rc);
	vc.mult(1 / MassC, vc);
	return MassC;
}
/*************************************************************************************************
*							CoordinateTrans	伽利略坐标变换
*************************************************************************************************/
void CoordinateTrans(Mat<double>* r, Mat<double>* v, Mat<double>* a,int N, Mat<double>& coord) {
	// 平动
	Mat<double> tmp;
	for (int i = 0; i < N; i++) {
		r.add(r, coord.negative(tmp));
		v.add(v, coord.negative(tmp));
		a.add(a, coord.negative(tmp));
	}
	// 转动
	Mat<double> mat(3);
	mat[0 * mat.cols + 0] = cos(coord.Angle[0]); mat[0 * mat.cols + 1] = -1 * sin(coord.Angle[0]);
	mat[1 * mat.cols + 0] = sin(coord.Angle[0]); mat[1 * mat.cols + 1] = cos(coord.Angle[0]);
	ans.r.mult(mat, ans.r);
}
/*************************************************************************************************
*							Centrifugal Potential 离心势
*	      0      2             1       2  2
*	U = ∫  Omega  * r dr = - --- Omega  R
*	      R                    2
*************************************************************************************************/
double CentrifugalPotential(double Omega, Mat<double>& center, Mat<double>& position) {
	Mat<double> tmp;
	double r = tmp.add(center, position.negative(tmp)).norm();
	return  -(Omega * r) * (Omega * r) / 2;
}
/*************************************************************************************************
*							Gravitation	万有引力
*	[定义]: -> 	       m1·m2  ^
			F = - G·-------- r12
						r12²
*	[微分方程]: 方程数 = Point Num × (3 + 3)
			[ xi' = (xi')  yi' = (yi')  zi' = (zi')
			|	令 Ai = - G·Σj mj / rij²				// F/m
			[ (xi')' = Ai_x  (yi')' = Ai_y  (zi')' = Ai_z
*	[算法]: Runge Kutta
*	[引力势]: U = - G * m1 / r12      
*************************************************************************************************/
/*----------------[ 万有引力 ]----------------*/
void GravitationAcceleration(Mat<double>* r, Mat<double>* v, double* mass, Mat<double>* accele, int N) {
	static const double G = 6.67408E-11;
	Mat<double> tmp, dr;
	for (int i = 0; i < N; i++) {
		accele[i].clean();
		for (int j = 0; j < N; j++) {
			if (j == i || mass[j] == 0)continue;
			dr.add(r[i], r[j].negative(tmp));
			double s = dr.norm(), accele_A = -G * mass[j] / (s * s);
			accele[i].add(accele[i], tmp.mult(accele_A / s, dr));
		}
	}
}
/*----------------[ 引力运动 ]----------------*/
void Gravitation(Mat<double>* r, Mat<double>* v, double* mass, int N, double dt, int enpoch) {
	static Mat<double>* rt = (Mat<double>*)calloc(N, sizeof(Mat<double>))
		, * vt = (Mat<double>*)calloc(N, sizeof(Mat<double>))
		, * kiv, * kir
		, * kr = (Mat<double>*)calloc(N, sizeof(Mat<double>))
		, * kv = (Mat<double>*)calloc(N, sizeof(Mat<double>))
		, * accele = (Mat<double>*)calloc(N, sizeof(Mat<double>))
		, tmp;
	for (int i = 0; i < N; i++)accele[i].zero(r[0].rows, 1);
	while (enpoch--) {
		// Runge Kutta  // k1,k2,k3,k4
		// k1  //r_diff = v;  v_diff = a
		kir = v;
		GravitationAcceleration(r, v, mass, accele, N);
		kiv = accele;
		for (int i = 0; i < N; i++) {
			kr[i] = kir[i]; kv[i] = kiv[i];
			rt[i].add(r[i], tmp.mult(dt / 2, kir[i]));
			vt[i].add(v[i], tmp.mult(dt / 2, kiv[i]));
		}
		// k2
		kir = vt;
		GravitationAcceleration(rt, vt, mass, accele, N);
		kiv = accele;
		for (int i = 0; i < N; i++) {
			kr[i].add(kr[i], tmp.mult(2, kir[i]));
			kv[i].add(kv[i], tmp.mult(2, kiv[i]));
			rt[i].add(r[i], tmp.mult(dt / 2, kir[i]));
			vt[i].add(v[i], tmp.mult(dt / 2, kiv[i]));
		}
		// k3
		kir = vt;
		GravitationAcceleration(rt, vt, mass, accele, N);
		kiv = accele;
		for (int i = 0; i < N; i++) {
			kr[i].add(kr[i], tmp.mult(2, kir[i]));
			kv[i].add(kv[i], tmp.mult(2, kiv[i]));
			rt[i].add(r[i], tmp.mult(dt, kir[i]));
			vt[i].add(v[i], tmp.mult(dt, kiv[i]));
		}
		// k4
		kir = vt;
		GravitationAcceleration(rt, vt, mass, accele, N);
		kiv = accele;
		for (int i = 0; i < N; i++) {
			kr[i].add(kr[i], kir[i]);
			kv[i].add(kv[i], kiv[i]);
			// y[n+1] = y[n] + h/6·(k1 + 2·k2 + 2·k3 + k4)
			r[i].add(r[i], kr[i].mult(dt / 6, kr[i]));
			v[i].add(v[i], kv[i].mult(dt / 6, kv[i]));
		}
	}
}
/*----------------[ 引力势 ]----------------*/
double* GravitatePotential(Mat<double>* r, double* mass, int N) {
	static const double G = 6.67408E-11;
	double* U = (double*)calloc(N, sizeof(double));
	Mat<double> dr;
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			if (j == i || mass[j] == 0)continue;
			dr.add(r[i], r[j].negative(tmp));
			U[i] = -G * mass[j] / dr.norm();
		}
	}return U;
}
/*************************************************************************************************

*************************************************************************************************/
/*---------------- 二体引力各轴长度 ----------------
*	A: 半长轴		C: 焦距
*	角动量 L = m r vθ
*	能量  E = Ek + Ep = 1/2 m v²  + G M m / r
** ---------------------------------------- 
void TwobodyGravitationAxis(Mat<double> PM, Particle& Pm, double& A, double& C) {
	double M = PM.mass, m = Pm.mass, r2 = 0, v2 = 0;						//r: PM Pm距离//v: PM Pm相对速度
	for (int dim = 0; dim < 3; dim++) {
		r2 += (PM.r[dim] - Pm.r[dim]) * (PM.r[dim] - Pm.r[dim]);
		v2 += (PM.v[dim] - Pm.v[dim]) * (PM.v[dim] - Pm.v[dim]);
	}
	double r = sqrt(r2), vtheta = sqrt(v2);					//bug:vtheta
	double L = m * r * vtheta;								//L:角动量
	double E = 1.0 / 2 * m * v2 - G * M * m / r;			//E:能量
	double p = (r2 * vtheta * vtheta) / (G * M);
	double epsi2 = 1 + (2 * E * p) / (G * M * m);
	A = p / (1 - epsi2);
	C = sqrt(epsi2) * A;
}*/
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
** ----------------------------------------
double TwobodyGravitationPeriod(Mat<double>& M, Mat<double>& m) {
	double A, C;
	TwobodyGravitationAxis(M, m, A, C);
	double Omega = sqrt(G * (M.mass + m.mass) / (A * A * A));		//旋转周期的平均角速度
	return Omega;
} */


/*
int main() {
	Graphics3D g(1000, 1000);
	Mat<double> tmp, tmp2;
	{ tmp.zero(3, 1); double a[] = { 5e8,5e8,5e8 }; tmp.getData(a); }
	g.setAxisLim(tmp.negative(tmp2), tmp);

	Mat<double>* r = (Mat<double>*)calloc(5, sizeof(Mat<double>));
	Mat<double>* v = (Mat<double>*)calloc(5, sizeof(Mat<double>));
	double* mass = (double*)calloc(5, sizeof(double));

	//---- 太阳 ----
	{ r[0].zero(3, 1); double t[] = { 0,0,0 }; r[0].getData(t); }
	{ v[0].zero(3, 1); double t[] = { 0,0,0 }; v[0].getData(t); }
	mass[0] = 0;//1.9891E30;
	//---- 地球 ----//1.496E11 29.783E3
	{ r[1].zero(3, 1); double t[] = { 0,0,0 }; r[1].getData(t); }
	{ v[1].zero(3, 1); double t[] = { 0,0,0 }; v[1].getData(t); }
	mass[1] = 5.965E24;
	//---- 月球 ----
	{ r[2].zero(3, 1); double t[] = { 3.844039E8,0,0 }; r[2].getData(t); }
	{ v[2].zero(3, 1); double t[] = { 0,1.01775E3,0 }; v[2].getData(t); }
	mass[2] = 7.349E22;
	//---- 木星 ----
	{ r[3].zero(3, 1); double t[] = { -8.165208E11,0,0 }; r[3].getData(t); }
	{ v[3].zero(3, 1); double t[] = { 0,-13069.7,0 }; v[3].getData(t); }
	mass[3] = 0; //1.90E27;
	//---- 人造卫星 ----
	{ r[4].zero(3, 1); double t[] = { 4E8,0,0 }; r[4].getData(t); }
	{ v[4].zero(3, 1); double t[] = { 0,0E3,0 }; v[4].getData(t); }
	mass[4] = 11e3;

	Mat<double>rc, vc;
	MassCentre(r, v, mass, 5, rc, vc);
	for (int i = 0; i < 5; i++) { r[i].add(r[i], rc.negative(tmp)); v[i].add(v[i], vc.negative(tmp)); }

	while (true) {
		Gravitation(r, v, mass, 5, 1, 1000);
		for (int i = 0; i < 5; i++)g.drawPoint(r[i]);
		g.g->PicWrite("D:/LIGU.ppm");
	}
}
*/
