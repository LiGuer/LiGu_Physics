#include "Dynamics.h"
/******************************************************************************
*                    Dynamics	力学类
******************************************************************************/
/*---------------- 伽利略坐标变换 ----------------*/
void Dynamics::CoordinateTrans(Particle& point, Coordinate& coord, Particle& ans) {
	// 平动
	for (int dim = 0; dim < 3; dim++) {
		ans.r[dim] = point.r[dim] - coord.r[dim];
		ans.v[dim] = point.v[dim] - coord.v[dim];
		ans.a[dim] = point.a[dim] - coord.a[dim];
	}
	// 转动
	Mat<double> mat(3);
	mat[0 * mat.cols + 0] = cos(coord.Angle[0]); mat[0 * mat.cols + 1] = -1 * sin(coord.Angle[0]);
	mat[1 * mat.cols + 0] = sin(coord.Angle[0]); mat[1 * mat.cols + 1] = cos(coord.Angle[0]);
	ans.r.mult(mat, ans.r);
}
/*---------------- MassCentre 求质心 ---------------- 
*	Mc = ΣMi
*	Mc * rc = Σ Mi * ri
*	Mc * vc = Σ Mi * vi
*	Mc * ac = Σ Mi * ai
** ---------------------------------------- */
Particle* Dynamics::MassCentre(Particle* point[], int n){
	Particle* ans = new Particle;
	for (int i = 0; i < n; i++) {					//对于每个质点
		for (int j = 0; j < 3; j++) {		//对于每个维度
			ans->r[j] += point[i]->r[j] * point[i]->mass;
			ans->v[j] += point[i]->v[j] * point[i]->mass;
			ans->a[j] += point[i]->a[j] * point[i]->mass;
		}
		ans->mass += point[i]->mass;
	}
	ans->r.mult(1.0 / ans->mass, ans->r);
	ans->v.mult(1.0 / ans->mass, ans->v);
	ans->a.mult(1.0 / ans->mass, ans->a);
	return ans;
}
/*---------------- Gravitation 万有引力 ----------------
*	*万有引力:
*	-> 	       m1 * m2    ^
*	F = - G * ---------- r12
*	            r12^2
** ---------------------------------------- */
void Dynamics::Gravitation(Particle *point[], int n, Vector& position, Vector& acceleration)
{
	for (int i = 0; i < n; i++) {
		if (point[i]->mass == 0)continue;
		// 距离s
		double s2 = 0;								//距离s平方
		for (int dim = 0; dim < 3; dim++)
			s2 += (point[i]->r[dim] - position[dim]) * (point[i]->r[dim] - position[dim]);
		// 加速度大小|a|
		double a = G * point[i]->mass / s2;					//stars对position的加速度
		// 分加速度ai
		for (int dim = 0; dim < 3; dim++)
			acceleration[dim] += a * (point[i]->r[dim] - position[dim]) / sqrt(s2);
	}
}
/*---------------- Gravitate Potential 引力势 ----------------
*	*引力势:
*	  		     m1
*	U = -G * ----------
*	            r12
** ---------------------------------------- */
double Dynamics::GravitatePotential(Particle point[], int n, Vector& position)
{
	double U = 0;
	for (int i = 0; i < n; i++) {
		double r = 0;
		for (int dim = 0; dim < 3; dim++) {
			r += (point[i].r[dim] - position[dim]) * (point[i].r[dim] - position[dim]);
		}
		r = sqrt(r);
		U += -G * (point[i].mass / r);
	}
	return U;
}
/*---------------- Centrifugal Potential 离心势 ----------------
*	      0      2             1       2  2
*	U = ∫  Omega  * r dr = - --- Omega  R 
*	      R                    2
** ---------------------------------------- */
double Dynamics::CentrifugalPotential(double Omega, Vector& center, Vector& position)
{
	double r2 = 0;
	for (int dim = 0; dim < 3; dim++) {
		r2 += (center[dim] - position[dim]) * (center[dim] - position[dim]);
	}
	return  -Omega * Omega * r2 / 2;
}
/*---------------- 二体引力各轴长度 ----------------
A: 半长轴		C: 焦距
角动量 L = m r vθ
					 1     2    G M m
能量  E = Ek + Ep = --- m v  + -------
					 2            r
参量
** ---------------------------------------- */
void Dynamics::TwobodyGravitationAxis(Particle& PM, Particle& Pm, double& A, double& C) {
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
** ---------------------------------------- */
double Dynamics::TwobodyGravitationPeriod(Particle& M, Particle& m) {
	double A, C;
	TwobodyGravitationAxis(M, m, A, C);
	double Omega = sqrt(G * (M.mass + m.mass) / (A * A * A));		//旋转周期的平均角速度
	return Omega;
}