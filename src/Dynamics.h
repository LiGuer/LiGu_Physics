#ifndef DYNAMICS_H
#define DYNAMICS_H
#include <stdlib.h>
#include <math.h>
#include "PhysicsBase.h"
#include "Ligu_Graphics/Mat.h"
/******************************************************************************
*                    Coordinate	坐标系
******************************************************************************/
class Coordinate
{
public:
	Vector r = Vector(3, 1); Vector v = Vector(3, 1); Vector a = Vector(3, 1);
	Vector Angle = Vector(3, 1), AngularVeloc = Vector(3, 1), AngularAccele = Vector(3, 1);//角度//角速度//角加速度
	void play(double dt) {
		for (int dim = 0; dim < 3; dim++) {
			r[dim] += v[dim] * dt;
			v[dim] += a[dim] * dt;
			Angle[dim] += AngularVeloc[dim] * dt;
			AngularVeloc[dim] += AngularAccele[dim] * dt;
		}
	}
};
/******************************************************************************
*                    Dynamics	力学类
******************************************************************************/
class Dynamics : public PhysicsBase
{
public:
	/*---------------- 函数 ----------------*/
	void play(Particle& point, double dt) {
		for (int dim = 0; dim < 3; dim++) {
			point.v[dim] += point.a[dim] * dt;//反一下为什么不行？
			point.r[dim] += point.v[dim] * dt;
		}
	}
	void play(Particle* point[], int n, double dt) {
		for (int i = 0; i < n; i++) {
			for (int dim = 0; dim < 3; dim++) {
				point[i]->v[dim] += point[i]->a[dim] * dt;//反一下为什么不行？
				point[i]->r[dim] += point[i]->v[dim] * dt;
			}
		}
	}
	/*---------------- 伽利略坐标变换 ----------------*/
	void CoordinateTrans(Particle& point, Coordinate& coord, Particle& ans) {
		for (int dim = 0; dim < 3; dim++) {
			ans.r[dim] = point.r[dim] - coord.r[dim];
			ans.v[dim] = point.r[dim] - coord.v[dim];
			ans.a[dim] = point.r[dim] - coord.a[dim];
		}
		Mat<double> mat;
		mat.E(3);
		mat[0 * mat.cols + 0] = cos(coord.Angle[0]); mat[0 * mat.cols + 1] = -1 * sin(coord.Angle[0]);
		mat[1 * mat.cols + 0] = sin(coord.Angle[0]); mat[1 * mat.cols + 1] = cos(coord.Angle[0]);
		ans.r.mult(mat, ans.r);
	}
	Particle* MassCentre(Particle* point[], int n);					//求质心
	/*---------------- 二体引力各轴长度 ----------------
	A: 半长轴		C: 焦距
	角动量 L = m r vθ
	                     1     2    G M m
	能量  E = Ek + Ep = --- m v  + -------
	                     2            r
	参量	
	** ---------------------------------------- */
	void TwobodyGravitationAxis(Particle& PM, Particle& Pm, double& A, double& C) {
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
	double TwobodyGravitationPeriod(Particle& M, Particle& m) {
		double A, C;
		TwobodyGravitationAxis(M, m, A, C);
		double Omega = sqrt(G * (M.mass + m.mass) / (A * A * A));		//旋转周期的平均角速度
		return Omega;
	}
	/*---------------- 万有引力 ----------------*/
	void Gravitation(Particle* point[], int n, Vector& position, Vector& acceleration);	//万有引力
	double GravitatePotential(Particle point[], int n, Vector& position);	//引力势
	double CentrifugalPotential(double Omega, Vector& center, Vector& position);	//引力势
};
#endif