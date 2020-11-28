#include "Dynamics.h"
/******************************************************************************
*                    Dynamics
******************************************************************************/
/*---------------- MassCentre 求质心 ---------------- 
*	Mc = ΣMi
*	Mc * rc = Σ Mi * ri
*	Mc * vc = Σ Mi * vi
*	Mc * ac = Σ Mi * ai
** ---------------------------------------- */
Particle* Dynamics::MassCentre(Particle* point[], int n)
{
	Particle* ans = new Particle;
	for (int i = 0; i < n; i++) {					//对于每个质点
		for (int j = 0; j < 3; j++) {		//对于每个维度
			ans->r[j] += point[i]->r[j] * point[i]->mass;
			ans->v[j] += point[i]->v[j] * point[i]->mass;
			ans->a[j] += point[i]->a[j] * point[i]->mass;
		}
		ans->mass += point[i]->mass;
	}
	for (int j = 0; j < 3; j++) {			//对于每个维度
		ans->r[j] /= ans->mass;
		ans->v[j] /= ans->mass;
		ans->a[j] /= ans->mass;
	}
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
		/*------ 距离s ------*/
		double s2 = 0;								//距离s平方
		for (int dim = 0; dim < 3; dim++)
			s2 += (point[i]->r[dim] - position[dim]) * (point[i]->r[dim] - position[dim]);
		/*------ 加速度大小|a| ------*/
		double a = G * point[i]->mass / s2;					//stars对position的加速度
		/*------ 分加速度ai ------*/
		for (int dim = 0; dim < 3; dim++) {
			acceleration[dim] += a * (point[i]->r[dim] - position[dim]) / sqrt(s2);
		}
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