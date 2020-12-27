#ifndef UNIVERSE_H
#define UNIVERSE_H
#include <stdlib.h>
#include <stdio.h>
#include "../LiGu_AlgorithmLib/Mat.h"
#include "../LiGu_Graphics/Plot.h"
#include "PhysicsBase.h"
/******************************************************************************
*                    Star	星球
******************************************************************************/
class Star : public RigidBody {
public:
	RGB color;
	double R = 0;
};
/******************************************************************************
*                    Universe
******************************************************************************/
class Universe
{
public:
	/*---------------- 物理常数 ----------------*/
	const double pi = 3.141592653589,					//圆周率
		G = 6.67408E-11;								//万有引力常数
	/*---------------- ���� ----------------*/
	static const int MaxStars = 200;
	double dt = 1;										// 时间分辨率
	Star Stars[MaxStars];
	int StarNum = 0;
	Plot plot;
	int times = 0;
	/*----------------[  ]----------------*/
	Universe() {
		plot.setAxisRange(-1E12, -1E12, 1E12, 1E12);
		InitStars();
		Star* masscentre = MassCentre(Stars, StarNum);
		for (int i = 0; i < StarNum; i++)
			CoordinateTrans(Stars[i], *masscentre, Stars[i]);
	}
	void InitStars();
	void Simulation() {
		GravitationSimulation();
		play();
		if (times++ > 3600) { DrawIt(); times = 0; }
	}
	void DrawIt() {
		for (int i = 0; i < StarNum; i++) {
			plot.g->PaintColor = Stars[i].color;
			plot.plotCircle(Stars[i].r[0], Stars[i].r[1], Stars[i].R);
		}
	}
	void play() {
		for (int i = 0; i < StarNum; i++) {
			for (int dim = 0; dim < 3; dim++) {
				Stars[i].v[dim] += Stars[i].a[dim] * dt;//反一下为什么不行？
				Stars[i].r[dim] += Stars[i].v[dim] * dt;
			}
		}
	}
	void GravitationSimulation() {
		for (int i = 0; i < StarNum; i++) {				//更新每个星球的加速度
			Stars[i].a.clean();
			double tmass = Stars[i].mass; Stars[i].mass = 0;
			Gravitation(Stars, StarNum, Stars[i].r, Stars[i].a);
			Stars[i].mass = tmass;
		}
	}
	/*----------------[ Gravitation 万有引力 ]----------------
	*	*万有引力:
	*	-> 	       m1 * m2    ^
	*	F = - G * ---------- r12
	*	            r12^2
	** ---------------------------------------- */
	void Gravitation(Star point[], int n, Vector& position, Vector& acceleration){
		for (int i = 0; i < n; i++) {
			if (point[i].mass == 0)continue;
			double s2 = 0;								//距离s平方
			for (int dim = 0; dim < 3; dim++)
				s2 += (point[i].r[dim] - position[dim]) * (point[i].r[dim] - position[dim]);
			double a = G * point[i].mass / s2;			//stars对position的加速度|a|
			for (int dim = 0; dim < 3; dim++)			// 分加速度ai
				acceleration[dim] += a * (point[i].r[dim] - position[dim]) / sqrt(s2);
		}
	}
	/*---------------- MassCentre 求质心 ----------------
	*	Mc = ΣMi
	*	Mc * rc = Σ Mi * ri
	*	Mc * vc = Σ Mi * vi
	*	Mc * ac = Σ Mi * ai
	** ---------------------------------------- */
	Star* MassCentre(Star point[], int n) {
		Star* ans = new Star;
		for (int i = 0; i < n; i++) {					//对于每个质点
			for (int j = 0; j < 3; j++) {		//对于每个维度
				ans->r[j] += point[i].r[j] * point[i].mass;
				ans->v[j] += point[i].v[j] * point[i].mass;
				ans->a[j] += point[i].a[j] * point[i].mass;
			}
			ans->mass += point[i].mass;
		}
		ans->r.mult(1.0 / ans->mass, ans->r);
		ans->v.mult(1.0 / ans->mass, ans->v);
		ans->a.mult(1.0 / ans->mass, ans->a);
		return ans;
	}
	/*---------------- 伽利略坐标变换 ----------------*/
	void CoordinateTrans(Star& point, Star& coord, Star& ans) {
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
};
#endif