#ifndef DYNAMICS_H
#define DYNAMICS_H
#include <stdlib.h>
#include <math.h>
#include "PhysicsBase.h"
#include "../LiGu_AlgorithmLib/Mat.h"
/******************************************************************************
*                    Coordinate	坐标系
******************************************************************************/
class Coordinate{
public:
	Vector r{ 3, 1 }, v{ 3, 1 }, a{ 3, 1 };					//位置 速度
	Vector Angle{ 3, 1 }, AngularVeloc{ 3, 1 }, AngularAccele{ 3, 1 };//角度//角速度//角加速度
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
class Dynamics : public PhysicsBase{
public:
	/*---------------- 函数 ----------------*/
	void play(Particle* point[], int n, double dt) {
		for (int i = 0; i < n; i++) {
			for (int dim = 0; dim < 3; dim++) {
				point[i]->v[dim] += point[i]->a[dim] * dt;//反一下为什么不行？
				point[i]->r[dim] += point[i]->v[dim] * dt;
			}
		}
	}
	void CoordinateTrans(Particle& point, Coordinate& coord, Particle& ans);		//伽利略坐标变换
	Particle* MassCentre(Particle* point[], int n);									//质心
	/*----------------[ 万有引力 ]----------------*/
	void Gravitation(Particle* point[], int n, Vector& position, Vector& acceleration);	//万有引力
	double GravitatePotential(Particle point[], int n, Vector& position);			//引力势
	void TwobodyGravitationAxis(Particle& PM, Particle& Pm, double& A, double& C);	//二体引力各轴长度
	double TwobodyGravitationPeriod(Particle& M, Particle& m);						//二体引力旋转周期

	double CentrifugalPotential(double Omega, Vector& center, Vector& position);	//离心势
};
#endif