#ifndef PHYSICSBASE_H
#define PHYSICSBASE_H
#include "Ligu_Graphics/Mat.h"
/******************************************************************************
*                    Particle	质点
******************************************************************************/
typedef Mat<double> Vector;
class Particle {
public:
	Vector r = Vector(3, 1), v = Vector(3, 1), a = Vector(3, 1);//位置//速度
	double mass = 0, eleccharge = 0;				//质量
	/*---------------- 动量 ----------------*/
	void Momentum(Vector& ans) {
		ans.mult(mass, v);
	}
};
/******************************************************************************
*                    RigidBody	刚体
******************************************************************************/
class RigidBody:public Particle
{
public:
	Vector Angle = Vector(3, 1), AngularVeloc = Vector(3, 1), AngularAccele = Vector(3, 1);//角度//角速度//角加速度
	Mat<double> Inertia = Mat<double>(3, 3);		//转动惯量
	/*---------------- 角动量 ----------------*/
	void AngularMomentum(Vector& ans) {
		ans.mult(Inertia, AngularVeloc);
	}

};
/******************************************************************************
*                    PhysicsBase	物理学基类
******************************************************************************/
class PhysicsBase
{
public:
	/*---------------- 常数 ----------------*/
	const double pi = 3.141592653589;				//圆周率
	const double G = 6.67408E-11;					//万有引力常数
	const double e = 1.602176565E-19;				//电子电荷
	const double epsilon0 = 8.85418781762E-12;		//真空电容率
	const double mu0 = 4 * pi * 1E-7;				//真空磁导率
};
#endif
