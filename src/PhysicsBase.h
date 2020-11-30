#ifndef PHYSICSBASE_H
#define PHYSICSBASE_H
#include "../LiGu_AlgorithmLib/Mat.h"
typedef Mat<double> Vector;
/******************************************************************************
*                    Particle	质点
******************************************************************************/
class Particle {
public:
	Vector r{ 3, 1 }, v{ 3, 1 }, a{ 3, 1 };					//位置 速度
	double mass = 0;								//质量
	void Momentum(Vector& ans) { ans.mult(mass, v); }	//动量 
};
/******************************************************************************
*                    RigidBody	刚体
******************************************************************************/
class RigidBody:public Particle{
public:
	Vector Angle{ 3, 1 }, AngularVeloc{ 3, 1 }, AngularAccele{ 3, 1 };	//角度//角速度//角加速度
	Mat<double> Inertia{ 3, 3 };				//转动惯量
	void AngularMomentum(Vector& ans) { ans.mult(Inertia, AngularVeloc); }//角动量

};
/******************************************************************************
*                    PhysicsBase	物理学基类
******************************************************************************/
class PhysicsBase
{
public:
	/*---------------- 物理常数 ----------------*/
	const double pi = 3.141592653589,					//圆周率
					G = 6.67408E-11,					//万有引力常数
					e = 1.602176565E-19,				//电子电荷
					epsilon0 = 8.85418781762E-12,		//真空电容率
					mu0 = 4 * pi * 1E-7;				//真空磁导率
};
#endif
