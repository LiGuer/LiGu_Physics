#ifndef PHYSICSBASE_H
#define PHYSICSBASE_H
#include "Ligu_Graphics/Mat.h"
/******************************************************************************
*                    Particle	�ʵ�
******************************************************************************/
typedef Mat<double> Vector;
class Particle {
public:
	Vector r = Vector(3, 1), v = Vector(3, 1), a = Vector(3, 1);//λ��//�ٶ�
	double mass = 0, eleccharge = 0;				//����
	/*---------------- ���� ----------------*/
	void Momentum(Vector& ans) {
		ans.mult(mass, v);
	}
};
/******************************************************************************
*                    RigidBody	����
******************************************************************************/
class RigidBody:public Particle
{
public:
	Vector Angle = Vector(3, 1), AngularVeloc = Vector(3, 1), AngularAccele = Vector(3, 1);//�Ƕ�//���ٶ�//�Ǽ��ٶ�
	Mat<double> Inertia = Mat<double>(3, 3);		//ת������
	/*---------------- �Ƕ��� ----------------*/
	void AngularMomentum(Vector& ans) {
		ans.mult(Inertia, AngularVeloc);
	}

};
/******************************************************************************
*                    PhysicsBase	����ѧ����
******************************************************************************/
class PhysicsBase
{
public:
	/*---------------- ���� ----------------*/
	const double pi = 3.141592653589;				//Բ����
	const double G = 6.67408E-11;					//������������
	const double e = 1.602176565E-19;				//���ӵ��
	const double epsilon0 = 8.85418781762E-12;		//��յ�����
	const double mu0 = 4 * pi * 1E-7;				//��մŵ���
};
#endif
