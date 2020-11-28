#ifndef DYNAMICS_H
#define DYNAMICS_H
#include <stdlib.h>
#include <math.h>
#include "PhysicsBase.h"
#include "Ligu_Graphics/Mat.h"
/******************************************************************************
*                    Coordinate	����ϵ
******************************************************************************/
class Coordinate
{
public:
	Vector r = Vector(3, 1); Vector v = Vector(3, 1); Vector a = Vector(3, 1);
	Vector Angle = Vector(3, 1), AngularVeloc = Vector(3, 1), AngularAccele = Vector(3, 1);//�Ƕ�//���ٶ�//�Ǽ��ٶ�
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
*                    Dynamics	��ѧ��
******************************************************************************/
class Dynamics : public PhysicsBase
{
public:
	/*---------------- ���� ----------------*/
	void play(Particle& point, double dt) {
		for (int dim = 0; dim < 3; dim++) {
			point.v[dim] += point.a[dim] * dt;//��һ��Ϊʲô���У�
			point.r[dim] += point.v[dim] * dt;
		}
	}
	void play(Particle* point[], int n, double dt) {
		for (int i = 0; i < n; i++) {
			for (int dim = 0; dim < 3; dim++) {
				point[i]->v[dim] += point[i]->a[dim] * dt;//��һ��Ϊʲô���У�
				point[i]->r[dim] += point[i]->v[dim] * dt;
			}
		}
	}
	/*---------------- ٤��������任 ----------------*/
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
	Particle* MassCentre(Particle* point[], int n);					//������
	/*---------------- �����������᳤�� ----------------
	A: �볤��		C: ����
	�Ƕ��� L = m r v��
	                     1     2    G M m
	����  E = Ek + Ep = --- m v  + -------
	                     2            r
	����	
	** ---------------------------------------- */
	void TwobodyGravitationAxis(Particle& PM, Particle& Pm, double& A, double& C) {
		double M = PM.mass, m = Pm.mass, r2 = 0, v2 = 0;						//r: PM Pm����//v: PM Pm����ٶ�
		for (int dim = 0; dim < 3; dim++) {
			r2 += (PM.r[dim] - Pm.r[dim]) * (PM.r[dim] - Pm.r[dim]);
			v2 += (PM.v[dim] - Pm.v[dim]) * (PM.v[dim] - Pm.v[dim]);
		}
		double r = sqrt(r2), vtheta = sqrt(v2);					//bug:vtheta
		double L = m * r * vtheta;								//L:�Ƕ���
		double E = 1.0 / 2 * m * v2 - G * M * m / r;			//E:����
		double p = (r2 * vtheta * vtheta) / (G * M);
		double epsi2 = 1 + (2 * E * p) / (G * M * m);
		A = p / (1 - epsi2);
		C = sqrt(epsi2) * A;
	}
	/*---------------- ����������ת���� ----------------
	����m��M������Ӱ�죬��
	Լ������: ��m = (M m) / (M +m)
	Fm = ��m am
	��
		 (M + m) m   ^     ->
	- G ------------ r = m am
			 r^2
				   3
	T = 2�� sqrt( A  / G (M + m) )
				 3    2                              2
	�����ն���: A  / T  = k ,  ���� k = (G M) / (4�� ) * (1 + m / M)
	** ---------------------------------------- */
	double TwobodyGravitationPeriod(Particle& M, Particle& m) {
		double A, C;
		TwobodyGravitationAxis(M, m, A, C);
		double Omega = sqrt(G * (M.mass + m.mass) / (A * A * A));		//��ת���ڵ�ƽ�����ٶ�
		return Omega;
	}
	/*---------------- �������� ----------------*/
	void Gravitation(Particle* point[], int n, Vector& position, Vector& acceleration);	//��������
	double GravitatePotential(Particle point[], int n, Vector& position);	//������
	double CentrifugalPotential(double Omega, Vector& center, Vector& position);	//������
};
#endif