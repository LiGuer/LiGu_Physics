#include "stdio.h"
#include "Universe.h"
#include "LiGu_GRAPHICS/Plot.h"
#include "PhysicsBase.h"
#include <iostream>

int main()
{
	Plot p = Plot();
	const double hx = 5E8;
	p.setAxisRange(-hx, -hx, hx, hx);
	Universe univer;
	Dynamics Dyna;
	Star* newStar;
	/*---------------- 地球 ----------------*/
	newStar = new Star;
	{
		double r[] = { 0,0 * hx,0 * hx };
		double v[] = { 0,0,0 };
		double m = 5.965E24;
		for (int i = 0; i < 3; i++) {
			newStar->r[i] = r[i];	newStar->v[i] = v[i];
		}
		newStar->mass = m;
	}univer.addStar(newStar);
	/*---------------- 月球 ----------------*/
	newStar = new Star;
	{
		double r[] = { 3.844039E8,0 * hx,0 * hx };
		double v[] = { 0,1.01775E3,0 };
		double m = 7.349E22;
		for (int i = 0; i < 3; i++) {
			newStar->r[i] = r[i];	newStar->v[i] = v[i];
		}
		newStar->mass = m;
	}univer.addStar(newStar);
	/*---------------- 人造卫星 ----------------*/
	newStar = new Star;
	{
		double r[] = { 3.5E8,0 * hx,0 * hx };
		double v[] = { 0,0.5E3,0 };
		double m = 11E3;
		for (int i = 0; i < 3; i++) {
			newStar->r[i] = r[i];	newStar->v[i] = v[i];
		}
		newStar->mass = m;
	}univer.addStar(newStar);
	/*---------------- 参考系 ----------------*/
	Particle* masscentre = Dyna.MassCentre(univer.Stars, 2);
	Coordinate coord;
	{
		double r[] = { masscentre->r[0],masscentre->r[1],masscentre->r[2] };
		double v[] = { masscentre->v[0],masscentre->v[1],masscentre->v[2] };
		double w = Dyna.TwobodyGravitationPeriod(*univer.Stars[0], *univer.Stars[1]);
		w += 0.5e-7;
		double AngularVeloc[] = { -w,0,0 };
		for (int i = 0; i < 3; i++) {
			coord.r[i] = r[i]; coord.v[i] = v[i];	coord.AngularVeloc[i] = AngularVeloc[i];
		}
	}
	/*---------------- 模拟 ----------------*/
	double dt = 10;  
	Particle tParti;
	for (int i = 0; i < 3600 * 24 * 365 / dt; i++) {
		univer.GravitationSimulation();
		for (int j = 0; j < univer.StarNum; j++) {
			Dyna.CoordinateTrans(*univer.Stars[j], coord, tParti);
			switch (j){
			case 0:p.g->PaintColor = 0x0000ff; break;
			case 1:p.g->PaintColor = 0xffffff; break;
			case 2:p.g->PaintColor = 0xffff00; break;
			}
			p.plotPoint(tParti.r[0], tParti.r[1]);
		}
		Dyna.play(univer.Stars, univer.StarNum, dt);
		coord.play(dt);
	}
	p.g->PicWrite("D:/LIGU.ppm");
}
