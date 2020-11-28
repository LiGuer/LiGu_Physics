#include "stdio.h"
#include "Universe.h"
#include "LiGu_GRAPHICS/Graphics.h"
#include "LiGu_GRAPHICS/Plot.h"
#include "PhysicsBase.h"
#include <iostream>

int main()
{
	Graphics* g = new Graphics;
	Plot p = Plot(g);
	p.setAxisRange(-10, -10, 10, 10);

	Particle a;
	{
		double r[3] = { -7,-7,0 };
		double v[3] = { 1,1,0 };
		double mass = 5;
		for (int i = 0; i < 3; i++) {
			a.r[i] = r[i];
			a.v[i] = v[i];
		}
		a.mass = mass;
	}
	Particle b;
	{
		double r[3] = { 7,7,0 };
		double v[3] = { -1,-1,0 };
		double mass = 5;
		for (int i = 0; i < 3; i++) {
			b.r[i] = r[i];
			b.v[i] = v[i];
		}
		b.mass = mass;
	}
	Kinematics kinem;
	Coordinate coord;
	{
		double r[] = { 1,1,0 };
		double v[] = { 0,0,0 };
		double Angle[] = { 0.1,0,0 };
		double AngularVeloc[] = { 0.1,0,0 };
		double mass = 5;
		for (int i = 0; i < 3; i++) {
			coord.r[i] = r[i];
			coord.v[i] = v[i];
			coord.Angle[i] = Angle[i];
			coord.AngularVeloc[i] = AngularVeloc[i];
		}
	}
	g->PaintSize = 2;
	g->PaintColor = 0xffffff;
	double dt = 0.01;
	for (int i = 0; i < 10000; i++) {
		Particle at;
		kinem.CoordinateTrans(a, coord, at);
		p.plotPoint(at.r[0],at.r[1]);
		Particle bt;
		kinem.CoordinateTrans(b, coord, bt);
		p.plotPoint(bt.r[0], bt.r[1]);
		coord.play(dt);
		kinem.play(a, dt);
		kinem.play(b, dt);
	}
	g->PaintSize = 5;
	g->PaintColor = 0xff0000;
	p.plotPoint(0, 0);
	g->PicWrite("D:/LIGU.ppm");
}
