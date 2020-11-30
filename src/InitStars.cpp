#include "LiGu_Physics/Universe.h"
void Universe::InitStars() {
	/*---------------- Ì«Ñô ----------------*/
	{
		double r[] = { 0,0,0 };
		double v[] = { 0,0,0 };
		Stars[StarNum].mass = 1.9891E30;
		Stars[StarNum].color = 0xFF0000;
		for (int i = 0; i < 3; i++) {
			Stars[StarNum].r[i] = r[i];	Stars[StarNum].v[i] = v[i];
		}
	}StarNum++;
	/*---------------- µØÇò ----------------*/
	{
		double r[] = { 1.496E11,0,0 };
		double v[] = { 0,29.783E3,0 };
		Stars[StarNum].mass = 5.965E24;
		Stars[StarNum].R = 6378.137E3;
		Stars[StarNum].color = 0x0000FF;
		for (int i = 0; i < 3; i++) {
			Stars[StarNum].r[i] = r[i];	Stars[StarNum].v[i] = v[i];
		}
	}StarNum++;
	/*---------------- ÔÂÇò ----------------*/
	{
		double r[] = { 1.496E11 + 3.844039E8,0,0 };
		double v[] = { 0,29.783E3 + 1.01775E3,0 };
		Stars[StarNum].mass = 7.349E22;
		Stars[StarNum].R = 1738E3;
		Stars[StarNum].color = 0xFFFFFF;
		for (int i = 0; i < 3; i++) {
			Stars[StarNum].r[i] = r[i];	Stars[StarNum].v[i] = v[i];
		}
	}StarNum++;
	/*---------------- Ä¾ÐÇ ----------------*/
	{
		double r[] = { -8165208E5,0,0 };
		double v[] = { 0,-13069.7,0 };
		Stars[StarNum].mass = 1.90E27;
		Stars[StarNum].R = 71492E3;
		Stars[StarNum].color = 0xAA0000;
		for (int i = 0; i < 3; i++) {
			Stars[StarNum].r[i] = r[i];	Stars[StarNum].v[i] = v[i];
		}
	}StarNum++;
	/*---------------- ÈËÔìÎÀÐÇ ----------------*/
	{
		double r[] = { 3.5E10,0,0 };
		double v[] = { 0,5E4,0 };
		Stars[StarNum].mass = 1E8;
		Stars[StarNum].R = 10;
		Stars[StarNum].color = 0xFFFF00;
		Stars[StarNum].mass = 11E3;
		for (int i = 0; i < 3; i++) {
			Stars[StarNum].r[i] = r[i];	Stars[StarNum].v[i] = v[i];
		}
	}StarNum++;
}