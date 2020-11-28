#ifndef UNIVERSE_H
#define UNIVERSE_H
#include <stdlib.h>
#include <stdio.h>
#include "Dynamics.h"
/******************************************************************************
*                    Star
******************************************************************************/
typedef Particle Star;
/******************************************************************************
*                    Universe
******************************************************************************/
class Universe
{
public:
	/*---------------- ���� ----------------*/
	static const int MaxStars = 2000;
	double dt = 1;										//ʱ��ֱ���
	int StarNum = 0;
	Star* Stars[MaxStars];
	const int StarDim = 3;
	Dynamics* dynamics = new Dynamics;
	/*---------------- ���� ----------------*/
	void addStar(Star* newstar){
	Stars[StarNum++] = newstar;
	}
	void GravitationSimulation() {
		for (int i = 0; i < StarNum; i++) {
			Stars[i]->a.clean();
			/*---------------- ���� ----------------*/
			double tmass = Stars[i]->mass; Stars[i]->mass = 0;
			dynamics->Gravitation(Stars, StarNum, Stars[i]->r, Stars[i]->a);
			Stars[i]->mass = tmass;
		}
	}
};
#endif