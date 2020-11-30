#include<stdio.h>
#include "LiGu_Physics/Universe.h"
Universe universe;
int main() {
	int time = 0;
	while (true) {
		universe.Simulation();
		if (time++ % 10000 == 0)universe.showIt();
	}
}