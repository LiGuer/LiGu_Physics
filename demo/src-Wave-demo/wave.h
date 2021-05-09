#ifndef WAVE_H
#define WAVE_H
#include<string.h>
#include<stdio.h>
#include<math.h>
struct boundCondPoint{
    int x,y;
    double amplitude;
    enum type{height,veloc};
};
class wave
{
private:
public:
    static const int MAXN = 500;
    double dt = 1;
    int N = 500;
    int boundCondPoint_N = 0;
    static boundCondPoint boundCondPoint_list[MAXN*MAXN];

    static double Map[MAXN][MAXN];
    static double veloc[MAXN][MAXN];
    static double accel[MAXN][MAXN];

    void init();
    void simulate();
    void addBoundCondPoint(int x,int y,double amplitude,int type);
    void clearBoundCondPoint();
    void boundCond(int time);
};
#endif
