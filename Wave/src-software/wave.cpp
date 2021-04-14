#include "wave.h"
double wave::Map[MAXN][MAXN];
double wave::veloc[MAXN][MAXN];
double wave::accel[MAXN][MAXN];
boundCondPoint wave::boundCondPoint_list[MAXN*MAXN];

void wave::init(){
    memset(Map, 0, sizeof(Map));
    memset(veloc, 0, sizeof(veloc));
}

void wave::simulate() {
    //波动方程：该点与周围点平均值的差值，便是该点目前的加速度。
    memset(accel, 0, sizeof(accel));
    int x_step[] = { 0, 0, 1, -1, 1, 1, -1, -1 };
    int y_step[] = { 1, -1, 0, 0, 1, -1, -1, 1 };
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < 8; k++) {
                int xt = i + x_step[k];
                int yt = j + y_step[k];
                if (xt > 0 && xt < N && yt > 0 && yt < N) {
                    accel[i][j] += Map[xt][yt];
                }
            }
            accel[i][j] = accel[i][j] / 8 - Map[i][j];
        }
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            veloc[i][j] = veloc[i][j] + accel[i][j] * dt;
            Map[i][j] = Map[i][j] + veloc[i][j] * dt;
        }
    }
}

void wave::addBoundCondPoint(int x,int y,double amplitude,int type){
    boundCondPoint_list[boundCondPoint_N].x = x;
    boundCondPoint_list[boundCondPoint_N].y = y;
    boundCondPoint_list[boundCondPoint_N].amplitude = amplitude;
    boundCondPoint_N++;
}

void wave::clearBoundCondPoint(){
    boundCondPoint_N = 0;
}

void wave::boundCond(int time)
{
    const double PI = 3.141592653589;
    for(int i=0; i<boundCondPoint_N; i++){
        int x = boundCondPoint_list[i].x;
        int y = boundCondPoint_list[i].y;
        double amplitude = boundCondPoint_list[i].amplitude;
        Map[x][y] = amplitude * sin(time / (2 * PI));
        veloc[x][y] = 0;
    }
}
