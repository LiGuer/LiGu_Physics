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
    boundCondPoint boundCondPoint_list[MAXN];
    int boundCondPoint_N = 0;

    static double Map[MAXN][MAXN];
    static double veloc[MAXN][MAXN];
    static double accel[MAXN][MAXN];
    void init(){
        memset(Map, 0, sizeof(Map));
        memset(veloc, 0, sizeof(veloc));
    }
    void simulate() {
        //波动方程：该点与周围点平均值的差值，便是该点目前的加速度。
        //for(int i=0;i<N;i++){
        //    for(int j=0;j<N;j++){
        //        printf("%f ",Map[i][j]);
        //    }printf("\n");
        //}printf("==\n");
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

    void addBoundCondPoint(int x,int y,double amplitude,int type){
        boundCondPoint_list[boundCondPoint_N].x = x;
        boundCondPoint_list[boundCondPoint_N].y = y;
        boundCondPoint_list[boundCondPoint_N].amplitude = amplitude;
        boundCondPoint_N++;
    }

    void clearBoundCondPoint(){
        boundCondPoint_N = 0;
    }

    void boundCond(int time)
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
};
#endif
