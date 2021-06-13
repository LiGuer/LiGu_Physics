#include "../LiGu_Codes/LiGu_Physics/src/Electromagnetics.h"
#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include <algorithm>

#define N 500
#define M 150
int main() {
	GraphicsND G(M, N);
	double dt = 0.1; int tt = 0, ttt = 0;
	Mat<> Water(N, M), P(N, M), velocX(N, M), velocY(N, M), velocZ(N, M), velocXTmp(N, M), velocYTmp(N, M), x(3), dx(3), g(3), v(3); dx.fill(1); 
	std::vector<double> dustX, dustY;
	{//初始化Water
		for (int i = 0; i < N * M; i++)
			if (fabs(velocX.i2x(i) - M / 2) < 20 && fabs(velocX.i2y(i) - M / 2) < 20) velocX[i] = 5;
	}
	for (int i = 0; i < 100000; i++) {
		dustX.push_back(rand() / (double)RAND_MAX * 40 + M / 2 - 20),
		dustY.push_back(rand() / (double)RAND_MAX * 40 + M / 2 - 20);
	}
	while (true) {
		for (int i = 0; i < velocX.size(); i++) {
			NavierStokesEquations(
				x.getData(velocX.i2x(i), velocX.i2y(i), 0), dx, dt,
				[&velocX](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= M ? M - 1 : y; y = y < 0 ? 0 : y;
					return velocX(x, y);
				},
				[&velocY](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= M ? M - 1 : y; y = y < 0 ? 0 : y;
					return velocY(x, y);
				},
				[](Mat<>& x) {return 0.0; },
				[&P](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= M ? M - 1 : y; y = y < 0 ? 0 : y;
					return P(x, y);
				}, g, v.zero()
				); velocXTmp[i] = v[0], velocYTmp[i] = v[1];
		}
		velocX = velocXTmp;
		velocY = velocYTmp;
		Mat<> divV(N, M), Ptmp(N, M);
		for (int i = 0; i < divV.size(); i++)
			divV[i] = Calculus::Div(x.getData(velocX.i2x(i), velocX.i2y(i), 0), dx,
				[&velocX](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= M ? M - 1 : y; y = y < 0 ? 0 : y;
					return velocX(x, y);
				},
				[&velocY](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= M ? M - 1 : y; y = y < 0 ? 0 : y;
					return velocY(x, y);
				},
				[](Mat<>& x) { return 0.0; }
			);
		for (int i = 0; i < P.size(); i++) {
			int x = velocX.i2x(i), 
				y = velocX.i2y(i);
			if (x >= N - 1 || x < 1 || y >= M - 1 || y < 1) continue;
			Ptmp[i] = (P(x + 1, y) + P(x - 1, y) + P(x, y + 1) + P(x, y - 1) - divV(x, y) / dt * dx[0] * dx[1]) / 4;
		} P = Ptmp;


		Water.zero();
		for (int i = 0; i < dustX.size(); i++) {
			int x = dustX[i]; x = std::min(std::max(x, 0), N - 1);
			int y = dustY[i]; y = std::min(std::max(y, 0), M - 1);
			dustX[i] += velocX(x, y) * dt;
			dustY[i] += velocY(x, y) * dt;
			x = dustX[i]; x = std::min(std::max(x, 0), N - 1);
			y = dustY[i]; y = std::min(std::max(y, 0), M - 1);
			if (Water(x, y) < 50) Water(x, y) += 1;
		}
		
		if (tt++ % 20 == 0) { ttt += 1;
			char filename[] = "D:/LiGuImg/liguWxxx.ppm"; filename[16] = ttt / 100 + '0'; filename[17] = (ttt % 100) / 10 + '0'; filename[18] = ttt % 10 + '0';
			G.contour(Water, Water, Water); G.g.writeImg("D:/LiGu.ppm"); G.g.writeImg(filename);
		}
	}
}