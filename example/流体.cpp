#include "../LiGu_Codes/LiGu_Physics/src/Electromagnetics.h"
#include "../LiGu_Codes/LiGu_Graphics/src/GraphicsND.h"
#include <algorithm>

#define N 500
int main() {
	GraphicsND G(N, N);
	double dt = 0.1; int tt = 0;
	Mat<> Water(N, N), P(N, N), velocX(N, N), velocY(N, N), velocZ(N, N), velocXTmp(N, N), velocYTmp(N, N), x(3), dx(3), g(3), v(3); dx.fill(1); 
	{//初始化Water
		for (int i = 0; i < N * N; i++)
			if (fabs(velocX.i2x(i) - N / 2) < 20 && fabs(velocX.i2y(i) - N / 2) < 20) velocX[i] = 5;
	}
	std::vector<double> dustX, dustY;
	for (int x = 0; x < N; x++)
		for (int y = 0; y < N; y++) 
			dustX.push_back(x),
			dustY.push_back(y);
	while (true) {
		for (int i = 0; i < velocX.size(); i++) {
			NavierStokesEquations(
				x.getData(velocX.i2x(i), velocX.i2y(i), 0), dx, dt,
				[&velocX](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
					return velocX(x, y);
				},
				[&velocY](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
					return velocY(x, y);
				},
				[](Mat<>& x) {return 0.0; },
				[&P](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
					return P(x, y);
				}, g, v.zero()
				); velocXTmp[i] = v[0], velocYTmp[i] = v[1];
		}
		velocX = velocXTmp;
		velocY = velocYTmp;
		Mat<> divV(N, N), Ptmp(N, N);
		for (int i = 0; i < divV.size(); i++)
			divV[i] = Calculus::Div(x.getData(velocX.i2x(i), velocX.i2y(i), 0), dx,
				[&velocX](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
					return velocX(x, y);
				},
				[&velocY](Mat<>& _x) {
					int x = _x[0]; x = x >= N ? N - 1 : x; x = x < 0 ? 0 : x;
					int y = _x[1]; y = y >= N ? N - 1 : y; y = y < 0 ? 0 : y;
					return velocY(x, y);
				},
				[](Mat<>& x) { return 0.0; }
			);
		for (int i = 0; i < P.size(); i++) {
			int x = velocX.i2x(i), 
				y = velocX.i2y(i);
			if (x >= N - 1 || x < 1 || y >= N - 1 || y < 1) continue;
			Ptmp[i] = (P(x + 1, y) + P(x - 1, y) + P(x, y + 1) + P(x, y - 1) - divV(x, y) / dt * dx[0] * dx[1]) / 4;
		} P = Ptmp;


		Water.zero();
		for (int i = 0; i < dustX.size(); i++) {
			int x = (int)dustX[i]; x = std::min(std::max(x, 0), N - 1);
			int y = (int)dustY[i]; y = std::min(std::max(y, 0), N - 1);
			dustX[i] += velocX(x, y) * dt;
			dustY[i] += velocY(x, y) * dt;
			x = (int)dustX[i]; x = std::min(std::max(x, 0), N - 1);
			y = (int)dustY[i]; y = std::min(std::max(y, 0), N - 1);
			Water(x, y) += 1;
		}
		
		if (tt++ % 20 == 0) {
			char filename[] = "D:/LiGuImg/liguXxxx.ppm"; filename[16] = tt / 100 + '0'; filename[17] = (tt % 100) / 10 + '0'; filename[18] = tt % 10 + '0';
			filename[15] = 'V'; G.contour(velocX);G.g.writeImg(filename); G.g.writeImg("D:/LiGu.ppm");
			filename[15] = 'X'; G.contour(velocX, velocY, velocZ);  G.g.writeImg(filename);
			filename[15] = 'P'; G.contour(P);	  G.g.writeImg(filename);
			filename[15] = 'W'; G.contour(Water); G.g.writeImg(filename);
		}
	}
}