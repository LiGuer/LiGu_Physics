# LiGu·Physics 物理学类
## Codes List
```
<Calculus.h>				微积分 / 微分方程
<Dynamics.h>			动力学类
<Electromagnetics.h>		电动力学类
<QuantumMechanics.h>	量子力学类
```

## <Calculus.h> 微积分 / 微分方程
```
double diff		(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);		\\导数 (N阶)
double diff_	(double x0, double(*f)(double x), int N = 1, double dx = 1E-3);
double Curvature(double x0, double(*y)(double x),			 double dx = 1E-3);		\\曲率
double PartiDeriv	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数
double PartiDeriv2	(Mat<>& x, int index, double dx, double(*func)(Mat<>& x));		\\偏导数 (2阶)
Mat<>& TaylorFormula(double x0, double(*f)(double x), Mat<>& Coeff, int N = 3);		\\Taylor展开
double Exp		(double x, int N = 18);												\\常用函数
double Sin		(double x, int N = 18);
double Cos		(double x, int N = 18);
double lnOneAdd	(double x, int N = 18);
double Arctan	(double x, int N = 18);
double PowOneAdd(double x, double p, int N = 18);
void   RungeKutta(Mat<>& y, double dt, double t0, int enpoch, Mat<>& (*derivY)(double t, Mat<>& y));
																					\\解微分方程: RungeKutta法
	   PoissonEquation();															\\Poisson's方程
	   WaveEquation();																\\波动方程
```

## <Dynamics.h> 动力学类
```
void run(Mat<>& x, Mat<>&dx, double dt, int enpoch, void(*AcceleFun)(Mat<>& x, Mat<>& dx, Mat<>& ddx));		//运动 (Runge Kutta)
Mat<>& Lagrange(Mat<>& x, Mat<>& dx, Mat<>& ddx, double(*T)(Mat<>& dx), double(*U)(Mat<>& x));				//Lagrange 方程
double MassCentre(Mat<>* r, Mat<>* v, double* m, int N, Mat<>& rc, Mat<>& vc);									//质心
double CentrifugalPotential(double angularVeloc, Mat<>& center, Mat<>& position);									//CentrifugalPotential
void GravitationAcceleration(Mat<>* r, double* mass, Mat<>* accele, int N);											//万有引力
double* GravitatePotential(Mat<>* r, double* mass, int N);															//引力势
```
* [Example]
```
[Example]: //N连杆摆动
	Mat<> angle(10), angleVeloc(10); angle.fill(PI / 2);
	Dynamics::run(
		angle, angleVeloc, 0.01, 1,
		[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
			Dynamics::Lagrange(
				x, dx, ddx,
				[](Mat<>& dx) {return 1.0 / 8 * (dx.dot(dx)); },
				[](Mat<>&  x) {
					double sum = 0;
					for (int i = 0; i < x.size(); i++) sum += (x.size() - i) * cos(x[i]);
					return 10.0 / 2 * (sum);
				}
			);
		}
	);
[Example] //重力 + 抛体 + 弹簧
	Mat<> x(2), v(2); v.fill(40);
	Dynamics::run(
		x, v, 0.1, 1,
		[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
			Dynamics::Lagrange(
				x, dx, ddx,
				[](Mat<>& dx) { return 1.0 / 2 * dx.dot(dx); },
				[](Mat<>& x)  { return x[1] * 10 
				+ 0.05 * (pow(x[0], 2) + pow(x[1], 2))
				; }
			);
		}
	);
[Example] 弹性碰撞
	#define N 20
	Mat<> x(2 * N), v(2 * N);
	x.rands(2 * N, 1, -200, 400);
	v.rands(2 * N, 1, -50, 50);
	Dynamics::run(
		x, v, 0.1, 1,
		[](Mat<>& x, Mat<>& dx, Mat<>& ddx) {
			Dynamics::Lagrange(
				x, dx, ddx,
				[](Mat<>& dx) { return 1.0 / 2 * dx.dot(dx); },
				[](Mat<>& x) {
					double sum = 0;
					for (int i = 0; i < N; i++) {
						for (int j = i + 1; j < N; j++) {
							double d = sqrt(pow(x[2 * i] - x[2 * j], 2) + pow(x[2 * i + 1] - x[2 * j + 1], 2));
							sum += d > 40 ? 0 : 1000 * (40 - d);
						}
					}return sum;
				}
			);
		}
	);
```
## <Electromagnetics.h> 电动力学类
```
void Electromagnetics(Tensor<ElecmagnCell>& Map, void (*setBoundaryValue) (Tensor<ElecmagnCell>& x, double time),
	double deltaTime, double deltaDim[], int EndTimes)							//ComputationalElectromagnetics -  FDTD
```
## <QuantumMechanics.h> 量子力学类
```
void Schrodinger_1D(double m, double dx, double xs, double xe, double(*U)(double x), Mat<double>& E, Mat<double>& Psi) ;//Schrödinger 方程
```