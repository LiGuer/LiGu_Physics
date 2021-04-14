#ifndef ELECTROMAGNETICS_H
#define ELECTROMAGNETICS_H
#include "../LiGu_AlgorithmLib/Mat.h"
#include "../LiGu_AlgorithmLib/Tensor.h"
/*----------------[ Electrostatic Field 静电场 ]----------------
*	[定义]: 
*	[静电场 麦克斯韦方程组]: ▽·E = ρ/ε0    ▽×E = 0
	[电势]: E = -▽φ
		=>	▽²φ = -ρ/ε0		Poisson's方程
		当ρ=0时, ▽²φ = 0	 Laplace's方程
		[三维直角坐标系]: ∂²φ/∂x² + ∂²φ/∂y² + ∂²φ/∂z² = -1/ε0·ρ(x,y,z)
	[算法]: 即 解Poisson's方程，格林函数，得 φ(r) = - 4π/ε0 ∫∫∫ f(rt) / |r-rt| d³rt
		//Tensor<double>* PoissonEquation(Mat<double>st, Mat<double>ed, Mat<double> delta, double (*f) (Mat<double>& x));
	[静电场唯一性定理]:
		对于各种边界条件，Poisson's方程可能有许多种解，但每个解的梯度相同. 
		静电场情况下, 意味着在边界条件下的满足泊松方程的势函数，所解得的电场确定且唯一.
**----------------------------------------------------------------------*/

/*----------------[ Magnetostatic Field 恒磁场 ]----------------
*	[定义]: 
*	[恒磁场 麦克斯韦方程组]: ▽·H = 0    ▽×H = 4π/c·J
	[磁矢势]: H = -▽×A		(▽·A = 0)
		=>	▽²A = -4π/c·J		Poisson's方程
	[算法]: 即 解Poisson's方程，格林函数, 得 A = 1/c ∫∫∫ J/r dV
		//Tensor<double>* PoissonEquation(Mat<double>st, Mat<double>ed, Mat<double> delta, double (*f) (Mat<double>& x));
**----------------------------------------------------------------------*/


/*----------------[ ElecmagnCell 电磁场基本结构体 ]----------------*/
struct ElecmagnCell
{
	double E[3] = { 0 }, H[3] = { 0 };			// 电场,磁场(三维)
	double sigmaE = 0, sigmaM = 0;				// 电导率,磁导率	// J =σE  Jm =σm H 
	double epsilon = 8.85E-12, mu = 1.2567E-6;	// ε介电常数,μ磁导系数 // D =εE  B =μH
};
/*----------------[ ComputationalElectromagnetics -  FDTD ]----------------
*	[方程]: 麦克斯韦方程组 (旋度那两个)
	[1] ▽ × H = ∂D/∂y + J
		∂Hz/∂y - ∂Hy/∂z = ε∂Ex/∂t + σEx
		∂Hx/∂z - ∂Hz/∂x = ε∂Ey/∂t + σEy
		∂Hy/∂x - ∂Hx/∂y = ε∂Ez/∂t + σEz
	[2] ▽ × E = - ∂B/∂y - Jm
		∂Ez/∂y - ∂Ey/∂z =  - μ∂Hx/∂t - σHx
		∂Ex/∂z - ∂Ez/∂x =  - μ∂Hy/∂t - σHy
		∂Ey/∂x - ∂Ex/∂y =  - μ∂Hz/∂t - σHz
**----------------------------------------------------------------------*/
void Electromagnetics(Tensor<ElecmagnCell>& Map, void (*setBoundaryValue) (Tensor<ElecmagnCell>& x, double time),
	double deltaTime, double deltaDim[], int EndTimes) {
	Tensor<ElecmagnCell> Mapprev(Map);
	for (int time = 0; time < EndTimes; time++) {
		for (int z = 1; z < Map.dimension[2] - 1; z++) {
			for (int y = 1; y < Map.dimension[1] - 1; y++) {
				for (int x = 1; x < Map.dimension[0] - 1; x++) {
					int r[] = { x,y,z };
					int rl[3][3] = { { x - 1,y,z },{ x,y - 1,z },{ x,y,z - 1 } };
					int rr[3][3] = { { x + 1,y,z },{ x,y + 1,z },{ x,y,z + 1 } };
					// CA,CB
					double CA = 1 - Map(r).sigmaE * deltaTime / 2 / Map(r).epsilon;
					double CB = deltaTime / Map(r).epsilon;
					CA /= 1 + Map(r).sigmaE * deltaTime / 2 / Map(r).epsilon;
					CB /= 1 + Map(r).sigmaE * deltaTime / 2 / Map(r).epsilon;
					// CP,CQ
					double CP = 1 - Map(r).sigmaM * deltaTime / 2 / Map(r).mu;
					double CQ =  - deltaTime / Map(r).mu;
					CP /= 1 + Map(r).sigmaM * deltaTime / 2 / Map(r).mu;
					CQ /= 1 + Map(r).sigmaM * deltaTime / 2 / Map(r).mu;
					for (int i = 0; i < 3; i++) {					//先算H 后算E, H要比E相差1/2
						// H(t+1)	// ∂Ez/∂y - ∂Ey/∂z		∂Ex/∂z - ∂Ez/∂x		∂Ey/∂x - ∂Ex/∂y
						Map(r).H[i] = CP * Map(r).H[i] + CQ * (
							(Mapprev(rr[(i + 1) % 3]).E[(i + 2) % 3] - Mapprev(r).E[(i + 2) % 3]) / deltaDim[(i + 1) % 3]
							- (Mapprev(rr[(i + 2) % 3]).E[(i + 1) % 3] - Mapprev(r).E[(i + 1) % 3]) / deltaDim[(i + 2) % 3]);
					}
					for (int i = 0; i < 3; i++) {
						// E(t+1)	// ∂Hz/∂y - ∂Hy/∂z		∂Hx/∂z - ∂Hz/∂x		∂Hy/∂x - ∂Hx/∂y
						Map(r).E[i] = CA * Map(r).E[i] + CB * (
							(Mapprev(r).H[(i + 2) % 3] - Mapprev(rl[(i + 1) % 3]).H[(i + 2) % 3]) / deltaDim[(i + 1) % 3]
							- (Mapprev(r).H[(i + 1) % 3] - Mapprev(rl[(i + 2) % 3]).H[(i + 1) % 3]) / deltaDim[(i + 2) % 3]);
					}
				}
			}
		}setBoundaryValue(Map, time* deltaTime);
		Mapprev = Map;
	}
}


#endif
