# LiGu·Physics 物理学类
* <Dynamics.h>          动力学类
* <Electromagnetics>    电磁场类

* 注: 项目所需的<Mat.h><Tensor.h>基础, 见本作者<LiGu_AlgorithmLib>仓库内.  

## <Dynamics.h> 动力学类
```
double MassCentre(Mat<double>* r, Mat<double>* v, double* m, int N, Mat<double>& rc, Mat<double>& vc);  //质心
void CoordinateTrans(Mat<double>* r, Mat<double>* v, Mat<double>* a,int N, Mat<double>& coord);         //伽利略坐标变换
double CentrifugalPotential(double Omega, Mat<double>& center, Mat<double>& position);                  //离心势
void GravitationAcceleration(Mat<double>* r, Mat<double>* v, double* mass, Mat<double>* accele, int N); //万有引力
void Gravitation(Mat<double>* r, Mat<double>* v, double* mass, int N, double dt, int enpoch);           //引力运动
double* GravitatePotential(Mat<double>* r, double* mass, int N);                                        //引力势
```

## LiguのWave·波动方程·水波物理引擎·v1.0
* 任意水波仿真。
* 通过设置边界条件，可以模拟[干涉],[衍射]等物理教学情境。
* 软件功能：  
    [1]. 暂停  
    [2]. 初始波面  
    [3]. 设置边界条件  
        [3.1] 格式:[坐标X, 坐标Y, 幅度A]  
        [3.2] 当幅度设"0"，相当于墙，可以模拟衍射  
![img1.干涉演示](hhttps://github.com/LiGuer/Ligu_Physics/example/001.png)  
![img2.衍涉演示](hhttps://github.com/LiGuer/Ligu_Physics/example/img2.png)  