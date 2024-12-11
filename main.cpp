//#include <iostream>
//#include <math.h>
//#include <cmath>  
//#include <vector>
//
//using namespace std;
//
//const double M_PI = 3.14159265358979323846;
//// 克拉索夫斯基椭球参数  
//const double a = 6378245.0; // 克拉索夫斯基椭球的长半轴  
//const double f = 1.0 / 298.3; // 扁率
//const double omega = 7.292115e-5; // 地球自转角速度
//const double fM = 3.98620e14; // 地心引力常数
//const double b = a * (1 - f);  // 短半轴  
//const double e2 = (2 - f) * f;  // 第一偏心率的平方  
//
//
//
//vector<double> crossProduct(const vector<double>& a, const vector<double>& b);
//double edcz(int n, double x, double a[], double b[]);
//double calculateDistanceToCenter(double latitudeDegrees);
//void toECEF(double lat, double lon, double Height, double& X_Dinxin, double& Y_Dinxin, double& Z_Dinxin);
//
//int main0(){
//	
//	double Time_all = 0.0; // 初始时间
//	double Time_fly = 0.0; // 飞行时间
//
//	// 发射点天文经纬度、瞄准方位角及发射点高度
//	double lambda_T = 35.0;	// 天文经度 35 度
//	double B_T = 120.0;		// 天文纬度 120 度
//	double A_T = 45.0;		// 瞄准方位角 45 度
//	double Height = 1000.0; // 发射高度1000米
//
//
//
//	// 经纬度转地心大地直角坐标系下位置坐标
//	double X_Dinxin, Y_Dinxin, Z_Dinxin;
//	double X_Dinxin0, Y_Dinxin0, Z_Dinxin0;	//地心大地直角坐标系下初始位置坐标
//	toECEF(lambda_T, B_T, Height, X_Dinxin0, Y_Dinxin0, Z_Dinxin0);
//	std::vector<double> C_Dinxin = { X_Dinxin0, Y_Dinxin0, Z_Dinxin0 };	// 保存位置
//	return 0;
//	//// 发射系下当前位置
//	//// double X_Launch = 0.0, Y_Launch = 0.0, Z_Launch = 0.0;
//	//vector<double> C_Launch = { 0.0, 0.0, 0.0 };
//
//	//// 发射系下速度m/s
//	//// double VX_Launch = 0.0, VY_Launch = 0.0, VZ_Launch = 0.0;
//	//vector<double> V_Launch = { 0.0, 0.0, 0.0 };
//
//	//// 初始重量KG
//	//double Weight = 44000 + 18800 + 5700 + 1000;
//
//	//// 所需变量
//	////int n = 8;	// y[]元素数量
//	////double y[] = { Time_fly, Height, C_Launch[0],C_Launch[1],C_Launch[2], V_Launch[0], V_Launch[1], V_Launch[2] };
//
//
//
//	//// -------------------开始主循环-------------------
//	//double dt = 0.02;		// 步长
//	//double P0;				// 发动机额定推力
//	//double M_per_s;			// 秒耗量
//	//double fai_now;			// 飞行程序角
//	//double Cx, Cy;			// 阻力系数、升力系数
//	//double X, Y;			// 阻力、升力
//	//double r, r0;			// 地心矢径、初始地心矢径
//	//double g;				// 引力加速度大小
//	//double omega_x, omega_y, omega_z;	// 自转角速度分量
//	//// omega_x = omega * cos();
//	//vector<double> Acce_Yinli(3);	// 引力加速度
//	//vector<double> Body_Direction = { 1,0,0 };	// 发射系下弹体方向矢量，用于表示推力方向
//	//vector<double> X_Direction = { 1,0,0 };		// 阻力方向
//	//vector<double> Y_Direction = { 0,1,0 };		// 升力方向
//
//	//r0 = sqrt(C_Dinxin[0] * C_Dinxin[0] + C_Dinxin[1] * C_Dinxin[01] + C_Dinxin[2] * C_Dinxin[2]);    // 地心矢径长
//	//
//	//// 未起飞阶段，一级发动机工作
//	//while (SegmentLable == 0) {
//	//	// 时间累加
//	//	Time_all += 0.02;
//	//	// 计算额定推力、秒耗量
//	//	P0 = edcz((sizeof(P0_1_Time)/sizeof(P0_1_Time[0])), Time_all, P0_1_Time, P0_1) * 1000;
//	//	M_per_s = edcz((sizeof(P0_1_Time) / sizeof(P0_1_Time[0])), Time_all, P0_1_Time, M_per_s_1);
//	//	// 更新rocket质量
//	//	Weight -= dt * M_per_s;
//	//	// 计算引力加速度
//	//	g = -fM / (r0 * r0);               // 引力加速度大小
//	//	for (size_t i = 0; i < C_Dinxin.size(); ) {
//	//		Acce_Yinli[i] = g * C_Dinxin[i] / r0;
//	//		i++;
//	//	}
//
//	//	if (P0 / Weight >= sqrt(Acce_Yinli[0] * Acce_Yinli[0] + Acce_Yinli[1] * Acce_Yinli[1] + Acce_Yinli[2] * Acce_Yinli[2])) {
//	//		SegmentLable = 1;
//	//		cout << Time_all << "s\t开始起飞\n" << endl;
//	//	}
//	//}
//
//	//// xp为大气压力，ru为大气密度, a为音速
//	//// double* xp;
//	//// double* ru;
//	//// double* a;
//	//double h0 = 1000;
//	//double p, ru, a, p0, ru0, a0;
//	//double S_hengjie, S_penkou;
//	//vector<double> cross_Keshi, cross_QianLian, Acce_QianLian;
//	//daqi_US(h0, &p0, &ru0, &a0);
//	//double v;
//	//double Ma_Now, P;
//	//double A_P;
//	//// 计算发动机推力加速度
//
//	//vector<double> Acce_X(3), Acce_Y(3);
//	//vector<double> Acce_P(3);
//	//vector<double> Acce_All(3);
//	//// 起飞后一级发动机工作阶段		
//	//S_hengjie = 3.2;
//	//S_penkou = 3.1;
//	//while (SegmentLable == 1) {
//	//	// 时间累加
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_all, P0_1_Time, M_per_s_1);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// 计算引力加速度
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // 地心矢径长
//	//	X_Dinxin = X_Dinxin0 + C_Launch[0];
//	//	Y_Dinxin = Y_Dinxin0 + C_Launch[1];
//	//	Z_Dinxin = Z_Dinxin0 + C_Launch[2];
//	//	xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Lat, Lon);
//	//	Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Lat);
//	//	// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // 引力加速度大小
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//	//	// 计算科氏加速度
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// 计算牵连加速度
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// 计算发动机推力加速度
//
//	//	// 根据段标选定横截面积、喷口面积
//
//	//	// 根据当前位置高度计算大气压强p、空气密度ru、音速a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	P0 = edcz(n1, Time_all, P0_1_Time, P0_1) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// 推力大小转推力加速度
//	//	A_P = P / Weight;
//	//	// 发射系下推力方向
//	//	// 
//	//	fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
//	//	if (Time_all > 68) fai_now = 35;
//	//	fai_now = fai_now * M_PI / 180.0;
//	//	Body_Direction = { cos(fai_now) ,0, sin(fai_now) };
//	//	for (size_t i = 0; i < Acce_P.size(); ++i) {
//	//		Acce_P[i] = A_P * Body_Direction[i];
//	//	}
//	//	/*Acce_P[0] = F_P * cos(fai_now);
//	//	Acce_P[1] = F_P * sin(fai_now);
//	//	Acce_P[2] = 0;*/
//
//
//	//	// 计算空气动力
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// 阻力、升力大小
//	//	Ma_Now = v / a;	// 马赫数大小
//	//	Cx = edcz(sizeof(Cx_1) / sizeof(Cx_1[0]), Ma_Now, Ma_1, Cx_1);	// 
//	//	Cy = Cy_1;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// 转加速度
//
//	//	for (size_t i = 0; i < X_Direction.size(); ++i) {
//	//		X_Direction[i] = -V_Launch[i] / (v * sqrt(3));
//	//	}
//	//	//// y_direction
//	//	Y_Direction[0] = -V_Launch[1] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = V_Launch[0] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = 0;
//	//	for (size_t i = 0; i < Acce_X.size(); ++i) {
//	//		Acce_X[i] = X * X_Direction[i];
//	//		Acce_Y[i] = Y * Y_Direction[i];
//	//		Acce_X[i] = 0;
//	//		Acce_Y[i] = 0;
//	//	}
//	//	// 合加速度
//	//	for (size_t i = 0; i < Acce_All.size(); ++i) {
//	//		Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i] + Acce_P[i] + Acce_X[i] + Acce_Y[i];
//	//		V_Launch[i] += Acce_All[i] * dt;
//	//		C_Launch[i] += V_Launch[i] * dt;
//	//	}
//
//
//	//	//double y[] = { Time_all,Height,C_Launch[0],C_Launch[1] ,C_Launch[2],V_Launch[0],V_Launch[1],V_Launch[2] };
//	//	//double dy[] = { dt,Height,V_Launch[0],V_Launch[1],V_Launch[2],Acce_All[0],Acce_All[1] ,Acce_All[2] };
//
//	//	// 速度大于0 累计飞行时间
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75) {
//	//		SegmentLable = 2;
//	//		Weight = 18800 + 5700 + 1000;
//	//		cout << "一级发动机工作结束，当前飞行高度：" << (int)Height << "m\n" << "累计时间：" << Time_all << "s" << endl;
//	//		cout << "当前发射系下坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 2) {
//	//	// 时间累加
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_fly, P0_2_Time, M_per_s_2);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// 计算引力加速度
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // 地心矢径长
//	//	X_Dinxin = X_Dinxin0 + C_Launch[0];
//	//	Y_Dinxin = Y_Dinxin0 + C_Launch[1];
//	//	Z_Dinxin = Z_Dinxin0 + C_Launch[2];
//	//	xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Lat, Lon);
//	//	Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Lat);
//	//	// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // 引力加速度大小
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// 计算科氏加速度
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// 计算牵连加速度
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// 计算发动机推力加速度
//
//	//	// 根据段标选定横截面积、喷口面积
//	//	S_hengjie = 2.8;
//	//	S_penkou = 2.5;
//	//	// 根据当前位置高度计算大气压强p、空气密度ru、音速a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	// 
//	//	// 
//	//	P0 = edcz(n1, Time_all - 65.75, P0_2_Time, P0_2) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// 推力大小转推力加速度
//	//	A_P = P / Weight;
//	//	// 发射系下推力方向
//	//	// 
//	//	fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
//	//	if (Time_all > 68) fai_now = 35;
//	//	fai_now = fai_now * 3.1415926 / 180.0;
//	//	Body_Direction = { cos(fai_now) ,sin(fai_now) ,0 };
//	//	for (size_t i = 0; i < Acce_P.size(); ++i) {
//	//		Acce_P[i] = A_P * Body_Direction[i];
//	//	}
//
//	//	// 计算空气动力
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// 阻力、升力大小
//	//	Ma_Now = v / a;	// 马赫数大小
//	//	Cx = edcz(sizeof(Cx_2) / sizeof(Cx_2[0]), Ma_Now, Ma_2, Cx_2);	// 
//	//	Cy = Cy_2;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// 转加速度
//
//	//	for (size_t i = 0; i < X_Direction.size(); ++i) {
//	//		X_Direction[i] = -V_Launch[i] / (v * sqrt(3));
//	//	}
//	//	//// y_direction
//	//	Y_Direction[0] = -V_Launch[1] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = V_Launch[0] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = 0;
//	//	for (size_t i = 0; i < Acce_X.size(); ++i) {
//	//		Acce_X[i] = X * X_Direction[i];
//	//		Acce_Y[i] = Y * Y_Direction[i];
//	//		Acce_X[i] = 0;
//	//		Acce_Y[i] = 0;
//	//	}
//
//	//	// 合加速度
//	//	for (size_t i = 0; i < Acce_All.size(); ++i) {
//	//		Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i] + Acce_P[i] + Acce_X[i] + Acce_Y[i];
//	//		V_Launch[i] += Acce_All[i] * dt;
//	//		C_Launch[i] += V_Launch[i] * dt;
//	//	}
//
//
//	//	//double y[] = { Time_all,Height,C_Launch[0],C_Launch[1] ,C_Launch[2],V_Launch[0],V_Launch[1],V_Launch[2] };
//	//	double dy[] = { dt,Height,V_Launch[0],V_Launch[1],V_Launch[2],Acce_All[0],Acce_All[1] ,Acce_All[2] };
//
//	//	// 速度大于0 累计飞行时间
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75 + 53.15) {
//	//		SegmentLable = 3;
//	//		Weight = 5700 + 1000;
//	//		cout << "二级发动机工作结束，当前飞行高度：" << (int)Height << "m\n" << "累计时间：" << Time_all << "s" << endl;
//	//		cout << "当前发射系下坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 3) {
//	//	// 时间累加
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_fly, P0_3_Time, M_per_s_3);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// 计算引力加速度
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // 地心矢径长
//	//	Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // 引力加速度大小
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// 计算科氏加速度
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// 计算牵连加速度
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// 计算发动机推力加速度
//
//	//	// 根据段标选定横截面积、喷口面积
//	//	S_hengjie = 3.2;
//	//	S_penkou = 3.1;
//	//	// 根据当前位置高度计算大气压强p、空气密度ru、音速a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	// 
//	//	// 
//	//	P0 = edcz(n1, Time_all - 65.75 - 53.15, P0_3_Time, P0_3) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// 推力大小转推力加速度
//	//	A_P = P / Weight;
//	//	// 发射系下推力方向
//	//	// 
//	//	fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
//	//	if (Time_all > 68) fai_now = 35;
//	//	fai_now = fai_now * 3.1415926 / 180.0;
//	//	Body_Direction = { cos(fai_now) ,sin(fai_now) ,0 };
//	//	for (size_t i = 0; i < Acce_P.size(); ++i) {
//
//	//		Acce_P[i] = A_P * Body_Direction[i];
//	//	}
//	//	/*Acce_P[0] = F_P * cos(fai_now);
//	//	Acce_P[1] = F_P * sin(fai_now);
//	//	Acce_P[2] = 0;*/
//
//
//	//	// 计算空气动力
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// 阻力、升力大小
//	//	Ma_Now = v / a;	// 马赫数大小
//	//	Cx = edcz(sizeof(Cx_3) / sizeof(Cx_3[0]), Ma_Now, Ma_3, Cx_3);	// 
//	//	Cy = Cy_3;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// 转加速度
//
//	//	for (size_t i = 0; i < X_Direction.size(); ++i) {
//	//		X_Direction[i] = -V_Launch[i] / (v * sqrt(3));
//	//	}
//	//	//// y_direction
//	//	Y_Direction[0] = -V_Launch[1] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = V_Launch[0] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = 0;
//	//	for (size_t i = 0; i < Acce_X.size(); ++i) {
//	//		Acce_X[i] = X * X_Direction[i];
//	//		Acce_Y[i] = Y * Y_Direction[i];
//	//		Acce_X[i] = 0;
//	//		Acce_Y[i] = 0;
//	//	}
//	//	// 合加速度
//	//	for (size_t i = 0; i < Acce_All.size(); ++i) {
//	//		Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i] + Acce_P[i] + Acce_X[i] + Acce_Y[i];
//	//		V_Launch[i] += Acce_All[i] * dt;
//	//		C_Launch[i] += V_Launch[i] * dt;
//	//	}
//
//
//	//	double y[] = { Time_all,Height,C_Launch[0],C_Launch[1] ,C_Launch[2],V_Launch[0],V_Launch[1],V_Launch[2] };
//	//	double dy[] = { dt,Height,V_Launch[0],V_Launch[1],V_Launch[2],Acce_All[0],Acce_All[1] ,Acce_All[2] };
//
//	//	// 速度大于0 累计飞行时间
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75 + 53.15 + 47.5) {
//	//		SegmentLable = 4;
//	//		Weight = 1000;
//	//		cout << "三级发动机工作结束，当前飞行高度：" << (int)Height << "m\n" << "累计时间：" << Time_all << "s" << endl;
//	//		cout << "当前发射系下坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 4) {
//	//	// 时间累加
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//
//	//	// 计算引力加速度
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // 地心矢径长
//	//	Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // 引力加速度大小
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// 计算科氏加速度
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// 计算牵连加速度
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//
//	//	// 合加速度
//	//	for (size_t i = 0; i < Acce_All.size(); ++i) {
//	//		Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i];
//	//		V_Launch[i] += Acce_All[i] * dt;
//	//		C_Launch[i] += V_Launch[i] * dt;
//	//	}
//
//
//	//	//double y[] = { Time_all,Height,C_Launch[0],C_Launch[1] ,C_Launch[2],V_Launch[0],V_Launch[1],V_Launch[2] };
//	//	//double dy[] = { dt,Height,V_Launch[0],V_Launch[1],V_Launch[2],Acce_All[0],Acce_All[1] ,Acce_All[2] };
//
//	//	// 速度大于0 累计飞行时间
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Height <= 80000) {
//	//		SegmentLable = 5;
//	//		cout << "即将再入大气层飞行，当前飞行高度：" << (int)Height << "m\n" << "累计时间：" << Time_all << "s" << endl;
//	//		cout << "当前发射系下坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//
//	////// 计算微分（速度、加速度）
//	////double dy[] = { 0,0,0,0,0,0,0,0 };
//	////double yc[] = {0,0,0,0,0,0,0,0};
//	////double y1[] = { 0,0,0,0,0,0,0,0 };
//
//	////// 隆格库塔积分
//	////runk(n, Height, y, dy, yc, y1);
//
//	////// 更新发射系下位置、速度
//
//
//	////// 更新当前飞行高度
//
//	//cout << "发射系下落点坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//
//
//	system("pause");
//
//}
//
//
//
