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
//const double f = 1.0 / 298.3; // 扁率( f = (a-b) / a )
//const double b = a * (1 - f);  // 短半轴  
//const double e2 = (2 - f) * f;  // 第一偏心率的平方  ( e = pow(a*a-b*b , 0.5) / b )
//const double omega = 7.292115e-5; // 地球自转角速度
//const double fM = 3.98620e14; // 地心引力常数
//const double miu = 26.33281e24; // 椭球体变率常数(m^5 / s^2)
//
//// rocket结构体
//struct rocket {
//	double work_time = 0; // 工作时间，总时间
//	double fly_time = 0; // 飞行时间
//	double Height = 0; // 飞行高度
//	double Height0 = 0; // 飞行高度
//	double Weight = 0; // 重量
//	vector<double> xyz_fs{ 0, 0 , 0 }; // 发射系下坐标
//	vector<double> v_xyz_fs{ 0, 0 , 0}; // 发射系下速度
//	vector<vector<double>> xyz_fs_storage = {xyz_fs}; // 发射系下坐标存储
//	vector<double> xyz_dixin{ 0,0,0 };// 地心直角坐标系下rocket的坐标
//	vector<double> xyz_dixin0{ 0,0,0 };// 地心直角坐标系下发射点的坐标
//	vector<double> r_fai_lambda = { 0,0,0 }; // 地心系经纬度与到地心距离
//};
//
//// 额定推力参数P_parameter结构体
//struct P_parameter {
//	vector<double> Time; // 时间
//	vector<double> P0; // 额定推力大小 
//	vector<double> M_per_s; // 秒耗量
//	int n = 0; // 参数数量
//	double S_hengjie; // 横截面积
//	double S_penkou; // 喷口截面积
//};
//
//// 调用函数
//double crossProduct(double a[], double b[], double result[]);
//double edcz(int n, double x, vector<double> a, vector<double> b);
//double calculateDistanceToCenter(double latitudeDegrees);
//void toECEF(double lat, double lon, double Height, double& X_Dinxin, double& Y_Dinxin, double& Z_Dinxin);
//void runk(int n, double h, double y[], double dy[], double yc[], double y1[]);
//double ECEFtoWGS84(double& Longitude, double& Latitude, double& Height, double x, double y, double z);
////大气    xp为大气压力，ru为大气密度,a为音速
//void daqi_US(double h_daqi, double* xp, double* ru, double* a);
//
//// 加速度的大小计算
//double Acce_count(vector<double> Acce) {
//	double count = 0;
//	for (int i = 0; i < 3; i++) {
//		count += Acce[i] * Acce[i];
//	}
//	return pow(count, 0.5);
//}
//
//// 发射系与地心直角坐标系的转换
//void g_s(vector<double>& g, vector<double>& s, int flag) {
//	// g：发射系坐标  s：地心系坐标  flag：1--发射到地心，0--地心到发射，3--omega
//
//	double x_g = g[0], y_g = g[1], z_g = g[2]; // 发射系坐标
//	double x_s = s[0], y_s = s[1], z_s = s[2]; // 地心系坐标
//
//	// 发射点在地心系下坐标
//	double lambda_T = 120.0;	// 天文经度 120 度
//	double B_T = 35.0;		// 天文纬度 35 度
//	double A_T = 45.0;		// 瞄准方位角 45 度
//	double Height = 1000.0; // 发射高度1000米
//	double x0, y0, z0;	//地心大地直角坐标系下初始位置坐标
//	// 角度转弧度
//	lambda_T = lambda_T * M_PI / 180.0;
//	B_T = B_T * M_PI / 180.0;
//	A_T = A_T * M_PI / 180.0;
//	toECEF(lambda_T, B_T, Height, x0, y0, z0);
//
//	// 转换矩阵
//	double d11 = -sin(lambda_T) * sin(A_T) - sin(B_T) * cos(A_T) * cos(lambda_T);
//	double d12 = cos(lambda_T) * cos(B_T);
//	double d13 = -sin(lambda_T) * cos(A_T) + cos(lambda_T) * sin(B_T) * sin(A_T);
//	double d21 = sin(A_T) * cos(lambda_T) - sin(lambda_T) * sin(B_T) * cos(A_T);
//	double d22 = sin(lambda_T) * cos(B_T);
//	double d23 = cos(lambda_T) * cos(A_T) + sin(lambda_T) * sin(B_T) * sin(A_T);
//	double d31 = cos(B_T) * cos(A_T);
//	double d32 = sin(B_T);
//	double d33 = -sin(A_T) * cos(B_T);
//
//	if (flag == 1) { // 发射到地心
//		x_s = d11 * x_g + d12 * y_g + d13 * z_g + x0;
//		y_s = d21 * x_g + d22 * y_g + d23 * z_g + y0;
//		z_s = d31 * x_g + d32 * y_g + d33 * z_g + z0;
//
//		s = { x_s ,y_s ,z_s };
//	}
//	else if(flag == 0) {// 地心到发射
//		x_g = d11 * (x_s - x0) + d21 * (y_s - x0) + d31 * (z_s - x0);
//		y_g = d12 * (x_s - x0) + d22 * (y_s - x0) + d32 * (z_s - x0);
//		z_g = d13 * (x_s - x0) + d23 * (y_s - x0) + d33 * (z_s - x0);
//		g = { x_g ,y_g ,z_g };
//	}
//	else if(flag == 3) {
//		x_g = d11 * (x_s) + d21 * (y_s) + d31 * (z_s);
//		y_g = d12 * (x_s) + d22 * (y_s) + d32 * (z_s);
//		z_g = d13 * (x_s) + d23 * (y_s) + d33 * (z_s);
//		g = { x_g ,y_g ,z_g };
//	}
//}
//
//// 引力加速度
//vector<double> Acce_g(rocket& temp) {
//	// 计算地心经纬度
//	double r = 0, fai = 0, lambda = 0;
//	// 计算rocket到地心距离
//	for (size_t i = 0; i < temp.xyz_dixin.size(); i++) {
//		r += temp.xyz_dixin[i] * temp.xyz_dixin[i];
//	}
//	r = pow(r, 0.5);
//	// 计算地心纬度
//	fai = acos(temp.xyz_dixin[2] / r);
//	// 计算地心经度
//	lambda = temp.xyz_dixin[1] / abs(temp.xyz_dixin[1]) * (M_PI / 2 - asin(temp.xyz_dixin[0] / pow(temp.xyz_dixin[0] * temp.xyz_dixin[0] + temp.xyz_dixin[1] * temp.xyz_dixin[1], 0.5)));
//
//	// 计算引力
//	double g_r = -fM / pow(r, 2) + miu / pow(r, 4) * (5 * pow(sin(fai), 2) - 1);
//	double g_omega = -2 * miu / pow(r, 4) * sin(fai);
//	vector<double> R_0_xyz = { 0,0,0 };
//	for (size_t i = 0; i < R_0_xyz.size(); i++) {
//		R_0_xyz[i] = temp.xyz_dixin0[i] + temp.xyz_fs[i];
//	}
//
//	// 地球自转角速度在地心直角坐标系和发射系下的分量
//	vector<double> omega_s = { 0,0,omega };
//	vector<double> omega_g = { 0, 0, 0 };
//	g_s(omega_g, omega_s, 3);
//
//	vector<double> g = { 0,0,0 };
//	for (size_t i = 0; i < R_0_xyz.size(); i++) {
//		g[i] = g_r / r * R_0_xyz[i] + g_omega / omega * omega_g[i];
//	}
//
//	double N = a / sqrt(1 - e2 * sin(fai) * sin(fai));
//	// test.Height = Acce_count(test.xyz_dixin) - N;
//
//	// 更新rocket的Height
//	temp.Height = Acce_count(temp.xyz_dixin) - N;
//	return g;
//}
//
//// 飞行程序角
//// 时间
//vector<double> fai_time = { 0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 64.0, 66.5, 68.0 }; // 后续等程序角飞行
//// 飞行程序角
//vector<double> fai = { 90.0, 90.0, 69.0, 59.0, 50.0, 40.0, 36.0, 36.0, 35.0 };
//
//// 推力加速度
//P_parameter P1, P2, P3;
//vector<double> Acce_p(rocket& temp, P_parameter P_segment) {
//
//	double h = 0.02; // 步长
//
//	// 插值计算额定推力
//	double P0 = edcz(P_segment.n, temp.work_time, P_segment.Time, P_segment.P0) * 1000;
//
//	// 根据当前位置高度计算大气压强p、空气密度ru、音速a
//	double p, ru, a, p0;
//	daqi_US(temp.Height0, &p0, &ru, &a);
//	daqi_US(temp.Height, &p, &ru, &a);
//	double P = P0 + P_segment.S_penkou * (1 - p / p0);
//	// 推力大小转推力加速度
//	double A_P = P / temp.Weight;
//
//	// 秒耗量
//	double M_per_s = edcz(P_segment.n, temp.work_time, P_segment.Time, P_segment.M_per_s);
//	temp.Weight -= M_per_s * h;
//
//	// 发射系下推力方向
//	// 飞行姿态角
//	double fai_now = edcz(fai.size(), temp.fly_time, fai_time, fai);
//	if (temp.fly_time > 68) fai_now = 35;
//	fai_now = fai_now * M_PI / 180.0;
//	// 推力在发射系下矢量
//	vector<double> Body_Direction{0,0,0}, Acce_P{ 0,0,0 };
// 	Body_Direction[0] = cos(fai_now);
//	Body_Direction[1] = 0;
//	Body_Direction[2] = sin(fai_now);
//
//	for (size_t i = 0; i < 3; ++i) {
//		Acce_P[i] = A_P * Body_Direction[i];
//	}
//
//	return Acce_P;
//}
//
//
//
//
//// 发射点天文经纬度、瞄准方位角及发射点高度
//const double lambda_T = 120.0;	// 天文经度 120 度
//const double B_T = 35.0;		// 天文纬度 35 度
//const double A_T = 45.0;		// 瞄准方位角 45 度
//const double Height = 1000.0; // 发射高度1000米
//
//// 定义rocket
//rocket test;
//
//int main() {
//	// 经纬度转地心大地直角坐标系下位置坐标
//	// double X_Dinxin, Y_Dinxin, Z_Dinxin;
//	double X_Dinxin0, Y_Dinxin0, Z_Dinxin0;	//地心大地直角坐标系下初始位置坐标
//	toECEF(lambda_T, B_T, Height, X_Dinxin0, Y_Dinxin0, Z_Dinxin0);
//	// double C_Dinxin[] = { X_Dinxin0, Y_Dinxin0, Z_Dinxin0 };	// 保存位置
//	double r = pow(X_Dinxin0 * X_Dinxin0 + Y_Dinxin0 * Y_Dinxin0 + Z_Dinxin0 * Z_Dinxin0, 0.5);
//	// 定义rocket
//	// rocket test;
//	test.Height = test.Height0 = Height; 
//	test.r_fai_lambda = { r, B_T, lambda_T };
//	test.xyz_dixin = test.xyz_dixin0 = { X_Dinxin0, Y_Dinxin0, Z_Dinxin0 };
//	// 初始重量KG
//	test.Weight = 44000 + 18800 + 5700 + 1000;
//
//	// 地球自转角速度在地心直角坐标系和发射系下的分量
//	vector<double> omega_s = { 0,0,omega };
//	vector<double> omega_g = { 0, 0, 0 };
//	g_s(omega_g, omega_s, 0);
//
//	// 额定推力参数
//
//	// 参数赋值
//	P1.Time = { 0.0, 1.25, 10.75, 30.75, 60.75, 65.75 };
//	P2.Time = { 0.0, 1.25, 20.75, 40.75, 53.15 };
//	for (size_t i = 0; i < P2.Time.size(); i++) {
//		P2.Time[i] += 65.75;
//		// P2.Time[i] += *(P1.Time.end() - 1);
//	}
//	P3.Time = { 0.0, 1.25, 10.75, 30.75, 47.5 };
//	for (size_t i = 0; i < P3.Time.size(); i++) {
//		P2.Time[i] += 65.75 + 53.15;
//	}
//	// 额定推力1-2-3
//	P1.P0 = { 0.0, 610.0, 860.0, 950.0, 820.0, 6.0 };
//	P2.P0 = { 0.0, 450.0, 740.0, 610.0, 10.0 };
//	P3.P0 = { 0.0, 160.0, 250.0, 300.0, 4.0 };
//	// 秒耗量1-2-3
//	P1.M_per_s = { 0.0, 250.0, 340.0, 380.0, 330.0, 6.0 };
//	P2.M_per_s = { 0.0, 250.0, 250.0, 210.0, 5.0 };
//	P3.M_per_s = { 0.0, 50.0, 80.0, 100.0, 2.0 };
//	// 参数数量
//	P1.n = 6;
//	P2.n = 5;
//	P3.n = 5;
//	// 横截面积、喷口截面积
//	P1.S_hengjie = 3.2; P1.S_penkou = 3.1;
//	P2.S_hengjie = 2.8; P2.S_penkou = 2.5;
//	P3.S_hengjie = 3.2; P3.S_penkou = 3.1;
//
//	double h = 0.02;	// 步长
//
//	// 所需变量
//	int n = 8;	// y[]元素数量
//	double y[] = { test.fly_time, test.Weight, test.xyz_fs[0],test.xyz_fs[1],test.xyz_fs[2], test.v_xyz_fs[0], test.v_xyz_fs[1], test.v_xyz_fs[2] };
//	double dy[] = { 0,0,0,0,0,0,0,0 };
//	double yc[8],y1[8]; // 过程变量
//	double Longitude, Latitude;
//
//	// 起飞前发动机工作
//	while (Acce_count(Acce_p(test, P1)) < Acce_count(Acce_g(test))) {
//		double guance1 = Acce_count(Acce_p(test, P1));
//		double guance2 = Acce_count(Acce_g(test));
//		test.work_time += h;
//	}
//	printf("%fs时开始起飞\n", test.work_time);
//
//	while (test.Height > 10) {
//
//		//// P_parameter P;
//		//if (test.work_time <= P1.Time[P1.n - 1]) {
//
//
//		//}
//		//else if (P1.Time[P1.n - 1] <= test.work_time <= P2.Time[P2.n - 1]) {
//		//	// P = P2;
//		//}
//		//else {
//
//		//}
//		//
//
//		runk(n, h, y, dy, yc, y1);
//
//		// 更新时间
//		test.fly_time += h;
//		test.work_time += h;
//
//		// 更新rocket状态
//		test.xyz_fs = { y[2],y[3],y[4] };
//		test.v_xyz_fs = { y[5],y[6],y[7] };
//
//		g_s(test.v_xyz_fs, test.xyz_dixin, 1);
//		test.Weight = y[1];
//
//		//double N = a / sqrt(1 - e2 * sin(test.r_fai_lambda[1]) * sin(test.r_fai_lambda[1]));
//		//test.Height = Acce_count(test.xyz_dixin) - N;
//
//
//		// test.Height = Acce_count(test.xyz_dixin) - 
//		// ECEFtoWGS84(Longitude, Latitude, Height, yc[2], yc[3], yc[4]);
//
//	}
//
//
//	//cout << "发射系下落点坐标为：(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//
//
//	system("pause");
//	return 0;
//}
//
//
//// 地球自转角速度  
//double A_zizhuan[] = { 0, 0, omega };
//
//// 额定推力
//// 时间1-2-3
//double P0_1_Time[] = { 0.0, 1.25, 10.75, 30.75, 60.75, 65.75 };
//double P0_2_Time[] = { 0.0, 1.25, 20.75, 40.75, 53.15 };
//double P0_3_Time[] = { 0.0, 1.25, 10.75, 30.75, 47.5 };
//// 额定推力1-2-3
//double P0_1[] = { 0.0, 610.0, 860.0, 950.0, 820.0, 6.0 };
//double P0_2[] = { 0.0, 450.0, 740.0, 610.0, 10.0 };
//double P0_3[] = { 0.0, 160.0, 250.0, 300.0, 4.0 };
//// 秒耗量1-2-3
//double M_per_s_1[] = { 0.0, 250.0, 340.0, 380.0, 330.0, 6.0 };
//double M_per_s_2[] = { 0.0, 250.0, 250.0, 210.0, 5.0 };
//double M_per_s_3[] = { 0.0, 50.0, 80.0, 100.0, 2.0 };
//// 数量
//int n1 = 6;
//int n2 = 5;
//int n3 = 5;
//
////// 飞行程序角
////// 时间1-2
////double fai_time[] = { 0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 64.0, 66.5, 68.0 };
////// double alpha_1_time[] = { 0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 64.0 };
////// double alpha_2_time[] = { 66.5, 68.0 };  // 后续等程序角飞行
////// 飞行程序角1-2
////double fai[] = { 90.0, 90.0, 69.0, 59.0, 50.0, 40.0, 36.0, 36.0, 35.0 };
//// double alpha_1[] = { 90.0, 90.0, 69.0, 59.0, 50.0, 40.0, 36.0 };
//// double alpha_2[] = { 36.0, 35.0 };
//
//// 空气动力系数
//// 马赫数1-2-3
//double Ma_1[] = { 0.2, 1.05, 4.0 };
//double Ma_2[] = { 3.0, 12.0 };
//double Ma_3[] = { 3.0, 12.0 };
//// Cx 1-2-3
//double Cx_1[] = { 0.15, 0.34, 0.16 };
//double Cx_2[] = { 0.18, 0.10 };
//double Cx_3[] = { 0.18, 0.10 };
//// Cy 1-2-3
//double Cy_1 = 0.039;
//double Cy_2 = 0.035;
//double Cy_3 = 0.035;
//// 再入阻力系数
//double Cx_zairu = 0.025;
//
//
//// 计算两个三维向量的叉乘  
//double crossProduct(double a[], double b[], double result[]) {
//	result[0] = a[1] * b[2] - a[2] * b[1];
//	result[1] = a[2] * b[0] - a[0] * b[2];
//	result[2] = a[0] * b[1] - a[1] * b[0];
//	return 0;
//}
//
//
////两点插值(增序)
//double edcz(int n, double x, vector<double> a, vector<double> b) {
//	// n--数组元素个数，x--带入点，a[]--x轴元素数组，b[]--y轴元素数组
//
//	//if (x < a[1]) {
//	//	return b[1];    
//	//}
//	int i;
//	for (i = 1; i < n; i++) {
//		if (x <= a[i]) {
//			break;
//		}
//	}
//	if (i == n) {
//		i--;  // 如果x大于所有给定的a值，则选择最后一个区间进行插值  
//	}
//	double fxy = b[i - 1] + (b[i] - b[i - 1]) * (x - a[i - 1]) / (a[i] - a[i - 1]);
//	return fxy;
//}
//
//// 椭球下天文经纬度-->地心大地直角坐标系
//void toECEF(double lambda, double B, double Height, double& X_Dinxin, double& Y_Dinxin, double& Z_Dinxin) {
//	// lambda天文经度，B天文纬度，
//
//	// 计算卯酉圈曲率半径N  
//	double N = a / sqrt(1 - e2 * sin(B) * sin(B));
//
//	// 计算X, Y, Z  
//	X_Dinxin = (N + Height) * cos(B) * cos(lambda);
//	Y_Dinxin = (N + Height) * cos(B) * sin(lambda);
//	Z_Dinxin = ((a * a / b / b) * N + Height) * sin(B);
//}
//
//double WGS84toECEF(double& latitude, double& longitude, double& height, double& x, double& y, double& z)
//{
//	double a = 6378137.0;
//	double b = 6356752.31424518;
//	double E = (a * a - b * b) / (a * a);
//	double COSLAT = cos(latitude * M_PI / 180);
//	double SINLAT = sin(latitude * M_PI / 180);
//	double COSLONG = cos(longitude * M_PI / 180);
//	double SINLONG = sin(longitude * M_PI / 180);
//	double N = a / (sqrt(1 - E * SINLAT * SINLAT));
//	double NH = N + height;
//	x = NH * COSLAT * COSLONG;
//	y = NH * COSLAT * SINLONG;
//	z = (b * b * N / (a * a) + height) * SINLAT;
//	return 0;
//}
//
//
//double ECEFtoWGS84(double& Longitude, double& Latitude, double& Height, double x, double y, double z)
//
//{
//	double a, b, c, d;
//	double p, q;
//	double N;
//	a = 6378137.0;
//	b = 6356752.31424518;
//	c = sqrt(((a * a) - (b * b)) / (a * a));
//	d = sqrt(((a * a) - (b * b)) / (b * b));
//	p = sqrt((x * x) + (y * y));
//	q = atan2((z * a), (p * b));
//	Longitude = atan2(y, x);
//	Latitude = atan2((z + (d * d) * b * pow(sin(q), 3)), (p - (c * c) * a * pow(cos(q), 3)));
//	N = a / sqrt(1 - ((c * c) * pow(sin(Latitude), 2)));
//	Height = (p / cos(Latitude)) - N;
//	Longitude = Longitude * 180.0 / M_PI;
//	Latitude = Latitude * 180.0 / M_PI;
//	return 0;
//}
//// 计算基于椭球经纬度的点到地心的距离  
//double calculateDistanceToCenter(double latitudeDegrees) {
//	// 输入为角度
//	double latitudeRadians = latitudeDegrees * M_PI / 180.0;	// 将纬度转换为弧度  
//	// 计算距离  
//	double numerator = a * b;
//	double denominator = sqrt(b * b * cos(latitudeRadians) * cos(latitudeRadians) + a * a * sin(latitudeRadians) * sin(latitudeRadians));
//	return numerator / denominator;
//}
//
//
//
//// 地心大地直角坐标-->椭球下天文经纬度  
//void xyz_to_latlon(double X_Dinxin, double Y_Dinxin, double Z_Dinxin, double& lat, double& lon) {
//	//lat = lat * M_PI / 180.0;
//	//lon = lon * M_PI / 180.0;
//	lon = atan2(Y_Dinxin, X_Dinxin);
//	double p = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin);
//	double theta = atan2(Z_Dinxin * a, p * b);
//	lat = atan2(Z_Dinxin + e2 * b * pow(sin(theta), 3), p - e2 * a * pow(cos(theta), 3));
//	lat = lat * 180.0 / M_PI;  // 转换为度  
//	lon = lon * 180.0 / M_PI;  // 转换为度  
//}
//
////大气    xp为大气压力，ru为大气密度,a为音速
//void daqi_US(double h_daqi, double* xp, double* ru, double* a)
//{
//	//多用一维,第一个为0
//	int i;
//	if (h_daqi < 0.0) h_daqi = 0.0;
//	//	if(flag_dlh && h_daqi<=31000.0) {              //干扰大气模型
//	//	  double zm=6356.766*h_daqi/(6356766.0+h_daqi);
//	//	  daqi3(zm);
//	//	  return;
//	//	}
//	double hb[9] = { 0.0,0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852 };
//	double tb[9] = { 0.0,288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.95 };
//	double pb[9] = { 0.0,1.01325E+5,2.2632E+4,5.4748E+3,8.6801E+2,1.10900E+2,6.6938E+1,3.9564E+0,3.7338E-1 };  //放大100倍
//	double td[8] = { 0.0,-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0 };
//	double zg = 6356.766 * h_daqi / (6356766.0 + h_daqi);
//	double work = 34.16319474;               //9.80665*28.9644/8.31432
//	double dt;
//	if (zg > 84.852e0)
//	{
//		dt = 186.95e0;
//		*xp = 0.e0;
//		*ru = 0.e0;
//		*a = 274.0992e0;
//		return;
//	}
//	for (i = 1; i <= 8; i++) {
//		if (zg >= hb[i] && zg < hb[i + 1]) {
//			dt = tb[i] + td[i] * (zg - hb[i]);
//			if (td[i] != 0.0) {
//				*xp = pb[i] * exp(work / td[i] * log(tb[i] / dt));
//			}
//			else {
//				*xp = pb[i] * exp(-work * (zg - hb[i]) / tb[i]);
//			}
//			double tmp = *xp;
//			*ru = tmp * 3.483676356e-3 / dt;
//			*a = 20.0468 * sqrt(dt);
//		}
//	}
//	return;
//}
//
//// 微分方程
//void rfun(int n, double yc[], double dy[]) {
//	double Time_all = yc[0];
//	double Weight = yc[1];
//	double C_Launch[] = { yc[2],yc[3],yc[4] };
//	double V_Launch[] = { yc[5],yc[6],yc[7] };
//	double h = 0.02;	// 步长
//
//	P_parameter P;
//	if (test.work_time <= 65.75) {
//		P = P1;
//	}
//	else if (65.75 < test.work_time <= 65.75 + 53.15) {
//		P = P2;
//	}
//	else if (65.75 + 53.15 < test.work_time <= 65.75 + 53.15 + 47.5) {
//		P = P3;
//	}
//	else {
//		P.M_per_s = { 0 };
//		P.n = 1;
//		P.P0 = { 0 };
//		P.S_hengjie = 0;
//		P.S_penkou = 0;
//		P.Time = { 0 };
//	}
//	vector<double> acce_g = Acce_g(test);
//
//	vector<double> acce_p = Acce_p(test, P);
//
//	for (int i = 0; i < 3; i++) {
//		dy[i + 2] = test.v_xyz_fs[i];
//		dy[i + 5] = acce_g[i] + acce_p[i];
//	}
//	dy[1] = edcz(P.n, test.work_time, P.Time, P.M_per_s);
//	// dy[0] = h;
//
//
//
//
//
//
//	//// int SegmentLable;	// 段标
//	///*飞行分段--根据不同阶段受力情况不同分段
//	//0-未起飞
//	//1-起飞后发动机一级工作阶段
//	//2-起飞后发动机二级工作阶段
//	//3-起飞后发动机三级工作阶段
//	//4-发动机关机后80km以下飞行阶段
//	//5-80km以上，无空气动力阶段
//	//6-80km以下，再入段
//	//7-落地终止*/
//
//
//	//// xp为大气压力，ru为大气密度, a为音速
// ////double* xp;
// ////double* ru;
// ////double* a;
//	//double h0 = 1000;
//	//double g, p, ru, a, p0, ru0, a0, P0;
//	//double S_hengjie, S_penkou;
//	//double cross_Keshi[3], cross_QianLian[3], Acce_QianLian[3];
//	//daqi_US(h0, &p0, &ru0, &a0);
//	//double v;
//	//double Ma_Now, P;
//	//double A_P;
//	//double Cx, Cy, X, Y, fai_now, Time_fly = 3.25;
//	//// 计算发动机推力加速度
//
//	//double Acce_X[3], Acce_Y[3], Acce_P[3], Acce_Yinli[3], Acce_All[3];
//	//double Body_Direction[3], X_Direction[3], Y_Direction[3];
//
//	//// 起飞后一级发动机工作阶段		
//	//S_hengjie = 3.2;
//	//S_penkou = 3.1;
//	//// 计算引力加速度
//	//double X_Dinxin0, Y_Dinxin0, Z_Dinxin0;	//地心大地直角坐标系下初始位置坐标
//	//toECEF(120, 35, 1000, X_Dinxin0, Y_Dinxin0, Z_Dinxin0);
//	//double Longitude, Latitude, Height, M_per_s;
//	//ECEFtoWGS84(Longitude, Latitude, Height, (X_Dinxin0 + yc[2]), (Y_Dinxin0 + yc[3]), (Z_Dinxin0 + yc[4]));
//
//	//double X_Dinxin, Y_Dinxin, Z_Dinxin;
//
//
//	//toECEF(Latitude, Longitude, Height, X_Dinxin, Y_Dinxin, Z_Dinxin);
//	//double C_Dinxin[] = { X_Dinxin, Y_Dinxin, Z_Dinxin };	// 保存位置
//	//double C_Dixin[] = { X_Dinxin, Y_Dinxin, Z_Dinxin };
//	//double r;
//	//r = sqrt(pow(C_Launch[0] + C_Dixin[0], 2) + pow(C_Launch[1] + C_Dixin[1], 2) + pow(C_Launch[2] + C_Dixin[2], 2));    // 地心矢径长
//	//X_Dinxin = X_Dinxin0 + C_Launch[0];
//	//Y_Dinxin = Y_Dinxin0 + C_Launch[1];
//	//Z_Dinxin = Z_Dinxin0 + C_Launch[2];
//	//xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Latitude, Longitude);
//	//Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Latitude);
//	//// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//g = -fM / (r * r);               // 引力加速度大小
//	//for (size_t i = 0; i < 3; i++) {
//	//	Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//}
//
//	//// 计算科氏加速度
//	//crossProduct(V_Launch, A_zizhuan, cross_Keshi);
//	//double Acce_Keshi[3];
//	//for (size_t i = 0; i < 3; ++i) {
//	//	Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//	if (v < 1) Acce_Keshi[i] = 0;
//	//}
//
//	//// 计算牵连加速度
//	//crossProduct(A_zizhuan, C_Dinxin, cross_QianLian);
//	//crossProduct(A_zizhuan, cross_QianLian, Acce_QianLian);
//	//for (size_t i = 0; i < 3; ++i) {
//	//	Acce_QianLian[i] = -cross_QianLian[i];
//	//}
//
//	//// 计算发动机推力加速度
//
//	//// 根据段标选定横截面积、喷口面积
//
//	//// 根据当前位置高度计算大气压强p、空气密度ru、音速a
//	//daqi_US(Height, &p, &ru, &a);
//	//if (Time_all < 180) {
//	//	P0 = edcz(n1, Time_all, P0_1_Time, P0_1) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// 推力大小转推力加速度
//	//	A_P = P / Weight;
//	//	// 发射系下推力方向
//	//	// 
//	//	M_per_s = edcz((sizeof(P0_1_Time) / sizeof(P0_1_Time[0])), Time_all, P0_1_Time, M_per_s_1);
//	//	fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
//	//	if (Time_all > 68) fai_now = 35;
//	//	fai_now = fai_now * M_PI / 180.0;
//	//	Body_Direction[0] = cos(fai_now);
//	//	Body_Direction[1] = 0;
//	//	Body_Direction[2] = sin(fai_now);
//
//	//	for (size_t i = 0; i < 3; ++i) {
//	//		Acce_P[i] = A_P * Body_Direction[i];
//	//	}
//	//	/*Acce_P[0] = F_P * cos(fai_now);
//	//	Acce_P[1] = F_P * sin(fai_now);
//	//	Acce_P[2] = 0;*/
//	//}
//	//else double Acce_P[] = { 0,0,0 };
//
//	//if (Height < 80000) {
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
//	//	for (size_t i = 0; i < 3; ++i) {
//	//		X_Direction[i] = -V_Launch[i] / (v * sqrt(3));
//	//	}
//	//	//// y_direction
//	//	Y_Direction[0] = -V_Launch[1] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = V_Launch[0] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
//	//	Y_Direction[1] = 0;
//	//	for (size_t i = 0; i < 3; ++i) {
//	//		Acce_X[i] = X * X_Direction[i];
//	//		Acce_Y[i] = Y * Y_Direction[i];
//	//		//Acce_X[i] = 0;
//	//		//Acce_Y[i] = 0;
//	//	}
//	//}
//	//else {
//	//	for (size_t i = 0; i < 3; ++i) { Acce_X[i] = 0; Acce_Y[i] = 0; }
//	//	M_per_s = 0;
//	//}
//
//	//// 合加速度
//	//for (size_t i = 0; i < 3; ++i) {
//	//	Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i] + Acce_P[i] + Acce_X[i] + Acce_Y[i];
//	//	V_Launch[i] += Acce_All[i] * h;
//	//	C_Launch[i] += V_Launch[i] * h;
//	//}
//
//	//dy[0] = h;
//	//dy[1] = M_per_s;
//	//dy[2] = yc[5];
//	//dy[3] = yc[6];
//	//dy[4] = yc[7];
//	//dy[5] = Acce_All[0];
//	//dy[6] = Acce_All[1];
//	//dy[7] = Acce_All[2];
//
//	return;
//}
//
////龙格库塔法
//void runk(int n, double h, double y[], double dy[], double yc[], double y1[])
//{
//	// n-元素个数，h-步长，y[]-函数值，dy[]-微分，yc[]、y1[]-过程量
//	//	extern void rf(int n,double y[25],double dy[25],double tj);
//	int i, k;
//	double a[5];
//	a[1] = h / 2.0;  a[2] = a[1];  a[3] = h; a[4] = h;
//	for (i = 0; i < n; i++)  y1[i] = y[i];	// 复制上一个步长的函数值
//
//	for (k = 1; k <= 3; k++)
//	{
//		for (i = 0; i < n; i++)
//		{
//			yc[i] = y1[i] + a[k] * dy[i];
//			y[i] = y[i] + a[k + 1] * dy[i] / 3.0;
//		}
//
//		rfun(n, yc, dy);
//	}
//	for (i = 0; i < n; i++)  y[i] = y[i] + a[1] * dy[i] / 3.0;
//	rfun(n, y, dy);
//	return;
//}
