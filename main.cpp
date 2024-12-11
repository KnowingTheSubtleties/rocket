//#include <iostream>
//#include <math.h>
//#include <cmath>  
//#include <vector>
//
//using namespace std;
//
//const double M_PI = 3.14159265358979323846;
//// ��������˹���������  
//const double a = 6378245.0; // ��������˹������ĳ�����  
//const double f = 1.0 / 298.3; // ����
//const double omega = 7.292115e-5; // ������ת���ٶ�
//const double fM = 3.98620e14; // ������������
//const double b = a * (1 - f);  // �̰���  
//const double e2 = (2 - f) * f;  // ��һƫ���ʵ�ƽ��  
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
//	double Time_all = 0.0; // ��ʼʱ��
//	double Time_fly = 0.0; // ����ʱ��
//
//	// ��������ľ�γ�ȡ���׼��λ�Ǽ������߶�
//	double lambda_T = 35.0;	// ���ľ��� 35 ��
//	double B_T = 120.0;		// ����γ�� 120 ��
//	double A_T = 45.0;		// ��׼��λ�� 45 ��
//	double Height = 1000.0; // ����߶�1000��
//
//
//
//	// ��γ��ת���Ĵ��ֱ������ϵ��λ������
//	double X_Dinxin, Y_Dinxin, Z_Dinxin;
//	double X_Dinxin0, Y_Dinxin0, Z_Dinxin0;	//���Ĵ��ֱ������ϵ�³�ʼλ������
//	toECEF(lambda_T, B_T, Height, X_Dinxin0, Y_Dinxin0, Z_Dinxin0);
//	std::vector<double> C_Dinxin = { X_Dinxin0, Y_Dinxin0, Z_Dinxin0 };	// ����λ��
//	return 0;
//	//// ����ϵ�µ�ǰλ��
//	//// double X_Launch = 0.0, Y_Launch = 0.0, Z_Launch = 0.0;
//	//vector<double> C_Launch = { 0.0, 0.0, 0.0 };
//
//	//// ����ϵ���ٶ�m/s
//	//// double VX_Launch = 0.0, VY_Launch = 0.0, VZ_Launch = 0.0;
//	//vector<double> V_Launch = { 0.0, 0.0, 0.0 };
//
//	//// ��ʼ����KG
//	//double Weight = 44000 + 18800 + 5700 + 1000;
//
//	//// �������
//	////int n = 8;	// y[]Ԫ������
//	////double y[] = { Time_fly, Height, C_Launch[0],C_Launch[1],C_Launch[2], V_Launch[0], V_Launch[1], V_Launch[2] };
//
//
//
//	//// -------------------��ʼ��ѭ��-------------------
//	//double dt = 0.02;		// ����
//	//double P0;				// �����������
//	//double M_per_s;			// �����
//	//double fai_now;			// ���г����
//	//double Cx, Cy;			// ����ϵ��������ϵ��
//	//double X, Y;			// ����������
//	//double r, r0;			// ����ʸ������ʼ����ʸ��
//	//double g;				// �������ٶȴ�С
//	//double omega_x, omega_y, omega_z;	// ��ת���ٶȷ���
//	//// omega_x = omega * cos();
//	//vector<double> Acce_Yinli(3);	// �������ٶ�
//	//vector<double> Body_Direction = { 1,0,0 };	// ����ϵ�µ��巽��ʸ�������ڱ�ʾ��������
//	//vector<double> X_Direction = { 1,0,0 };		// ��������
//	//vector<double> Y_Direction = { 0,1,0 };		// ��������
//
//	//r0 = sqrt(C_Dinxin[0] * C_Dinxin[0] + C_Dinxin[1] * C_Dinxin[01] + C_Dinxin[2] * C_Dinxin[2]);    // ����ʸ����
//	//
//	//// δ��ɽ׶Σ�һ������������
//	//while (SegmentLable == 0) {
//	//	// ʱ���ۼ�
//	//	Time_all += 0.02;
//	//	// ���������������
//	//	P0 = edcz((sizeof(P0_1_Time)/sizeof(P0_1_Time[0])), Time_all, P0_1_Time, P0_1) * 1000;
//	//	M_per_s = edcz((sizeof(P0_1_Time) / sizeof(P0_1_Time[0])), Time_all, P0_1_Time, M_per_s_1);
//	//	// ����rocket����
//	//	Weight -= dt * M_per_s;
//	//	// �����������ٶ�
//	//	g = -fM / (r0 * r0);               // �������ٶȴ�С
//	//	for (size_t i = 0; i < C_Dinxin.size(); ) {
//	//		Acce_Yinli[i] = g * C_Dinxin[i] / r0;
//	//		i++;
//	//	}
//
//	//	if (P0 / Weight >= sqrt(Acce_Yinli[0] * Acce_Yinli[0] + Acce_Yinli[1] * Acce_Yinli[1] + Acce_Yinli[2] * Acce_Yinli[2])) {
//	//		SegmentLable = 1;
//	//		cout << Time_all << "s\t��ʼ���\n" << endl;
//	//	}
//	//}
//
//	//// xpΪ����ѹ����ruΪ�����ܶ�, aΪ����
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
//	//// ���㷢�����������ٶ�
//
//	//vector<double> Acce_X(3), Acce_Y(3);
//	//vector<double> Acce_P(3);
//	//vector<double> Acce_All(3);
//	//// ��ɺ�һ�������������׶�		
//	//S_hengjie = 3.2;
//	//S_penkou = 3.1;
//	//while (SegmentLable == 1) {
//	//	// ʱ���ۼ�
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_all, P0_1_Time, M_per_s_1);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// �����������ٶ�
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // ����ʸ����
//	//	X_Dinxin = X_Dinxin0 + C_Launch[0];
//	//	Y_Dinxin = Y_Dinxin0 + C_Launch[1];
//	//	Z_Dinxin = Z_Dinxin0 + C_Launch[2];
//	//	xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Lat, Lon);
//	//	Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Lat);
//	//	// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // �������ٶȴ�С
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//	//	// ������ϼ��ٶ�
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// ����ǣ�����ٶ�
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// ���㷢�����������ٶ�
//
//	//	// ���ݶα�ѡ����������������
//
//	//	// ���ݵ�ǰλ�ø߶ȼ������ѹǿp�������ܶ�ru������a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	P0 = edcz(n1, Time_all, P0_1_Time, P0_1) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// ������Сת�������ٶ�
//	//	A_P = P / Weight;
//	//	// ����ϵ����������
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
//	//	// �����������
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// ������������С
//	//	Ma_Now = v / a;	// �������С
//	//	Cx = edcz(sizeof(Cx_1) / sizeof(Cx_1[0]), Ma_Now, Ma_1, Cx_1);	// 
//	//	Cy = Cy_1;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// ת���ٶ�
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
//	//	// �ϼ��ٶ�
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
//	//	// �ٶȴ���0 �ۼƷ���ʱ��
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75) {
//	//		SegmentLable = 2;
//	//		Weight = 18800 + 5700 + 1000;
//	//		cout << "һ��������������������ǰ���и߶ȣ�" << (int)Height << "m\n" << "�ۼ�ʱ�䣺" << Time_all << "s" << endl;
//	//		cout << "��ǰ����ϵ������Ϊ��(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 2) {
//	//	// ʱ���ۼ�
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_fly, P0_2_Time, M_per_s_2);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// �����������ٶ�
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // ����ʸ����
//	//	X_Dinxin = X_Dinxin0 + C_Launch[0];
//	//	Y_Dinxin = Y_Dinxin0 + C_Launch[1];
//	//	Z_Dinxin = Z_Dinxin0 + C_Launch[2];
//	//	xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Lat, Lon);
//	//	Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Lat);
//	//	// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // �������ٶȴ�С
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// ������ϼ��ٶ�
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// ����ǣ�����ٶ�
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// ���㷢�����������ٶ�
//
//	//	// ���ݶα�ѡ����������������
//	//	S_hengjie = 2.8;
//	//	S_penkou = 2.5;
//	//	// ���ݵ�ǰλ�ø߶ȼ������ѹǿp�������ܶ�ru������a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	// 
//	//	// 
//	//	P0 = edcz(n1, Time_all - 65.75, P0_2_Time, P0_2) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// ������Сת�������ٶ�
//	//	A_P = P / Weight;
//	//	// ����ϵ����������
//	//	// 
//	//	fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
//	//	if (Time_all > 68) fai_now = 35;
//	//	fai_now = fai_now * 3.1415926 / 180.0;
//	//	Body_Direction = { cos(fai_now) ,sin(fai_now) ,0 };
//	//	for (size_t i = 0; i < Acce_P.size(); ++i) {
//	//		Acce_P[i] = A_P * Body_Direction[i];
//	//	}
//
//	//	// �����������
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// ������������С
//	//	Ma_Now = v / a;	// �������С
//	//	Cx = edcz(sizeof(Cx_2) / sizeof(Cx_2[0]), Ma_Now, Ma_2, Cx_2);	// 
//	//	Cy = Cy_2;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// ת���ٶ�
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
//	//	// �ϼ��ٶ�
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
//	//	// �ٶȴ���0 �ۼƷ���ʱ��
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75 + 53.15) {
//	//		SegmentLable = 3;
//	//		Weight = 5700 + 1000;
//	//		cout << "����������������������ǰ���и߶ȣ�" << (int)Height << "m\n" << "�ۼ�ʱ�䣺" << Time_all << "s" << endl;
//	//		cout << "��ǰ����ϵ������Ϊ��(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 3) {
//	//	// ʱ���ۼ�
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//	//	M_per_s = edcz(n1, Time_fly, P0_3_Time, M_per_s_3);
//
//	//	Weight -= dt * M_per_s;
//
//	//	// �����������ٶ�
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // ����ʸ����
//	//	Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // �������ٶȴ�С
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// ������ϼ��ٶ�
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// ����ǣ�����ٶ�
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//	//	// ���㷢�����������ٶ�
//
//	//	// ���ݶα�ѡ����������������
//	//	S_hengjie = 3.2;
//	//	S_penkou = 3.1;
//	//	// ���ݵ�ǰλ�ø߶ȼ������ѹǿp�������ܶ�ru������a
//	//	daqi_US(Height, &p, &ru, &a);
//
//	//	// 
//	//	// 
//	//	P0 = edcz(n1, Time_all - 65.75 - 53.15, P0_3_Time, P0_3) * 1000;
//	//	P = P0 + S_penkou * (1 - p / p0);
//	//	// ������Сת�������ٶ�
//	//	A_P = P / Weight;
//	//	// ����ϵ����������
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
//	//	// �����������
//	//	// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// ������������С
//	//	Ma_Now = v / a;	// �������С
//	//	Cx = edcz(sizeof(Cx_3) / sizeof(Cx_3[0]), Ma_Now, Ma_3, Cx_3);	// 
//	//	Cy = Cy_3;	// 
//
//	//	X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
//	//	// ת���ٶ�
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
//	//	// �ϼ��ٶ�
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
//	//	// �ٶȴ���0 �ۼƷ���ʱ��
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Time_all >= 65.75 + 53.15 + 47.5) {
//	//		SegmentLable = 4;
//	//		Weight = 1000;
//	//		cout << "����������������������ǰ���и߶ȣ�" << (int)Height << "m\n" << "�ۼ�ʱ�䣺" << Time_all << "s" << endl;
//	//		cout << "��ǰ����ϵ������Ϊ��(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//	//while (SegmentLable == 4) {
//	//	// ʱ���ۼ�
//	//	Time_all += 0.02;
//	//	// V_Launch = { VX_Launch, VY_Launch, VZ_Launch };
//
//	//	// �����������ٶ�
//	//	///----------
//	//	r = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2));    // ����ʸ����
//	//	Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
//	//	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
//	//	g = -fM / (r * r);               // �������ٶȴ�С
//	//	for (size_t i = 0; i < C_Dinxin.size(); i++) {
//	//		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
//	//	}
//
//
//	//	// ������ϼ��ٶ�
//	//	cross_Keshi = crossProduct(V_Launch, A_zizhuan);
//	//	vector<double> Acce_Keshi(cross_Keshi.size());
//	//	for (size_t i = 0; i < cross_Keshi.size(); ++i) {
//	//		Acce_Keshi[i] = 2 * cross_Keshi[i];
//	//		if (v < 1) Acce_Keshi[i] = 0;
//	//	}
//
//	//	// ����ǣ�����ٶ�
//	//	cross_QianLian = crossProduct(A_zizhuan, C_Dinxin);
//	//	Acce_QianLian = crossProduct(A_zizhuan, cross_QianLian);
//	//	for (size_t i = 0; i < cross_QianLian.size(); ++i) {
//	//		Acce_QianLian[i] = -cross_QianLian[i];
//	//	}
//
//
//	//	// �ϼ��ٶ�
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
//	//	// �ٶȴ���0 �ۼƷ���ʱ��
//	//	if (v > 0) {
//	//		Time_fly += 0.02;
//	//	}
//
//	//	//if ();
//	//	if (Height <= 80000) {
//	//		SegmentLable = 5;
//	//		cout << "���������������У���ǰ���и߶ȣ�" << (int)Height << "m\n" << "�ۼ�ʱ�䣺" << Time_all << "s" << endl;
//	//		cout << "��ǰ����ϵ������Ϊ��(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//	//	}
//	//}
//
//
//	////// ����΢�֣��ٶȡ����ٶȣ�
//	////double dy[] = { 0,0,0,0,0,0,0,0 };
//	////double yc[] = {0,0,0,0,0,0,0,0};
//	////double y1[] = { 0,0,0,0,0,0,0,0 };
//
//	////// ¡���������
//	////runk(n, Height, y, dy, yc, y1);
//
//	////// ���·���ϵ��λ�á��ٶ�
//
//
//	////// ���µ�ǰ���и߶�
//
//	//cout << "����ϵ���������Ϊ��(" << (int)C_Launch[0] << "," << (int)C_Launch[1] << "," << (int)C_Launch[2] << ")" << endl;
//
//
//	system("pause");
//
//}
//
//
//
