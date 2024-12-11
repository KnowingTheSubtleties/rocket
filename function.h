#include <iostream>
#include <math.h>
#include <cmath>  

using namespace std;

const double M_PI = 3.14159265358979323846;
// ��������˹���������  
const double a = 6378245.0; // ��������˹������ĳ�����  
const double f = 1.0 / 298.3; // ����
const double omega = 7.292115e-5; // ������ת���ٶ�
const double fM = 3.98620e14; // ������������
const double b = a * (1 - f);  // �̰���  
const double e2 = (2 - f) * f;  // ��һƫ���ʵ�ƽ��  

// ������ת���ٶ�  
double A_zizhuan[] = {0, 0, omega};

// �����
// ʱ��1-2-3
double P0_1_Time[] = { 0.0, 1.25, 10.75, 30.75, 60.75, 65.75 };
double P0_2_Time[] = { 0.0, 1.25, 20.75, 40.75, 53.15 };
double P0_3_Time[] = { 0.0, 1.25, 10.75, 30.75, 47.5 };
// �����1-2-3
double P0_1[] = { 0.0, 610.0, 860.0, 950.0, 820.0, 6.0 };
double P0_2[] = { 0.0, 450.0, 740.0, 610.0, 10.0 };
double P0_3[] = { 0.0, 160.0, 250.0, 300.0, 4.0 };
// �����1-2-3
double M_per_s_1[] = { 0.0, 250.0, 340.0, 380.0, 330.0, 6.0 };
double M_per_s_2[] = { 0.0, 250.0, 250.0, 210.0, 5.0 };
double M_per_s_3[] = { 0.0, 50.0, 80.0, 100.0, 2.0 };
// ����
int n1 = 6;
int n2 = 5;
int n3 = 5;

// ���г����
// ʱ��1-2
double fai_time[] = { 0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 64.0, 66.5, 68.0 };
// double alpha_1_time[] = { 0.0, 5.0, 10.0, 20.0, 30.0, 50.0, 64.0 };
// double alpha_2_time[] = { 66.5, 68.0 };  // �����ȳ���Ƿ���
// ���г����1-2
double fai[] = { 90.0, 90.0, 69.0, 59.0, 50.0, 40.0, 36.0, 36.0, 35.0 };
// double alpha_1[] = { 90.0, 90.0, 69.0, 59.0, 50.0, 40.0, 36.0 };
// double alpha_2[] = { 36.0, 35.0 };

// ��������ϵ��
// �����1-2-3
double Ma_1[] = { 0.2, 1.05, 4.0 };
double Ma_2[] = { 3.0, 12.0 };
double Ma_3[] = { 3.0, 12.0 };
// Cx 1-2-3
double Cx_1[] = { 0.15, 0.34, 0.16 };
double Cx_2[] = { 0.18, 0.10 };
double Cx_3[] = { 0.18, 0.10 };
// Cy 1-2-3
double Cy_1 = 0.039;
double Cy_2 = 0.035;
double Cy_3 = 0.035;
// ��������ϵ��
double Cx_zairu = 0.025;


// ����������ά�����Ĳ��  
double crossProduct(double a[],  double b[], double result[]) {
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
	return 0;
}


//�����ֵ(����)
double edcz(int n, double x, double a[], double b[]) {
	// n--����Ԫ�ظ�����x--����㣬a[]--x��Ԫ�����飬b[]--y��Ԫ������

	//if (x < a[1]) {
	//	return b[1];    
	//}
	int i;
	for (i = 1; i < n; i++) {
		if (x <= a[i]) {
			break;
		}
	}
	if (i == n) {
		i--;  // ���x�������и�����aֵ����ѡ�����һ��������в�ֵ  
	}
	double fxy = b[i - 1] + (b[i] - b[i - 1]) * (x - a[i - 1]) / (a[i] - a[i - 1]);
	return fxy;
}

double WGS84toECEF(double& latitude, double& longitude, double& height, double& x, double& y, double& z)
{
	double a = 6378137.0;
	double b = 6356752.31424518;
	double E = (a * a - b * b) / (a * a);
	double COSLAT = cos(latitude * M_PI / 180);
	double SINLAT = sin(latitude * M_PI / 180);
	double COSLONG = cos(longitude * M_PI / 180);
	double SINLONG = sin(longitude * M_PI / 180);
	double N = a / (sqrt(1 - E * SINLAT * SINLAT));
	double NH = N + height;
	x = NH * COSLAT * COSLONG;
	y = NH * COSLAT * SINLONG;
	z = (b * b * N / (a * a) + height) * SINLAT;
	return 0;
}


double ECEFtoWGS84(double& Longitude, double& Latitude, double& Height, double x, double y, double z)

{
	double a, b, c, d;
	double p, q;
	double N;
	a = 6378137.0;
	b = 6356752.31424518;
	c = sqrt(((a * a) - (b * b)) / (a * a));
	d = sqrt(((a * a) - (b * b)) / (b * b));
	p = sqrt((x * x) + (y * y));
	q = atan2((z * a), (p * b));
	Longitude = atan2(y, x);
	Latitude = atan2((z + (d * d) * b * pow(sin(q), 3)), (p - (c * c) * a * pow(cos(q), 3)));
	N = a / sqrt(1 - ((c * c) * pow(sin(Latitude), 2)));
	Height = (p / cos(Latitude)) - N;
	Longitude = Longitude * 180.0 / M_PI;
	Latitude = Latitude * 180.0 / M_PI;
	return 0;
}
// �����������γ�ȵĵ㵽���ĵľ���  
double calculateDistanceToCenter(double latitudeDegrees) {
	// ����Ϊ�Ƕ�
	double latitudeRadians = latitudeDegrees * M_PI / 180.0;	// ��γ��ת��Ϊ����  
	// �������  
	double numerator = a * b;
	double denominator = sqrt(b * b * cos(latitudeRadians) * cos(latitudeRadians) + a * a * sin(latitudeRadians) * sin(latitudeRadians));
	return numerator / denominator;
}

// ���������ľ�γ��-->���Ĵ��ֱ������ϵ
void toECEF(double lambda, double B, double Height, double& X_Dinxin, double& Y_Dinxin, double& Z_Dinxin) {
	// lambda���ľ��ȣ�B����γ�ȣ�

		lambda = lambda * M_PI / 180.0;
		B = B * M_PI / 180.0;
		// ����î��Ȧ���ʰ뾶N  
		double N = a / sqrt(1 - e2 * sin(B) * sin(B));

		// ����X, Y, Z  
		X_Dinxin = (N + Height) * cos(B) * cos(lambda);
		Y_Dinxin = (N + Height) * cos(B) * sin(lambda);
		Z_Dinxin = ((a * a / b / b) * N + Height) * sin(B);
}


// ���Ĵ��ֱ������-->���������ľ�γ��  
void xyz_to_latlon(double X_Dinxin, double Y_Dinxin, double Z_Dinxin, double& lat, double& lon) {
	//lat = lat * M_PI / 180.0;
	//lon = lon * M_PI / 180.0;
	lon = atan2(Y_Dinxin, X_Dinxin);
	double p = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin);
	double theta = atan2(Z_Dinxin * a, p * b);
	lat = atan2(Z_Dinxin + e2 * b * pow(sin(theta), 3), p - e2 * a * pow(cos(theta), 3));
	lat = lat * 180.0 / M_PI;  // ת��Ϊ��  
	lon = lon * 180.0 / M_PI;  // ת��Ϊ��  
}


void body_launch(double launch[], double body[], double fai, int flag) {
	// ����ϵ���꣬����ϵ���꣬�����ǣ�ת��flag
	if (flag == 0) { // ����ת����ϵ
		launch[0] = body[0] * cos(fai) - body[1] * sin(fai);
		launch[1] = body[0] * sin(fai) + body[1] * cos(fai);
		launch[2] = body[2];
	}
	else {
		body[0] = launch[0] * cos(fai) + launch[1] * sin(fai);
		body[1] = launch[1] * sin(fai) - launch[0] * cos(fai);
		body[2] = launch[2];
	}
}

//����    xpΪ����ѹ����ruΪ�����ܶ�,aΪ����
void daqi_US(double h_daqi, double* xp, double* ru, double* a)
{
	//����һά,��һ��Ϊ0
	int i;
	if (h_daqi < 0.0) h_daqi = 0.0;
	//	if(flag_dlh && h_daqi<=31000.0) {              //���Ŵ���ģ��
	//	  double zm=6356.766*h_daqi/(6356766.0+h_daqi);
	//	  daqi3(zm);
	//	  return;
	//	}
	double hb[9] = { 0.0,0.0,11.0,20.0,32.0,47.0,51.0,71.0,84.852 };
	double tb[9] = { 0.0,288.15,216.65,216.65,228.65,270.65,270.65,214.65,186.95 };
	double pb[9] = { 0.0,1.01325E+5,2.2632E+4,5.4748E+3,8.6801E+2,1.10900E+2,6.6938E+1,3.9564E+0,3.7338E-1 };  //�Ŵ�100��
	double td[8] = { 0.0,-6.5,0.0,1.0,2.8,0.0,-2.8,-2.0 };
	double zg = 6356.766 * h_daqi / (6356766.0 + h_daqi);
	double work = 34.16319474;               //9.80665*28.9644/8.31432
	double dt;
	if (zg > 84.852e0)
	{
		dt = 186.95e0;
		*xp = 0.e0;
		*ru = 0.e0;
		*a = 274.0992e0;
		return;
	}
	for (i = 1; i <= 8; i++) {
		if (zg >= hb[i] && zg < hb[i + 1]) {
			dt = tb[i] + td[i] * (zg - hb[i]);
			if (td[i] != 0.0) {
				*xp = pb[i] * exp(work / td[i] * log(tb[i] / dt));
			}
			else {
				*xp = pb[i] * exp(-work * (zg - hb[i]) / tb[i]);
			}
			double tmp = *xp;
			*ru = tmp * 3.483676356e-3 / dt;
			*a = 20.0468 * sqrt(dt);
		}
	}
	return;
}

// ΢�ַ���
void rfun(int n, double yc[], double dy[]) {
	double Time_all = yc[0];
	double Weight = yc[1];
	double C_Launch[] = { yc[2],yc[3],yc[4] };
	double V_Launch[] = { yc[5],yc[6],yc[7] };
	double h = 0.02;	// ����

	// int SegmentLable;	// �α�
	/*���зֶ�--���ݲ�ͬ�׶����������ͬ�ֶ�
	0-δ���
	1-��ɺ󷢶���һ�������׶�
	2-��ɺ󷢶������������׶�
	3-��ɺ󷢶������������׶�
	4-�������ػ���80km���·��н׶�
	5-80km���ϣ��޿��������׶�
	6-80km���£������
	7-�����ֹ*/


	// xpΪ����ѹ����ruΪ�����ܶ�, aΪ����
 //double* xp;
 //double* ru;
 //double* a;
	double h0 = 1000;
	double g, p, ru, a, p0, ru0, a0, P0;
	double S_hengjie, S_penkou;
	double cross_Keshi[3], cross_QianLian[3], Acce_QianLian[3];
	daqi_US(h0, &p0, &ru0, &a0);
	double v;
	double Ma_Now, P;
	double A_P;
	double Cx, Cy, X, Y, fai_now, Time_fly = 3.25;
	// ���㷢�����������ٶ�

	double Acce_X[3], Acce_Y[3], Acce_P[3], Acce_Yinli[3], Acce_All[3];
	double Body_Direction[3], X_Direction[3], Y_Direction[3];

	// ��ɺ�һ�������������׶�		
	S_hengjie = 3.2;
	S_penkou = 3.1;
	// �����������ٶ�
	double X_Dinxin0, Y_Dinxin0, Z_Dinxin0;	//���Ĵ��ֱ������ϵ�³�ʼλ������
	toECEF(120, 35, 1000, X_Dinxin0, Y_Dinxin0, Z_Dinxin0);
	double Longitude, Latitude, Height, M_per_s;
	ECEFtoWGS84(Longitude, Latitude, Height, (X_Dinxin0+yc[2]), (Y_Dinxin0+yc[3]),( Z_Dinxin0+yc[4]));

	double X_Dinxin, Y_Dinxin, Z_Dinxin;


	toECEF(Latitude, Longitude, Height, X_Dinxin, Y_Dinxin, Z_Dinxin);
	double C_Dinxin[] = { X_Dinxin, Y_Dinxin, Z_Dinxin };	// ����λ��
	double C_Dixin[] = { X_Dinxin, Y_Dinxin, Z_Dinxin };
	double r;
	r = sqrt(pow(C_Launch[0] + C_Dixin[0], 2) + pow(C_Launch[1] + C_Dixin[1], 2) + pow(C_Launch[2] + C_Dixin[2], 2));    // ����ʸ����
	X_Dinxin = X_Dinxin0 + C_Launch[0];
	Y_Dinxin = Y_Dinxin0 + C_Launch[1];
	Z_Dinxin = Z_Dinxin0 + C_Launch[2];
	xyz_to_latlon(X_Dinxin, Y_Dinxin, Z_Dinxin, Latitude, Longitude);
	Height = sqrt(X_Dinxin * X_Dinxin + Y_Dinxin * Y_Dinxin + Z_Dinxin * Z_Dinxin) - calculateDistanceToCenter(Latitude);
	// Height = sqrt(pow(C_Launch[0] + C_Dinxin[0], 2) + pow(C_Launch[1] + C_Dinxin[1], 2) + pow(C_Launch[2] + C_Dinxin[2], 2)) - r0;
	v = sqrt(V_Launch[0] * V_Launch[0] + V_Launch[1] * V_Launch[1] + V_Launch[2] * V_Launch[2]);
	g = -fM / (r * r);               // �������ٶȴ�С
	for (size_t i = 0; i < 3; i++) {
		Acce_Yinli[i] = g * (C_Launch[i] + C_Dinxin[i]) / r;
	}

	// ������ϼ��ٶ�
	crossProduct(V_Launch, A_zizhuan, cross_Keshi);
	double Acce_Keshi[3];
	for (size_t i = 0; i < 3; ++i) {
		Acce_Keshi[i] = 2 * cross_Keshi[i];
		if (v < 1) Acce_Keshi[i] = 0;
	}

	// ����ǣ�����ٶ�
	crossProduct(A_zizhuan, C_Dinxin, cross_QianLian);
	crossProduct(A_zizhuan, cross_QianLian, Acce_QianLian);
	for (size_t i = 0; i < 3; ++i) {
		Acce_QianLian[i] = -cross_QianLian[i];
	}

	// ���㷢�����������ٶ�

	// ���ݶα�ѡ����������������

	// ���ݵ�ǰλ�ø߶ȼ������ѹǿp�������ܶ�ru������a
	daqi_US(Height, &p, &ru, &a);
	if (Time_all < 180) {
		P0 = edcz(n1, Time_all, P0_1_Time, P0_1) * 1000;
		P = P0 + S_penkou * (1 - p / p0);
		// ������Сת�������ٶ�
		A_P = P / Weight;
		// ����ϵ����������
		// 
		M_per_s = edcz((sizeof(P0_1_Time) / sizeof(P0_1_Time[0])), Time_all, P0_1_Time, M_per_s_1);
		fai_now = edcz(sizeof(fai_time) / sizeof(fai_time[0]), Time_fly, fai_time, fai);
		if (Time_all > 68) fai_now = 35;
		fai_now = fai_now * M_PI / 180.0;
		Body_Direction[0] = cos(fai_now);
		Body_Direction[1] = 0;
		Body_Direction[2] = sin(fai_now);

		for (size_t i = 0; i < 3; ++i) {
			Acce_P[i] = A_P * Body_Direction[i];
		}
		/*Acce_P[0] = F_P * cos(fai_now);
		Acce_P[1] = F_P * sin(fai_now);
		Acce_P[2] = 0;*/
	}
	else double Acce_P[] = {0,0,0};

	if (Height < 80000) {
		// �����������
		// F_Cx*q*S_hengjie		// q=1/2 * ru *  v^2	// ������������С
		Ma_Now = v / a;	// �������С
		Cx = edcz(sizeof(Cx_1) / sizeof(Cx_1[0]), Ma_Now, Ma_1, Cx_1);	// 
		Cy = Cy_1;	// 

		X = Cx * 0.5 * ru * pow(v, 2) * S_hengjie;
		Y = Cy * 0.5 * ru * pow(v, 2) * S_hengjie;
		// ת���ٶ�

		for (size_t i = 0; i < 3; ++i) {
			X_Direction[i] = -V_Launch[i] / (v * sqrt(3));
		}
		//// y_direction
		Y_Direction[0] = -V_Launch[1] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
		Y_Direction[1] = V_Launch[0] * V_Launch[2] / (pow(2, 0.5) * sqrt(pow(-V_Launch[1] * V_Launch[2], 2) + pow(V_Launch[0] * V_Launch[2], 2)));
		Y_Direction[1] = 0;
		for (size_t i = 0; i < 3; ++i) {
			Acce_X[i] = X * X_Direction[i];
			Acce_Y[i] = Y * Y_Direction[i];
			//Acce_X[i] = 0;
			//Acce_Y[i] = 0;
		}
	}
	else {
		for (size_t i = 0; i < 3; ++i) { Acce_X[i] = 0; Acce_Y[i] = 0; }
		M_per_s	= 0;
	}

	// �ϼ��ٶ�
	for (size_t i = 0; i < 3; ++i) {
		Acce_All[i] = Acce_Yinli[i] + Acce_Keshi[i] + Acce_QianLian[i] + Acce_P[i] + Acce_X[i] + Acce_Y[i];
		V_Launch[i] += Acce_All[i] * h;
		C_Launch[i] += V_Launch[i] * h;
	}

	dy[0] = h;
	dy[1] = M_per_s;
	dy[2] = yc[5];
	dy[3] = yc[6];
	dy[4] = yc[7];
	dy[5] = Acce_All[0];
	dy[6] = Acce_All[1];
	dy[7] = Acce_All[2];

	return;
}

//���������
void runk(int n, double h, double y[], double dy[], double yc[], double y1[])
{
	// n-Ԫ�ظ�����h-������y[]-����ֵ��dy[]-΢�֣�yc[]��y1[]-������
	//	extern void rf(int n,double y[25],double dy[25],double tj);
	int i, k;
	double a[5];
	a[1] = h / 2.0;  a[2] = a[1];  a[3] = h; a[4] = h;
	for (i = 0; i < n; i++)  y1[i] = y[i];	// ������һ�������ĺ���ֵ

	for (k = 1; k <= 3; k++)
	{
		for (i = 0; i < n; i++)
		{
			yc[i] = y1[i] + a[k] * dy[i];
			y[i] = y[i] + a[k + 1] * dy[i] / 3.0;
		}

		rfun(n, yc, dy);
	}
	for (i = 0; i < n; i++)  y[i] = y[i] + a[1] * dy[i] / 3.0;
	rfun(n, y, dy);
	return;
}
