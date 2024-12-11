//#include <iostream>  
//#include <cmath>  
//
//const double a = 6378245.0; // ��������˹������ĳ�����  
//const double f = 1.0 / 298.3; // ����  
//const double b = a * (1 - f);  // �̰���  
//const double e2 = (2 - f) * f;  // ��һƫ���ʵ�ƽ��  
//const double M_PI = 3.14159265358979323846;
//// ת���ȵ�����  
//double deg2rad(double deg) {
//    return deg * M_PI / 180.0;
//}
//
//// ���Ĵ��ֱ������ϵ����  
//double geodetic_to_cartesian(double lambda, double B, double height, double& X, double& Y, double& Z) {
//    double lambda_rad = deg2rad(lambda);
//    double B_rad = deg2rad(B);
//    double N = a / sqrt(1 - e2 * sin(B_rad) * sin(B_rad));
//
//    X = (N + height) * cos(B_rad) * cos(lambda_rad);
//    Y = (N + height) * cos(B_rad) * sin(lambda_rad);
//    Z = ((b * b) / (a * a) * N + height) * sin(B_rad);
//
//    return 0;
//}
//
//// �������ֱ������ϵ�����굽�ر�ĸ߶ȣ�ʵ��Ϊ��֤������  
//double cartesian_to_surface_height(double X, double Y, double Z) {
//    double lambda_rad = atan2(Y, X);
//    double B_rad = atan2(Z, sqrt(X * X + Y * Y));
//    double N = a / sqrt(1 - e2 * sin(B_rad) * sin(B_rad));
//
//    // ����Z�����ɷ���߶�����Ĳ��֣����������������ĸ߶�  
//    double h = (sqrt(X * X + Y * Y) / cos(B_rad)) - N;
//    return h;
//}
//
//int main3() {
//    double lambda_T = 120.0;
//    double B_T = 35.0;
//    double A_T = 45.0; // ע�⣺A_T�ڴ˼�����δֱ��ʹ��  
//    double Height = 1000.0;
//    double X, Y, Z;
//    geodetic_to_cartesian(lambda_T, B_T, Height, X, Y, Z);
//
//    std::cout << "����ֱ������ϵ����: (" << X << ", " << Y << ", " << Z << ")" << std::endl;
//    std::cout << "�ر�߶ȣ����Է���߶ȣ�: " << cartesian_to_surface_height(X, Y, Z) << " ��" << std::endl;
//
//    system("pause");
//    return 0;
//}