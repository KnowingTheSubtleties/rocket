//#include <iostream>  
//#include <cmath>  
//
//const double a = 6378245.0; // 克拉索夫斯基椭球的长半轴  
//const double f = 1.0 / 298.3; // 扁率  
//const double b = a * (1 - f);  // 短半轴  
//const double e2 = (2 - f) * f;  // 第一偏心率的平方  
//const double M_PI = 3.14159265358979323846;
//// 转换度到弧度  
//double deg2rad(double deg) {
//    return deg * M_PI / 180.0;
//}
//
//// 地心大地直角坐标系坐标  
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
//// 计算地心直角坐标系下坐标到地表的高度（实际为验证函数）  
//double cartesian_to_surface_height(double X, double Y, double Z) {
//    double lambda_rad = atan2(Y, X);
//    double B_rad = atan2(Z, sqrt(X * X + Y * Y));
//    double N = a / sqrt(1 - e2 * sin(B_rad) * sin(B_rad));
//
//    // 忽略Z轴中由发射高度引起的部分，仅计算基于椭球体的高度  
//    double h = (sqrt(X * X + Y * Y) / cos(B_rad)) - N;
//    return h;
//}
//
//int main3() {
//    double lambda_T = 120.0;
//    double B_T = 35.0;
//    double A_T = 45.0; // 注意：A_T在此计算中未直接使用  
//    double Height = 1000.0;
//    double X, Y, Z;
//    geodetic_to_cartesian(lambda_T, B_T, Height, X, Y, Z);
//
//    std::cout << "地心直角坐标系坐标: (" << X << ", " << Y << ", " << Z << ")" << std::endl;
//    std::cout << "地表高度（忽略发射高度）: " << cartesian_to_surface_height(X, Y, Z) << " 米" << std::endl;
//
//    system("pause");
//    return 0;
//}