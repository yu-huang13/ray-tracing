#include "camera.h"
#include <cstdlib> //rand()

Camera::Camera(){//设置参数
	H = 760 / 1;
	W = 1280 / 1;
	O = vector3(0, -1.8, 0.5);
	N = vector3(0, 1, 0);

	if (N == vector3(0, 0, 1) || N == vector3(0, 0, -1)){
		X = vector3(0, 1, 0);
		Y = vector3(1, 0, 0);
	}
	else {
		X = N * vector3(0, 0, 1);
		Y = N * X;
		X = X.unitVector();
		Y = Y.unitVector();
	}

	len_X = 1 * 1.25;
	len_Y = 0.59 * 1.25;
	dx = len_X / W;
	dy = len_Y / H;
	data = new Color*[H];
	for (int i = 0; i < H; ++i)
		data[i] = new Color[W];
}

Camera::~Camera()
{
	for (int i = 0; i < H; ++i)
		delete[]data[i];
	delete[]data;
}

//在当前像素点范围内随机发射光线
vector3 Camera::emit(int i, int j){
	vector3 DX = X * ((rand() % 100) / 100.0 * dx);
	vector3 DY = Y * ((rand() % 100) / 100.0 * dy);
	return (X * len_X * (2 * j  / (double)W - 1) + DX + Y * len_Y * (2 * i / (double)H - 1) + DY + N).unitVector();
}

//输出到.bmp文件
void Camera::output(string filename){
	Mat image(H, W, CV_8UC3);
	
	for (int i = 0; i < image.rows; ++i)
	for (int j = 0; j < image.cols; ++j)
		image.at<Vec3b>(i, j) = Vec3b(data[i][j].b * 255, data[i][j].g * 255, data[i][j].r * 255);

	imwrite(filename, image);
	imshow("picture", image);
	waitKey(0);
}

