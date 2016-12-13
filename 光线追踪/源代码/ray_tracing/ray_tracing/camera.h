#ifndef CAMERA_H
#define CAMERA_H
#include <string>

#include "vector3.h"
#include "color.h"
#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;
class Camera
{
public:
	Camera();
	~Camera();
	vector3 emit(int i, int j);//计算感光点到i行j列像素点的单位方向向量，在当前像素点范围内随机发射光线
	void output(string filename);//将data颜色矩阵输出到.bmp文件

	vector3 O, N, X, Y;//感光点，照相机单位法向，镜头宽单位法向，镜头长单位法向,默认感光点到镜头距离为1
	Color **data;//存储每个像素点的颜色
	double len_X, len_Y;//镜头长、宽
	double dx, dy;//每小格的长、宽： dx = len_X / W; dy = len_Y / H
	int W, H;//bmp图像的像素宽，像素长

};



#endif