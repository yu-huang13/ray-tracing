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
	vector3 emit(int i, int j);//����й�㵽i��j�����ص�ĵ�λ�����������ڵ�ǰ���ص㷶Χ������������
	void output(string filename);//��data��ɫ���������.bmp�ļ�

	vector3 O, N, X, Y;//�й�㣬�������λ���򣬾�ͷ��λ���򣬾�ͷ����λ����,Ĭ�ϸй�㵽��ͷ����Ϊ1
	Color **data;//�洢ÿ�����ص����ɫ
	double len_X, len_Y;//��ͷ������
	double dx, dy;//ÿС��ĳ����� dx = len_X / W; dy = len_Y / H
	int W, H;//bmpͼ������ؿ����س�

};



#endif