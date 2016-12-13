#ifndef RAYTRACER_H
#define RAYTRACER_H

#include "scene.h"

class rayTracer
{
public:
	rayTracer();
	~rayTracer();
	Scene *scene;

	void initScene();
	void cal_bmp(string filename);

private:
	Color rayTracing(vector3 ray_O, vector3 ray_V, double weight, double cur_n);//����׷�٣�cur_nΪ��ǰ������
	Color cal_local_light(object *obj, Light *light, vector3 ray_V, double cur_n);//����ֲ�����

	void initObject(vector<object*> &obj);//��ʼ�������е�����
	void initLight(vector<Light*> &light_array);//��ʼ�������еĹ�Դ
	Color depthOfField(vector3 ray_O, vector3 ray_V, vector3 focal, double f, double R, int depthOfField_time);//���㾰��
};



#endif