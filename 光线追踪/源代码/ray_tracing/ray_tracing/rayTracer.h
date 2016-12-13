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
	Color rayTracing(vector3 ray_O, vector3 ray_V, double weight, double cur_n);//光线追踪，cur_n为当前折射率
	Color cal_local_light(object *obj, Light *light, vector3 ray_V, double cur_n);//计算局部光照

	void initObject(vector<object*> &obj);//初始化场景中的物体
	void initLight(vector<Light*> &light_array);//初始化场景中的光源
	Color depthOfField(vector3 ray_O, vector3 ray_V, vector3 focal, double f, double R, int depthOfField_time);//计算景深
};



#endif