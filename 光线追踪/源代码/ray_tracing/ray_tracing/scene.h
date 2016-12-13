#ifndef SCENE_H
#define SCENE_H

#include "object.h"
#include "camera.h"
#include <vector>
#include <fstream>
using namespace std;
class Scene
{
public:
	Color evironmentLight;
	Camera *camera;
	vector <object*> obj;
	vector <Light*> light_array;

	Scene();
	~Scene();
	void add_object(object *pobj);//�����м�������
	void add_light(Light *l);//�����м����Դ
	void set_environmentLight(Color c);//���û�����
	object* getFirstIntersect(vector3 ray_O, vector3 ray_V);//���������壬����ֵ����Ϊ��

};


#endif