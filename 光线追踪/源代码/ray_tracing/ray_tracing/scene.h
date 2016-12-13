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
	void add_object(object *pobj);//场景中加入物体
	void add_light(Light *l);//场景中加入光源
	void set_environmentLight(Color c);//设置环境光
	object* getFirstIntersect(vector3 ray_O, vector3 ray_V);//检测最近物体，返回值可能为空

};


#endif