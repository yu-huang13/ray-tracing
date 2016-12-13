#include <iostream>

#include "scene.h"
#include "float.h"


Scene::Scene(){
	camera = new Camera;
}

Scene::~Scene() {
	object *temp;
	while (obj.size()){
		temp = obj.back();
		obj.pop_back();
		delete temp;
		temp = NULL;
	}
	Light *light_temp;
	while (light_array.size()){
		light_temp = light_array.back();
		light_array.pop_back();
		delete light_temp;
		light_temp = NULL;
	}
	delete camera;
}

void Scene::add_object(object *pobj){
	obj.push_back(pobj);
}
void Scene::add_light(Light *l){
	light_array.push_back(l);
}
void Scene::set_environmentLight(Color c){
	evironmentLight = c;
}
object* Scene::getFirstIntersect(vector3 ray_O, vector3 ray_V){
	double dist = DBL_MAX;
	object* nearest = NULL;
	int size = obj.size();
	for (int i = 0; i < size; ++i){
		if (obj[i]->intersect_test(ray_O, ray_V)){
			if (obj[i]->intersect->dist < dist){
				dist = obj[i]->intersect->dist;
				nearest = obj[i];
			}
		}
	}
	size = light_array.size();
	for (int i = 0; i < size; ++i){
		if (light_array[i]->intersect_test(ray_O, ray_V)){
			if (light_array[i]->intersect->dist < dist){
				dist = light_array[i]->intersect->dist;
				nearest = light_array[i];
			}
		}
	}
	return nearest;
}







