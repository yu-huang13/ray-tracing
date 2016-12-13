#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include "rayTracer.h"

const double EPS = 1e-6;
const double PI = 3.1415926535897932;

rayTracer::rayTracer(){
	scene = new Scene;
}
rayTracer::~rayTracer(){
	delete scene;
}

void rayTracer::initObject(vector<object*> &obj){
	Material m;
	object *obj_temp;

	//������
	m.setMaterial(Color(1, 1, 1), Color(0, 0, 0), 0.05, 0.1, 0, 0.85, 1.7);
	obj_temp = new Ball(m, vector3(-1.5, 4.5, -1.5), 0.5);
	obj.push_back(obj_temp);

	//���淴����
	m.setMaterial(Color(0.3, 0.3, 1), Color(0, 0, 0), 0, 0, 1, 0, 1);
	obj_temp = new Ball(m, vector3(0, 4.5, -1.4), 0.6);
	obj.push_back(obj_temp);

	//����ʯ��
	m.setMaterial(Color(1, 1, 1), Color(0, 0, 0), 0.4, 0.2, 0.3, 0, 1);
	obj_temp = new Ball(m, vector3(1.6, 4.5, -1.3), 0.7);
	obj_temp->setStripData("ball.bmp");//ball1.bmp
	obj.push_back(obj_temp);
	 
	//�ذ� 
	m.setMaterial(Color(1, 0.5, 0.3), Color(0, 0, 0), 0.6, 0, 0.4, 0, 1);
	obj_temp = new Plane(m, vector3(0, 0, -2), vector3(0, 0, 1));
	obj_temp->setStripData("floor.bmp");//floor1.bmp
	obj.push_back(obj_temp);
	
	//ǰ����ǽ
	m.setMaterial(Color(1, 1, 1), Color(0, 0, 0), 0.2, 0, 0.8, 0, 1);
	obj_temp = new Plane(m, vector3(0, 16, 0), vector3(0, -1, 0));
	obj.push_back(obj_temp);

	//�ҷ���ǽ
	m.setMaterial(Color(1, 0.5, 0.3), Color(0, 0, 0), 0.8, 0.2, 0, 0, 1);
	obj_temp = new Plane(m, vector3(10, 0, 0), vector3(-1, 0, 0));
	obj_temp->setStripData("wall.bmp");
	obj.push_back(obj_temp);

	//�󷽵�ǽ
	m.setMaterial(Color(1, 0.5, 0.3), Color(0, 0, 0), 0.8, 0, 0.2, 0, 1);
	obj_temp = new Plane(m, vector3(-10, 0, 0), vector3(1, 0, 0));
	obj_temp->setStripData("wall.bmp");
	obj.push_back(obj_temp);

	//�Ϸ���ǽ
	m.setMaterial(Color(1, 0.5, 0.3), Color(0, 0, 0), 0.8, 0.2, 0, 0, 1);
	obj_temp = new Plane(m, vector3(0, 0, 11), vector3(0, 0, -1));
	obj_temp->setStripData("wall.bmp");
	obj.push_back(obj_temp);

	//�󷽵�ǽ
	m.setMaterial(Color(1, 0.5, 0.3), Color(0, 0, 0), 0.8, 0.2, 0, 0, 1);
	obj_temp = new Plane(m, vector3(0, -6, 0), vector3(0, 1, 0));
	obj_temp->setStripData("wall.bmp");
	obj.push_back(obj_temp);
	
	//bunny
	m.setMaterial(Color(0, 1, 0.25), Color(0, 0, 0), 0.4, 0.2, 0.4, 0, 1);
	obj_temp = new OBJ(m, "bunny.obj", PI / 2, 0, 0, 15, vector3(5, 7, -2.6));
	obj.push_back(obj_temp);

	//dragon
	m.setMaterial(Color(1, 0, 0), Color(0, 0, 0), 0.4, 0.2, 0.4, 0, 1);
	obj_temp = new OBJ(m, "dragon.obj", PI / 2, 0, PI, 1.5, vector3(-5, 7, -1));
	obj.push_back(obj_temp);

	//dinosaur
	m.setMaterial(Color(1, 0.5, 0.5), Color(0, 0, 0), 0.4, 0.2, 0.4, 0, 1);
	obj_temp = new OBJ(m, "dinosaur.obj", 0, 0, PI / 2, 0.026, vector3(-4, 10, -0.75));
	obj.push_back(obj_temp);

	//kitten
	m.setMaterial(Color(1, 1, 1), Color(0, 0, 0), 0, 0.1, 0, 0.9, 1.7);
	obj_temp = new OBJ(m, "kitten.obj", PI / 2, 0, 0, 0.045, vector3(0, 8.5, -2));
	obj.push_back(obj_temp);

	//horse
	m.setMaterial(Color(0.5, 0.25, 0), Color(0, 0, 0), 0.4, 0.2, 0.4, 0, 1);
	obj_temp = new OBJ(m, "horse.obj", 0, 0, PI / 2, 1.5, vector3(4, 10, -0.9));
	obj.push_back(obj_temp);

}

void rayTracer::initLight(vector<Light*> &light_array){
	Material m;
	pointLight *plight;
	m.setMaterial(Color(1.4, 1.4, 1.4), Color(0, 0, 0), 0, 0, 0, 0, 1);

	//���Դ
	plight = new pointLight(m, vector3(0, 2, 8));//-3 3 3
	light_array.push_back(plight);

	//���Դ������Ӱ��
	/*planeLight *planel = new planeLight(m, vector3(0, 2, 8), vector3(0, 0, -1), 1, 1, 7, 7);
	light_array.push_back(planel);*/
	
}

void rayTracer::initScene(){
	srand((unsigned)time(NULL));
	initObject(scene->obj);
	initLight(scene->light_array);
	scene->set_environmentLight(Color(0.1, 0.1, 0.1));
}

//ray_V: �ӵ㵽��ײ��ĵ�λ����
Color rayTracer::cal_local_light(object *obj, Light *light, vector3 ray_V, double cur_n){
	Color color(0, 0, 0);

	int size = light->point_array.size();
	for (int i = 0; i < size; ++i){

		vector3 L = (light->point_array[i] - obj->intersect->pos).unitVector();

		double LN = L.dot(obj->intersect->N);

		if ((obj->islight() == false &&  cur_n != 1) || LN < EPS)//���ߴ��������ڲ� �� ��Դ�����䵽�õ�
			return color;

		object* temp_obj = scene->getFirstIntersect(obj->intersect->pos, L);

		if (temp_obj != NULL && temp_obj->islight() == false) {//�������ڵ���Դ
			if ((temp_obj->intersect->pos - obj->intersect->pos).module2() < (light->point_array[i] - obj->intersect->pos).module2())
				return color;
		}

		color += light->everyPointLight * LN * obj->material->diff * obj->getColor(obj->intersect->pos);//������

		double RV = (L.reflect(obj->intersect->N)).dot(-ray_V);
		if (RV > EPS)
			color += light->everyPointLight * pow(RV, 50) * obj->material->spec * obj->getColor(obj->intersect->pos);//���淴��
	}
	
	color.limit();
	return color;
}

Color rayTracer::rayTracing(vector3 ray_O, vector3 ray_V, double weight, double cur_n)
{
	Color color(0, 0, 0);
	weight *= 0.95;
	if (weight < 0.01) return color; //˥��
	
	object *obj = scene->getFirstIntersect(ray_O, ray_V);
	
	if (obj == NULL) return color; //û�н���
	else if (obj->islight())
		return obj->getColor(vector3(0, 0, 0));

	//�ֲ�����
	if (cur_n == 1){//��������ײ�ż���ֲ�����
		int size = scene->light_array.size();
		for (int i = 0; i < size; ++i) {
			color += cal_local_light(obj, scene->light_array[i], ray_V, cur_n);
		}	
	}

	//׷�ٷ������
	vector3 Refl_ray = (-ray_V).reflect(obj->intersect->N);
	double Wrefl = obj->material->refl;
	color += rayTracing(obj->intersect->pos, Refl_ray, Wrefl, cur_n) * obj->getColor(obj->intersect->pos)  * Wrefl;


	//׷���������
	vector3 Refr_ray(0, 0, 0);
	double n;//���������
	if (cur_n == 1){
		n = obj->material->n; cur_n = obj->material->n;
	}
	else {
		n = 1 / obj->material->n; cur_n = 1;
	}
	if ((-ray_V).refract(obj->intersect->N, n, Refr_ray)){//������߿��㣬���ǡ�ȫ���䡱
		double Wrefr = obj->material->refr;
		color += rayTracing(obj->intersect->pos, Refr_ray, Wrefr, cur_n) * obj->getColor(obj->intersect->pos) * Wrefr;
	}

	return color * weight;
}

vector3 randomPosInCircle(vector3 O, double R, vector3 X, vector3 Y){//X, YΪX��Y�᷽��ĵ�λ����
	double H, l;//�����������ĽǶȣ���Բ�ĵľ���
	H = (rand() % 10000) * 2 * PI / 10000;
	l = (rand() % 10000) * R / 10000;
	return O + X * l * cos(H) + Y * l * sin(H);
}

Color rayTracer::depthOfField(vector3 ray_O, vector3 ray_V, vector3 focal, double f, double R, int depthOfField_time){
	vector3 O = ray_O;
	focal = O + ray_V * f;
	Color result(0, 0, 0);
	for (int u = 0; u < depthOfField_time; ++u){
		ray_O = randomPosInCircle(O, R, scene->camera->X, scene->camera->Y);
		ray_V = (focal - ray_O).unitVector();
		result += rayTracing(ray_O, ray_V, 1, 1);
	}
	return result / depthOfField_time;
}

 void rayTracer::cal_bmp(string filename){
	vector3 ray_V;
	vector3 ray_O;
	Color color;
	Color AC(0, 0, 0);

	//double f = 6.3, R = 0.2;//���࣬��Ȧ�뾶                �翪�������򽫸��� ȡ��ע��
	//vector3 focal;//����                                   �翪�������򽫸��� ȡ��ע��
	//int depthOfField_time = 50;//��������ߵĴ���         �翪�������򽫸��� ȡ��ע��

	int MTC_time = 1;//����ݷ������ߵ���Ŀ��Ϊ1��ص������    �翪����������Ϊ1�����رտ����

	for (int i = 0; i < scene->camera->H; ++i){
		cout << i << endl;
		for (int j = 0; j < scene->camera->W; ++j){
			AC = Color(0, 0, 0);

			for (int k = 0; k < MTC_time; ++k){
				ray_V = scene->camera->emit(i, j);
				ray_O = scene->camera->O;

				//focal = ray_O + ray_V * f;                                            �翪�������򽫸��� ȡ��ע��
				//AC += depthOfField(ray_O, ray_V, focal, f, R, depthOfField_time);     �翪�������򽫸��� ȡ��ע��
				
				AC += rayTracing(ray_O, ray_V, 1, 1);//�翪�������򽫸��� ע��
			}
			color = scene->evironmentLight + AC / MTC_time;
			color.limit();
			scene->camera->data[i][j] = color;
		}
	}
	scene->camera->output(filename);
}



