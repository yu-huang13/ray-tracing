#ifndef OBJECT_H
#define OBJECT_H

#include "vector3.h"
#include "color.h"
#include <cstdlib>
#include <vector>

using namespace std;

//����
struct Material
{
	Material(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p);
	Material();
	Material(const Material &m);
	~Material() {}
	void setMaterial(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p);
	Color color, absor;//������ɫ������ʱ�����ղ���
	double diff, spec;//��������ռ�ٷֱȣ����淴����ռ�ٷֱȣ�phoneģ�ͣ�
	double refl, refr;//������ռ�ٷֱȣ�������ռ�ٷֱ�(����׷��)
	double n;//������
};

//������Ϣ
struct Intersect
{
	Intersect() {}
	Intersect(const Intersect &I) : N(I.N), pos(I.pos), dist(I.dist) {}
	~Intersect() {}
	vector3 N, pos;//��ײ��ķ���(�뷢����ߴ���ͬһ��)����ײ�������
	double dist;//��ײǰ�����߹��ľ���
};

//���������
class object
{
public:
	object(Material m);
	object(const object &obj);
	~object();
	
	//�洢�������õ�ͼƬ
	void setStripData(string name);

	//�����Ƿ��ཻ�����ཻ�򽻵���Ϣ�洢��intersect��
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V) = 0; //��������꣬�����ѵ�λ����
	//�Ƿ�Ϊ��Դ
	virtual bool islight() = 0;
	//����һ��λ�õ���ɫ
	virtual Color getColor(vector3 pos) = 0;

	Intersect *intersect;
	Material *material;

	//����
	bool has_strip;
	Color **data;
	int strip_W, strip_H;//���ؿ����ظ�
	double strip_len_X, strip_len_Y;
};

//��
class Ball: public object
{
public:
	Ball(Material m, vector3 OO, double RR);
	~Ball() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
	vector3 X, Y;//����������ͼ

private:
	vector3 O;//��������
	double R;//�뾶
};

//���޴�ƽ��
class Plane : public object
{
public:
	Plane(Material m, vector3 OO, vector3 NN);
	~Plane() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
private:
	vector3 O, N, X, Y;//ƽ����һ�㣬ƽ�淨��ƽ���ϵ�x��y�ᵥλ����
};

//���޴�ƽ��
class Plane_XY : public object
{
public:
	Plane_XY(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY);
	~Plane_XY() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
private:
	vector3 O, N, X, Y; //ƽ����һ�㣬ƽ�淨��ƽ���ϵ�x��y�ᵥλ����
	double len_X, len_Y;//ƽ�泤����

};

//������Ƭ
class Triangle
{
public:
	Triangle(vector3& A, vector3& B, vector3& C);
	Triangle(const Triangle& T);
	~Triangle() {}
	bool intersect_test(vector3 ray_O, vector3 ray_V, Intersect* intersect);

	vector3 center;//x��y��z���ֵ
	vector3 p[3];//�����������
	vector3 N;//����

	double det(vector3 v1, vector3 v2, vector3 v3);
};

//��Χ��
class BoundingBox
{
public:
	vector3 minP, maxP;//xyz�����Ϊ��С�ĵ㣬xyz�����Ϊ���ĵ�

	BoundingBox() {}
	BoundingBox(vector3 minPP, vector3 maxPP): minP(minPP), maxP(maxPP) {}
	bool intersect_test(vector3 ray_O, vector3 ray_V);

private:
	inline double cal_dis(vector3& O, vector3& N, vector3& ray_O, vector3& ray_V);//���������㵽�棨O, N���ľ���
	inline double max3(double& a, double& b, double& c);//����a��b��c�����������ֵ
	bool inBox(vector3& pos);//�ж�pos���Ƿ��ں��ӱ����Լ�������

};

//KD���ڵ�
struct KDnode
{
	BoundingBox box;
	KDnode *lc, *rc;
	int tri;
};

//�Ƚ���
struct Compare
{
	Compare(vector<Triangle>& T) : triangles(T){}
	vector<Triangle>& triangles;
	int flag;
	bool operator () (const int& a, const int& b)const{
		switch (flag){
		case 0: return triangles[a].center.x < triangles[b].center.x;
		case 1: return triangles[a].center.y < triangles[b].center.y;
		case 2: return triangles[a].center.z < triangles[b].center.z;
		default: cout << "compare error!" << endl; return false;
		}
	}
};

//��.obj�ļ����������
class OBJ : public object 
{
public: 
	OBJ(Material m, char* filename, double ax, double ay, double az, double rate, vector3 dir);
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
private:
	KDnode* Root;
	BoundingBox Box;
	vector<Triangle> triangles;//�±��0��ʼ
	long long gotoleaf_num;

	double& Map(vector3& v, int& flag);//v.x---0, v.y---1, v.z---2
	bool readFile(char* filename);//����Obj�ļ�
	void rotate(double ax, double ay, double az);//��i����תai�ȣ����ȣ�
	void coordinate_trans(double ax, double ay, double az, double rate, vector3 dir);//��i����תai�ȣ����ȣ����Ŵ������ƶ�����
	KDnode* buildKDTree(int* Array, int length, int deepth, BoundingBox box, Compare& compare);//����KD��
	bool KDSearch(KDnode* root, vector3& ray_O, vector3& ray_V);//����KD��
};

//��Դ����
class Light: public object
{
public:
	Light(Material m) : object(m) {}
	~Light() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V) { return false; }
	virtual bool islight();
	virtual Color getColor(vector3 pos);
	vector<vector3> point_array;//�Ӵ洢ÿ�����Դ��λ��
	Color everyPointLight;//�洢ÿ�����Դ����ɫ
};

//���Դ
class pointLight : public Light
{
public:
	pointLight(Material m, vector3 OO);
	~pointLight() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
};

//���Դ
class planeLight : public Light//�������Դ
{
public:
	planeLight(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY, int W, int H);//NN��������
	~planeLight() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);

	vector3 O, N, X, Y;//�������ģ����η��򣬳���λ���������θߵ�λ����
	int W, H;//�����ηֳ� W * H �����񣨵��Դ��
	double len_X, len_Y;//���γ�����
};

#endif

