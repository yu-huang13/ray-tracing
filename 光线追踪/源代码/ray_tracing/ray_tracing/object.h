#ifndef OBJECT_H
#define OBJECT_H

#include "vector3.h"
#include "color.h"
#include <cstdlib>
#include <vector>

using namespace std;

//材质
struct Material
{
	Material(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p);
	Material();
	Material(const Material &m);
	~Material() {}
	void setMaterial(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p);
	Color color, absor;//物体颜色，折射时的吸收参数
	double diff, spec;//漫反射所占百分比，镜面反射所占百分比（phone模型）
	double refl, refr;//反射所占百分比，折射所占百分比(光线追踪)
	double n;//折射率
};

//交点信息
struct Intersect
{
	Intersect() {}
	Intersect(const Intersect &I) : N(I.N), pos(I.pos), dist(I.dist) {}
	~Intersect() {}
	vector3 N, pos;//碰撞点的法向(与发射光线处于同一面)，碰撞点的坐标
	double dist;//碰撞前光线走过的距离
};

//物体类基类
class object
{
public:
	object(Material m);
	object(const object &obj);
	~object();
	
	//存储纹理所用的图片
	void setStripData(string name);

	//计算是否相交，如相交则交点信息存储于intersect中
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V) = 0; //发射点坐标，方向（已单位化）
	//是否为光源
	virtual bool islight() = 0;
	//返回一定位置的颜色
	virtual Color getColor(vector3 pos) = 0;

	Intersect *intersect;
	Material *material;

	//纹理
	bool has_strip;
	Color **data;
	int strip_W, strip_H;//像素宽、像素高
	double strip_len_X, strip_len_Y;
};

//球
class Ball: public object
{
public:
	Ball(Material m, vector3 OO, double RR);
	~Ball() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
	vector3 X, Y;//用于纹理贴图

private:
	vector3 O;//球心坐标
	double R;//半径
};

//无限大平面
class Plane : public object
{
public:
	Plane(Material m, vector3 OO, vector3 NN);
	~Plane() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
private:
	vector3 O, N, X, Y;//平面上一点，平面法向，平面上的x、y轴单位向量
};

//有限大平面
class Plane_XY : public object
{
public:
	Plane_XY(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY);
	~Plane_XY() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
private:
	vector3 O, N, X, Y; //平面上一点，平面法向，平面上的x、y轴单位向量
	double len_X, len_Y;//平面长、宽

};

//三角面片
class Triangle
{
public:
	Triangle(vector3& A, vector3& B, vector3& C);
	Triangle(const Triangle& T);
	~Triangle() {}
	bool intersect_test(vector3 ray_O, vector3 ray_V, Intersect* intersect);

	vector3 center;//x、y、z最大值
	vector3 p[3];//三个点的坐标
	vector3 N;//法向

	double det(vector3 v1, vector3 v2, vector3 v3);
};

//包围盒
class BoundingBox
{
public:
	vector3 minP, maxP;//xyz坐标均为最小的点，xyz坐标均为最大的点

	BoundingBox() {}
	BoundingBox(vector3 minPP, vector3 maxPP): minP(minPP), maxP(maxPP) {}
	bool intersect_test(vector3 ray_O, vector3 ray_V);

private:
	inline double cal_dis(vector3& O, vector3& N, vector3& ray_O, vector3& ray_V);//计算光线起点到面（O, N）的距离
	inline double max3(double& a, double& b, double& c);//返回a、b、c三个数的最大值
	bool inBox(vector3& pos);//判断pos点是否在盒子表面以及盒子内

};

//KD树节点
struct KDnode
{
	BoundingBox box;
	KDnode *lc, *rc;
	int tri;
};

//比较器
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

//由.obj文件读入的物体
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
	vector<Triangle> triangles;//下标从0开始
	long long gotoleaf_num;

	double& Map(vector3& v, int& flag);//v.x---0, v.y---1, v.z---2
	bool readFile(char* filename);//读入Obj文件
	void rotate(double ax, double ay, double az);//绕i轴旋转ai度（弧度）
	void coordinate_trans(double ax, double ay, double az, double rate, vector3 dir);//绕i轴旋转ai度（弧度），放大倍数，移动方向
	KDnode* buildKDTree(int* Array, int length, int deepth, BoundingBox box, Compare& compare);//建立KD树
	bool KDSearch(KDnode* root, vector3& ray_O, vector3& ray_V);//搜索KD树
};

//光源基类
class Light: public object
{
public:
	Light(Material m) : object(m) {}
	~Light() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V) { return false; }
	virtual bool islight();
	virtual Color getColor(vector3 pos);
	vector<vector3> point_array;//从存储每个点光源的位置
	Color everyPointLight;//存储每个点光源的颜色
};

//点光源
class pointLight : public Light
{
public:
	pointLight(Material m, vector3 OO);
	~pointLight() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);
};

//面光源
class planeLight : public Light//矩形面光源
{
public:
	planeLight(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY, int W, int H);//NN：法向量
	~planeLight() {}
	virtual bool intersect_test(vector3 ray_O, vector3 ray_V);
	virtual bool islight();
	virtual Color getColor(vector3 pos);

	vector3 O, N, X, Y;//矩形中心，矩形法向，长单位向量，矩形高单位向量
	int W, H;//将矩形分成 W * H 的网格（点光源）
	double len_X, len_Y;//矩形长，高
};

#endif

