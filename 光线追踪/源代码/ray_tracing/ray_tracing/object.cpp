#include "object.h"
#include "opencv2/opencv.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>//rand()
#include <sstream>
using namespace std;
using namespace cv;

const double EPS = 1e-6;
const double PI = 3.1415926535897932;

Material::Material(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p)
: color(c), absor(ab), diff(diff_p), spec(spec_p), refl(refl_p), refr(refr_p), n(n_p) {}

Material::Material(const Material &m){
	color = m.color; absor = m.absor; diff = m.diff; spec = m.spec; refl = m.refl; refr = m.refr; n = m.n;
}

Material::Material(){
	color = Color(0, 0, 0); absor =Color(0, 0, 0); diff = 0; spec = 0; refl = 0; refr = 0; n = 1;
}

void Material::setMaterial(Color c, Color ab, double diff_p, double spec_p, double refl_p, double refr_p, double n_p){
	color = c; absor = ab; diff = diff_p; spec = spec_p; refl = refl_p; refr = refr_p; n = n_p;
}

object::object(Material m){
	material = new Material(m);
	intersect = new Intersect;
	has_strip = false;
}

object::object(const object &obj){
	material = new Material(*(obj.material));
	intersect = new Intersect(*(obj.intersect));
	has_strip = false;
}

object::~object(){
	delete material;
	delete intersect;
	if (has_strip){
		for (int i = 0; i < strip_H; ++i)
			delete[]data[i];
		delete[]data;
	}
}

void object::setStripData(string filename){
	Mat image = imread(filename);
	strip_W = image.cols;
	strip_H = image.rows;
	strip_len_Y = strip_H * 0.01;
	strip_len_X = strip_W * 0.01;
	data = new Color*[strip_H];
	
	Vec3b temp;
	for (int i = 0; i < image.rows; ++i){
		data[i] = new Color[strip_W];
		for (int j = 0; j < image.cols; ++j){
			temp = image.at<Vec3b>(i, j);
			data[i][j] = Color(temp[2] / 255.0, temp[1] / 255.0, temp[0] / 255.0);
		}
	}
	has_strip = true;
}

Ball::Ball(Material m, vector3 OO, double RR) : object(m) {
	O = OO;
	R = RR;
	X = vector3(0, 0, 1);
	Y = vector3(1, 0, 0);
}

bool Ball::intersect_test(vector3 ray_O, vector3 ray_V)
{
	vector3 rayO_O = O - ray_O;
	double R2 = R * R;
	double module2_rayO_O = rayO_O.module2();
	double projection = ray_V.dot(rayO_O);
	double half_chord2 = R * R - module2_rayO_O + projection * projection;

	if (half_chord2 < EPS)//光线所在直线与球无交
		return false;
	if (module2_rayO_O < R * R - EPS){//光线起点在球内
		intersect->dist = sqrt(half_chord2) + projection;
		intersect->pos = ray_V * intersect->dist + ray_O;
		intersect->N = (O - intersect->pos).unitVector();
	}
	else{//光线起点在球外或球面上
		if (projection < EPS) 
			return false;
		double temp = projection - sqrt(half_chord2);
		if (temp > EPS){//在球外
			intersect->dist = temp;
			intersect->pos = ray_V * intersect->dist + ray_O;
			intersect->N = (intersect->pos - O).unitVector();
		}
		else{//在球面上
			intersect->dist = projection + sqrt(half_chord2);
			intersect->pos = ray_V * intersect->dist + ray_O;
			intersect->N = (O - intersect->pos).unitVector();
		}
	}
	return true;
}

bool Ball::islight(){
	return false;
}

Color Ball::getColor(vector3 pos){
	if (!has_strip)
		return material->color;
	else{
		double x = (pos - O).unitVector().dot(X);
		double y = (pos - O).unitVector().dot(Y);
		
		int i = (int)(acos(y) / PI * (strip_H - 1));
		int j = (int)(acos(x) / PI * (strip_W - 1));

		return data[i][j];
	}
}

Plane::Plane(Material m, vector3 OO, vector3 NN) : object(m)
{
	O = OO;
	N = NN;
	if (N == vector3(0, 0, 1) || N == vector3(0, 0, -1)){
		X = vector3(0, 1, 0);
		Y = vector3(1, 0, 0);
	}
	else {
		X = N * vector3(0, 0, 1);
		Y = N * X;
		X = X.unitVector();
		Y = Y.unitVector();
	}
}

bool Plane::intersect_test(vector3 ray_O, vector3 ray_V)
{
	double cosH = ray_V.dot(N);
	if (fabs(cosH) < EPS) return false;//平行
	intersect->dist = N.dot(N * (O.dot(N)) - ray_O) / cosH;
	if (intersect->dist < EPS) return false;
	intersect->pos = ray_V * intersect->dist + ray_O;
	intersect->N = cosH < EPS ? N : -N;
	return true;
}

bool Plane::islight(){
	return false;
}

Color Plane::getColor(vector3 pos){
	if (!has_strip)
		return material->color;
	else{
		double x = (pos - O).dot(X);
		double y = (pos - O).dot(Y);

		x = x / strip_len_X;
		y = y / strip_len_Y;

		x = x - (int)x;
		y = y - (int)y;

		if (x < 0) x += 1;
		if (y < 0) y += 1;

		int i = (int)(y * (strip_H - 1));
		int j = (int)(x * (strip_W - 1));

		return data[i][j];
	}
}

Plane_XY::Plane_XY(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY) : object(m)
{
	O = OO;
	N = NN;
	if (N == vector3(0, 0, 1) || N == vector3(0, 0, -1)){
		X = vector3(0, 1, 0);
		Y = vector3(1, 0, 0);
	}
	else {
		X = N * vector3(0, 0, 1);
		Y = N * X;
		X = X.unitVector();
		Y = Y.unitVector();
	}
	len_X = len_XX;	len_Y = len_YY;
}

bool Plane_XY::intersect_test(vector3 ray_O, vector3 ray_V)
{
	double cosH = ray_V.dot(N);
	if (fabs(cosH) < EPS) return false;//平行
	intersect->dist = N.dot(N * (O.dot(N)) - ray_O) / cosH;
	if (intersect->dist < EPS) return false;
	intersect->pos = ray_V * intersect->dist + ray_O;
	vector3 O_pos = intersect->pos - O;
	if (fabs(O_pos.dot(X)) > len_X || fabs(O_pos.dot(Y)) > len_Y)
		return false;
	intersect->N = cosH < EPS ? N : -N;
	return true;
}

bool Plane_XY::islight(){
	return false;
}

Color Plane_XY::getColor(vector3 pos){
	if (!has_strip)
		return material->color;
	else{
		double x = (pos - O).dot(X);
		double y = (pos - O).dot(Y);

		x = x / strip_len_X;
		y = y / strip_len_Y;

		x = x - (int)x;
		y = y - (int)y;

		if (x < 0) x += 1;
		if (y < 0) y += 1;

		int i = (int)(y * (strip_H - 1));
		int j = (int)(x * (strip_W - 1));
		
		return data[i][j];
	}
}

Triangle::Triangle(vector3& A, vector3& B, vector3& C){
	p[0] = A; p[1] = B; p[2] = C;
}

Triangle::Triangle(const Triangle& T){
	p[0] = T.p[0]; p[1] = T.p[1]; p[2] = T.p[2];
	center = T.center;
	N = T.N;
}

//计算行列式|v1 v2 v3|
double Triangle::det(vector3 v1, vector3 v2, vector3 v3){
	return v1.x * (v2.y * v3.z - v3.y * v2.z) - v2.x * (v1.y * v3.z - v3.y * v1.z) + v3.x * (v1.y * v2.z - v2.y * v1.z);
}

bool Triangle::intersect_test(vector3 ray_O, vector3 ray_V, Intersect* intersect){
	vector3 E1 = p[0] - p[1], E2 = p[0] - p[2], S = p[0] - ray_O;
	double cramer = det(ray_V, E1, E2);
	if (fabs(cramer) < EPS) return false;//方程组无解（线面平行）
	double t = det(S, E1, E2) / cramer, B = det(ray_V, S, E2) / cramer, r = det(ray_V, E1, S) / cramer;
	double a = 1 - B - r;
	
	if (t > EPS && a >= 0 && a <= 1 && B >= 0 && B <= 1 && r >= 0 && r <= 1){
		intersect->dist = t;
		intersect->pos = ray_O + ray_V * t;
		double sign = ray_V.dot(N);
		intersect->N = sign < 0 ? N : -N;
		return true;
	}
	else//反向或交点不在三角形内
		return false;
}


inline double BoundingBox::cal_dis(vector3& O, vector3& N, vector3& ray_O, vector3& ray_V){
	return N.dot(N * (O.dot(N)) - ray_O) / (ray_V.dot(N));
}

inline double BoundingBox::max3(double& a, double& b, double& c){
	if (a > b){
		return a > c ? a : c;
	}
	else{//a <= b
		return b > c ? b : c;
	}
}

bool BoundingBox::inBox(vector3& pos){
	return pos.x >= (minP.x - EPS) && pos.x <= (maxP.x + EPS)
		&& pos.y >= (minP.y - EPS) && pos.y <= (maxP.y + EPS)
		&& pos.z >= (minP.z - EPS) && pos.z <= (maxP.z + EPS);
}

bool BoundingBox::intersect_test(vector3 ray_O, vector3 ray_V){//Woo算法求交

	if (inBox(ray_O)) return true;//光线的起点在长方体内

	double dis_x1, dis_x2, dis_y1, dis_y2, dis_z1, dis_z2, cosH;
	vector3 N;
	//处理与x轴垂直的平面
	N = vector3(1, 0, 0);
	dis_x1 = cal_dis(minP, N, ray_O, ray_V);
	dis_x2 = cal_dis(maxP, N, ray_O, ray_V);
	if (dis_x1 > dis_x2) swap(dis_x1, dis_x2);

	//处理与y轴垂直的平面
	N = vector3(0, 1, 0);
	dis_y1 = cal_dis(minP, N, ray_O, ray_V);
	dis_y2 = cal_dis(maxP, N, ray_O, ray_V);
	if (dis_y1 > dis_y2) swap(dis_y1, dis_y2);

	//处理与z轴垂直的平面
	N = vector3(0, 0, 1);
	cosH = ray_V.dot(N);
	dis_z1 = cal_dis(minP, N, ray_O, ray_V);
	dis_z2 = cal_dis(maxP, N, ray_O, ray_V);
	if (dis_z1 > dis_z2) swap(dis_z1, dis_z2);

	double maxt = max3(dis_x1, dis_y1, dis_z1);
	if (maxt < 0) return false;//在光线的反向
	vector3 pos = ray_O + ray_V * maxt;
	return inBox(pos);
}

void new_matrix(double**& a, int n, int m){//n行m列
	a = new double*[n];
	for (int i = 0; i < n; ++i)
		a[i] = new double[m];
}

void init0_matrix(double**& A, int n, int m){
	for (int i = 0; i < n; ++i)
	for (int j = 0; j < m; ++j)
		A[i][j] = 0;
}

void delete_matrix(double**& a, int n, int m){
	for (int i = n; i < n; ++i)
		delete[]a[i];
	delete[]a;
}

void matrix_multi(double** A, int nA, int mAnB, double** B, int mB, double** result){//n行m列
	for (int i = 0; i < nA; ++i)
	for (int j = 0; j < mB; ++j){
		result[i][j] = 0;
		for (int k = 0; k < mAnB; ++k){
			result[i][j] += A[i][k] * B[k][j];
		}
	}
}

void print_matrix(double** A, int n, int m){
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < m; ++j){
			cout << A[i][j] << " ";
		}
		cout << endl;
	}
}

//旋转
void OBJ::rotate(double ax, double ay, double az){
	double** A, **B, **C;
	new_matrix(A, 3, 3); init0_matrix(A, 3, 3);
	new_matrix(B, 3, 3); init0_matrix(B, 3, 3);
	new_matrix(C, 3, 3); init0_matrix(C, 3, 3);
	double sinax = sin(ax), cosax = cos(ax), sinay = sin(ay), cosay = cos(ay), sinaz = sin(az), cosaz = cos(az);
	A[0][0] = 1; A[1][1] = A[2][2] = cosax; A[1][2] = sinax; A[2][1] = -sinax;
	B[1][1] = 1; B[0][0] = B[2][2] = cosay; B[0][2] = -sinay; B[2][0] = sinay;
	C[2][2] = 1; C[0][0] = C[1][1] = cosaz; C[0][1] = sinaz; C[1][0] = -sinaz;

	double** temp, **rotate;
	new_matrix(temp, 3, 3); new_matrix(rotate, 3, 3);
	matrix_multi(B, 3, 3, C, 3, temp);
	matrix_multi(A, 3, 3, temp, 3, rotate);

	double** pos, **result; new_matrix(pos, 1, 3); new_matrix(result, 1, 3);
	for (vector<Triangle>::iterator it = triangles.begin(); it != triangles.end(); ++it){
		for (int k = 0; k < 3; ++k){
			pos[0][0] = (*it).p[k].x; pos[0][1] = (*it).p[k].y; pos[0][2] = (*it).p[k].z;
			matrix_multi(pos, 1, 3, rotate, 3, result);
			(*it).p[k].x = result[0][0]; (*it).p[k].y = result[0][1]; (*it).p[k].z = result[0][2];//旋转
		}
	}
	delete_matrix(A, 3, 3); delete_matrix(B, 3, 3); delete_matrix(C, 3, 3); delete_matrix(temp, 3, 3); delete_matrix(rotate, 3, 3);
	delete_matrix(pos, 1, 3); delete_matrix(result, 1, 3);
}

//记得计算包围盒
void OBJ::coordinate_trans(double ax, double ay, double az, double rate, vector3 dir){
	//旋转
	rotate(ax, ay, az);

	//缩放、平移，计算包围盒
	vector3 minP(DBL_MAX, DBL_MAX, DBL_MAX), maxP(-DBL_MAX, -DBL_MAX, -DBL_MAX);
	for (vector<Triangle>::iterator it = triangles.begin(); it != triangles.end(); ++it){
		for (int k = 0; k < 3; ++k){
			(*it).p[k] = (*it).p[k] * rate + dir;//缩放+平移
			if ((*it).p[k].x < minP.x) minP.x = (*it).p[k].x; if ((*it).p[k].x > maxP.x) maxP.x = (*it).p[k].x;
			if ((*it).p[k].y < minP.y) minP.y = (*it).p[k].y; if ((*it).p[k].y > maxP.y) maxP.y = (*it).p[k].y;
			if ((*it).p[k].z < minP.z) minP.z = (*it).p[k].z; if ((*it).p[k].z > maxP.z) maxP.z = (*it).p[k].z;
		}	
		(*it).center = ((*it).p[0] + (*it).p[1] + (*it).p[2]) / 3;
		(*it).N = (((*it).p[1] - (*it).p[0]) * ((*it).p[2] - (*it).p[0])).unitVector();
	}
	Box.maxP = maxP; Box.minP = minP;
	/*cout << "box: " << Box.minP << " " << Box.maxP << endl;
	system("pause");*/
	
}

OBJ::OBJ(Material m, char* filename, double ax, double ay, double az, double rate, vector3 dir) :object(m) {
	readFile(filename);
	coordinate_trans(ax, ay, az, rate, dir);
	int size = triangles.size();
	int* Array = new int[size];
	for (int i = 0; i < size; ++i)
		Array[i] = i;
	Compare compare(triangles);
	Root = buildKDTree(Array, size, 0, Box, compare);
	delete[]Array;
}

bool OBJ::readFile(char* filename){
	ifstream fin(filename);
	if (!fin){
		cout << "fail to open " << filename << endl;
		return false;
	}
	string line;
	string mark;
	stringstream fin2;
	vector3 input;
	int T[3];
	vector<vector3> points;
	points.push_back(vector3(0, 0, 0));//为使points的下标从1开始有效，预先加入

	while (getline(fin, line)){
		fin2.clear();
		fin2.str(line);
		fin2 >> mark;
		if (mark == "v"){
			fin2 >> input.x >> input.y >> input.z;
			points.push_back(input);
		}
		if (mark == "f"){
			for (int i = 0; i < 3; ++i)
				fin2 >> T[i];
			triangles.push_back(Triangle(points[T[0]], points[T[1]], points[T[2]]));
		}
	}
	fin.close();
	return true;
}

double& OBJ::Map(vector3& v, int& flag){
	switch (flag){
	case 0: return v.x;
	case 1: return v.y;
	case 2: return v.z;
	default: cout << "Map error" << endl; return v.x;
	}
}

//建立KD树，左右孩子节点包围盒可能重叠
KDnode* OBJ::buildKDTree(int* Array, int length, int deepth, BoundingBox box, Compare& compare){
	KDnode* root = new KDnode;
	root->box = box;
	if (length == 1){
		root->lc = NULL; root->rc = NULL;
		root->tri = Array[0];
		return root;
	}

	int flag = deepth % 3;
	compare.flag = flag;
	sort(Array, Array + length, compare);
	int s_num = length / 2;//左孩子节点中三角形的个数
	int b_num = length - s_num;//右孩子节点中三角形的个数

	//三角形子集划分
	int *smaller = new int[s_num];
	int *bigger = new int[b_num];
	memcpy(smaller, Array, s_num * sizeof(int));
	memcpy(bigger, Array + s_num, b_num * sizeof(int));

	//确定左右包围盒大小,包围盒可能重叠
	double smaller_max = -DBL_MAX, bigger_min = DBL_MAX;
	double temp;
	for (int i = 0; i < s_num; ++i){
		for (int k = 0; k < 3; ++k){
			temp = Map(triangles[smaller[i]].p[k], flag);
			if (temp > smaller_max)
				smaller_max = temp;
		}
	}
	for (int i = 0; i < b_num; ++i){
		for (int k = 0; k < 3; ++k){
			temp = Map(triangles[bigger[i]].p[k], flag);
			if (temp < bigger_min)
				bigger_min = temp;
		}
	}
	BoundingBox smallBox = box, bigBox = box;
	Map(smallBox.maxP, flag) = smaller_max; Map(bigBox.minP, flag) = bigger_min;

	//递归建树
	root->lc = buildKDTree(smaller, s_num, deepth + 1, smallBox, compare);
	root->rc = buildKDTree(bigger, b_num, deepth + 1, bigBox, compare);

	delete[]smaller;
	delete[]bigger;
	return root;
}

//判断光线与物体是否有交，有交则交点信息存在intersect中
bool OBJ::KDSearch(KDnode* root, vector3& ray_O, vector3& ray_V){
	if (root->lc == NULL && root->rc == NULL){//叶子节点
		gotoleaf_num++;
		Intersect* temp = new Intersect;
		bool inter = false;
		if (triangles[root->tri].intersect_test(ray_O, ray_V, temp)){//和三角面片相交
			inter = true;
			if (temp->dist < intersect->dist){//交点到光线起点距离更短
				intersect->dist = temp->dist; intersect->pos = temp->pos; intersect->N = temp->N;
			}
		}
		delete temp;
		return inter;
	}
	if (!root->box.intersect_test(ray_O, ray_V))//若和当前包围盒无交
		return false;

	bool l_inter = KDSearch(root->lc, ray_O, ray_V);
	bool r_inter = KDSearch(root->rc, ray_O, ray_V);

	return l_inter || r_inter;

}

bool OBJ::intersect_test(vector3 ray_O, vector3 ray_V){
	intersect->dist = DBL_MAX;
	gotoleaf_num = 0;
	return KDSearch(Root, ray_O, ray_V);
}

bool OBJ::islight(){
	return false;
}
Color OBJ::getColor(vector3 pos){
	return material->color;
}

bool Light::islight(){
	return true;
}

Color Light::getColor(vector3 pos){
	return material->color;
}

pointLight::pointLight(Material m, vector3 OO): Light(m){
	everyPointLight = m.color;
	point_array.push_back(OO);
}

bool pointLight::intersect_test(vector3 ray_O, vector3 ray_V){
	return (point_array[0] - ray_O) == ray_V;
}

bool pointLight::islight(){
	return true;
}

Color pointLight::getColor(vector3 pos){
	return material->color;
}

planeLight::planeLight(Material m, vector3 OO, vector3 NN, double len_XX, double len_YY, int WW, int HH): Light(m) {
	vector3 t = vector3(0, 0, 1) * vector3(0, 0, 1);
	O = OO;
	N = NN;
	if (N == vector3(0, 0, 1) || N == vector3(0, 0, -1)){
		X = vector3(0, 1, 0);
		Y = vector3(1, 0, 0);
	}
	else {
		X = N * vector3(0, 0, 1);
		Y = N * X;
		X = X.unitVector();
		Y = Y.unitVector();
	}
	W = WW;	H = HH;

	len_X = len_XX;	len_Y = len_YY;
	double dx = len_X / W, dy = len_Y / H;
	vector3 DX, DY;
	vector3 point;
	for (int i = 0; i < H; ++i){
		for (int j = 0; j < W; ++j){
			DX = X * ((rand() % 100) / 100.0 * dx);
			DY = Y * ((rand() % 100) / 100.0 * dy);
			point = O + X * len_X * (2 * j / (double)W - 1) + DX + Y * len_Y * (2 * i / (double)H - 1) + DY;
			point_array.push_back(point);
		}
	}
	everyPointLight = m.color / (W * H);
}

bool planeLight::intersect_test(vector3 ray_O, vector3 ray_V){
	double cosH = ray_V.dot(N);
	if (fabs(cosH) < EPS) return false;//平行
	intersect->dist = N.dot(N * (O.dot(N)) - ray_O) / cosH;
	if (intersect->dist < EPS) return false;
	intersect->pos = ray_V * intersect->dist + ray_O;
	vector3 O_pos = intersect->pos - O;
	if (fabs(O_pos.dot(X)) > len_X || fabs(O_pos.dot(Y)) > len_Y)
		return false;
	intersect->N = cosH < EPS ? N : -N;
	return true;
}

bool planeLight::islight(){
	return true;
}

Color planeLight::getColor(vector3 pos){
	return material->color;
}