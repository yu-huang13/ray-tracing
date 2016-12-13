#include <iostream>

#include "rayTracer.h"

using namespace std;

int main()
{
	rayTracer *rt = new rayTracer;
	rt->initScene();
	rt->cal_bmp("picture.bmp");
	delete rt;
	return 0;
}

