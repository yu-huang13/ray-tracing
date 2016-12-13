#ifndef COLOR_H
#define COLOR_H

#include <iostream>
using namespace std;

struct Color
{
	double r, g, b;
	Color(double rr = 0, double gg = 0, double bb = 0);
	Color(const Color& c);

	friend Color operator + (const Color &A, const Color &B);
	friend Color operator - (const Color &A, const Color &B);
	friend Color operator * (const Color &A, const double &num);
	friend Color operator / (const Color &A, const double &num);
	friend Color operator * (const Color &A, const Color &B);//µã³Ë

	friend Color& operator += (Color &A, const Color &B);
	friend Color& operator -= (Color &A, const Color &B);
	friend Color& operator *= (Color &A, const double &num);
	friend Color& operator /= (Color &A, const double &num);
	friend Color& operator *= (Color &A, const Color &B);//µã³Ë

	friend ostream& operator << (ostream &output, Color &A);

	void limit();
};

#endif