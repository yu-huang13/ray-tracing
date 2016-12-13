#include "color.h"

Color::Color(double rr, double gg, double bb) : r(rr), g(gg), b(bb) {}
Color::Color(const Color& c){
	r = c.r; g = c.g; b = c.b;
}
Color operator + (const Color &A, const Color &B){
	return Color(A.r + B.r, A.g + B.g, A.b + B.b);
}
Color operator - (const Color &A, const Color &B){
	return Color(A.r - B.r, A.g - B.g, A.b - B.b);
}
Color operator * (const Color &A, const double &num){
	return Color(A.r * num, A.g * num, A.b * num);
}
Color operator / (const Color &A, const double &num){
	return Color(A.r / num, A.g / num, A.b / num);
}
Color operator * (const Color &A, const Color &B){
	return Color(A.r * B.r, A.g * B.g, A.b * B.b);
}

Color& operator += (Color &A, const Color &B){
	A = A + B;
	return A;
}
Color& operator -= (Color &A, const Color &B){
	A = A - B;
	return A;
}
Color& operator *= (Color &A, const double &num){
	A = A * num;
	return A;
}
Color& operator /= (Color &A, const double &num){
	A = A / num;
	return A;
}
Color& operator *= (Color &A, const Color &B){
	A = A * B;
	return A;
}

void Color::limit(){
	if (r > 1) r = 1;
	if (g > 1) g = 1;
	if (b > 1) b = 1;
}

ostream& operator << (ostream &output, Color &A){
	output << "(" << A.r << ", " << A.g << ", " << A.b << ")" << endl;
	return output;
}
