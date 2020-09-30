

#ifndef MYCOMMON_H_
#define MYCOMMON_H_
#include <CImg.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <limits.h>

using namespace std;
using namespace cimg_library;

double beta = 2;

class Point
{
public:
  Point() {}
  Point(int _col, int _row) : row(_row), col(_col) {}
  int row, col;
};


double multivariateGaussian(CImg<double> img, CImg<double> mean, CImg<double> variance)
{
	CImg<double> temp = img-mean;
	CImg<double> temp2 = temp * variance.invert(); //* temp.transpose();
	temp2 = temp2 * temp.transpose();
	double val = temp2(0,0);
	val = exp(-0.5 * val);
	val = (1 / (pow(2.0*M_PI, 3.0/2.0) * sqrt(variance.det()))) * val;
	return val;
}

bool isIn(vector<Point> points, Point point)
{
	for(int i = 0; i < points.size(); i++)
	{
		if(point.col == points[i].col && point.row == points[i].row)
			return true;
	}
	return false;
}
#endif /* MYCOMMON_H_ */
