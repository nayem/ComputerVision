/*
 * loopy_belief.h
 *
 *  Created on: Apr 27, 2017
 *
 *      inspired from http://nghiaho.com/?page_id=1366
 */

#include <iostream>
#include <vector>
#include <string>
#include "CImg.h"
#include <limits>
#include "utill.h"
using namespace cimg_library;
using namespace std;

enum DIRECTION {LEFT, RIGHT, UP, DOWN, DATA};

typedef double TYPE;

// parameters, specific to dataset
int BP_ITERATIONS = 1;
const int LABELS = 2;
const int LAMBDA = 20;
const int SMOOTHNESS_TRUNC = 2;

struct Pixel
{
    // Each pixel has 5 'message box' to store incoming data
    TYPE msg[5][LABELS];
    int best_assignment;
};

struct MRF2D
{
    std::vector <Pixel> grid;
    int width, height;
};

// Application specific code
void InitDataCost(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg, MRF2D &mrf);
CImg <double> DataCost(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg);
TYPE SmoothnessCost(int i, int j);

// Loppy belief propagation specific
void BP(MRF2D &mrf, DIRECTION direction);
void SendMsg(MRF2D &mrf, int x, int y, DIRECTION direction);
TYPE MAP(MRF2D &mrf);

CImg<double> getMRFSegmentation(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
    MRF2D mrf;

    InitDataCost(img, fg, bg, mrf);

    for(int i=0; i < BP_ITERATIONS; i++) {
        BP(mrf, RIGHT);
        BP(mrf, LEFT);
        BP(mrf, UP);
        BP(mrf, DOWN);

        TYPE energy = MAP(mrf);

        cout << "iteration " << (i+1) << "/" << BP_ITERATIONS << ", energy = " << energy << endl;
    }

    //cv::Mat output = cv::Mat::zeros(mrf.height, mrf.width, CV_8U);
    CImg<double> output(mrf.width, mrf.height);
    output.fill(0);

    for(int y=LABELS; y < mrf.height-LABELS; y++)
    {
        for(int x=LABELS; x < mrf.width-LABELS; x++)
        {
            // Increase the intensity so we can see it
            output(x,y,0,0) = mrf.grid[y*mrf.width+x].best_assignment;// * (256/LABELS);
            if((x == 247 && y == 198) || (x == 246 && y == 198))
			{
            	//cout << "x = " << x << " y = " << y << " label = " <<output(x,y,0,0)<< " mrf.data= " << mrf.grid[y*img.width()+x].msg[DATA][0]<< "," <<mrf.grid[y*img.width()+x].msg[DATA][1]<< endl;
            	//cout << "R = " << img(x,y,0,0) << " G="<< img(x,y,0,1) << " B = " << img(x,y,0,0)<<endl;
            	//cout<< " mrf.data= " << mrf.grid[y*img.width()+x].msg[DATA][0]<<endl;
			}
        }
    }

    return output;
}

CImg<double> DataCost(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
	cout << "beta = " << beta << endl;
	  CImg<double> mean(3,1);
	  mean.fill(0);
	  CImg<double> variance(3,3);
	  variance.fill(0);

		//cout << "1" << endl;
	  for(int i = 0; i < fg.size(); i++)
	  {
		  int col = fg[i].col;
		  int row = fg[i].row;
		  mean(0,0) += img(col, row, 0, 0);
		  mean(1,0) += img(col, row, 0, 1);
		  mean(2,0) += img(col, row, 0, 2);
	  }

	  mean(0,0)/=fg.size();
	  mean(1,0)/=fg.size();
	  mean(2,0)/=fg.size();

	  for(int i = 0; i < fg.size(); i++)
	  {
		  int col = fg[i].col;
		  int row = fg[i].row;
		  variance(0,0) += (img(col, row, 0, 0) - mean(0,0)) * (img(col, row, 0, 0) - mean(0,0));
		  variance(1,1) += (img(col, row, 0, 1) - mean(1,0)) * (img(col, row, 0, 1) - mean(1,0));
		  variance(2,2) += (img(col, row, 0, 2) - mean(2,0)) * (img(col, row, 0, 2) - mean(2,0));
	  }

	  variance(0,0)/=(fg.size()-1);
	  variance(1,1)/=(fg.size()-1);
	  variance(2,2)/=(fg.size()-1);

	  variance(0,0)=sqrt(variance(0,0));
	  variance(1,1)=sqrt(variance(1,1));
	  variance(2,2)=sqrt(variance(2,2));

	  cout << "variance(0,0) = " << variance(0,0) << " variance(1,1) = " << variance(1,1) << " variance(2,2) = " << variance(2,2)<<endl;
	  //cout << "2" << endl;

	  //double D[2][img.width() * img.height()];
	  CImg<double> D(img.width(), img.height(), 1, LABELS);

	  //int pixelCount = 0;
	  for(int i=0; i<img.height(); i++)
	    for(int j=0; j<img.width(); j++)
	    {
	    	for(int l = 0; l < LABELS; l++)
	    	{
				if(isIn(fg, Point(j,i)))
				{
					if(l)
					{
						D(j,i,0,l) = 0;
					}
					else
						D(j,i,0,l) = 9999999999;
				}
				else
					if(isIn(bg, Point(j,i)))
					{
						if(l)
						{
							D(j,i,0,l) = 9999999999;
						}
						else
						{
							D(j,i,0,l) = 0;
						}
					}
					else
					{
						if(l)
						{
							CImg<double> pixel(3,1);
							pixel(0,0) = img(j,i,0,0);
							pixel(1,0) = img(j,i,0,1);
							pixel(2,0) = img(j,i,0,2);
							D(j,i,0,l) = -log(multivariateGaussian(pixel, mean, variance));
						}
						else
						{
							D(j,i,0,l) = beta;
						}
					}
	    	}
	    }
    return D;
}

TYPE SmoothnessCost(int i, int j)
{
    //int d = i - j;

    //return LAMBDA*std::min(abs(d), SMOOTHNESS_TRUNC);
	//cout << "beta" << beta << endl;
	//return (i != j)*(beta/5);
	if(j ==1 &&i == 0)
		return beta*5;
	return (i != j)*(beta);
}

void InitDataCost(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg, MRF2D &mrf)
{
    // Cache the data cost results so we don't have to recompute it every time

    mrf.width = img.width();
    mrf.height = img.height();

    int total = mrf.width*mrf.height;

    mrf.grid.resize(total);

    // Initialise all messages to zero
    for(int i=0; i < total; i++) {
        for(int j=0; j < 5; j++) {
            for(int k=0; k < LABELS; k++) {
                mrf.grid[i].msg[j][k] = 0;
            }
        }
    }

    // Add a border around the image
    int border = LABELS;

    CImg<double> D = DataCost(img, fg, bg);

    for(int y=border; y < mrf.height-border; y++) {
        for(int x=border; x < mrf.width-border; x++) {
            for(int i=0; i < LABELS; i++) {
                mrf.grid[y*img.width()+x].msg[DATA][i] = D(x, y, 0, i);
            }
        }
    }
}

void SendMsg(MRF2D &mrf, int x, int y, DIRECTION direction)
{
    TYPE new_msg[LABELS];

    int width = mrf.width;

    for(int i=0; i < LABELS; i++)
    {
        TYPE min_val = INT_MAX;

        for(int j=0; j < LABELS; j++)
        {
            TYPE p = 0;

            p += SmoothnessCost(i,j);
            //cout << p;
            p += mrf.grid[y*width+x].msg[DATA][j];

            // Exclude the incoming message direction that we are sending to
            if(direction != LEFT) p += mrf.grid[y*width+x].msg[LEFT][j];
            if(direction != RIGHT) p += mrf.grid[y*width+x].msg[RIGHT][j];
            if(direction != UP) p += mrf.grid[y*width+x].msg[UP][j];
            if(direction != DOWN) p += mrf.grid[y*width+x].msg[DOWN][j];

            min_val = std::min(min_val, p);
        }

        new_msg[i] = min_val;
    }

    for(int i=0; i < LABELS; i++)
    {
        switch(direction)
        {
            case LEFT:
            mrf.grid[y*width + x-1].msg[RIGHT][i] = new_msg[i];
            break;

            case RIGHT:
            mrf.grid[y*width + x+1].msg[LEFT][i] = new_msg[i];
            break;

            case UP:
            mrf.grid[(y-1)*width + x].msg[DOWN][i] = new_msg[i];
            break;

            case DOWN:
            mrf.grid[(y+1)*width + x].msg[UP][i] = new_msg[i];
            break;

            default:
            assert(0);
            break;
        }
    }
}

void BP(MRF2D &mrf, DIRECTION direction)
{
    int width = mrf.width;
    int height = mrf.height;

    switch(direction)
    {
        case RIGHT:
        	//cout << "right" << endl;
        for(int y=0; y < height; y++)
        {
            for(int x=0; x < width-1; x++)
            {
                SendMsg(mrf, x, y, direction);
            }
        }
        break;

        case LEFT:
        	//cout << "left" << endl;
        for(int y=0; y < height; y++)
        {
            for(int x=width-1; x >= 1; x--)
            {
                SendMsg(mrf, x, y, direction);
            }
        }
        break;

        case DOWN:
        	//cout << "down" << endl;
        for(int x=0; x < width; x++)
        {
            for(int y=0; y < height-1; y++)
            {
                SendMsg(mrf, x, y, direction);
            }
        }
        break;

        case UP:
        	//cout << "up" << endl;
        for(int x=0; x < width; x++)
        {
            for(int y=height-1; y >= 1; y--)
            {
                SendMsg(mrf, x, y, direction);
            }
        }
        break;

        case DATA:
        assert(0);
        break;
    }
}

TYPE MAP(MRF2D &mrf)
{
    // Finds the MAP assignment as well as calculating the energy

    // MAP assignment
    for(size_t i=0; i < mrf.grid.size(); i++)
    {
        TYPE best = std::numeric_limits<TYPE>::max();
        for(int j=0; j < LABELS; j++) {
            TYPE cost = 0;

            cost += mrf.grid[i].msg[LEFT][j];
            cost += mrf.grid[i].msg[RIGHT][j];
            cost += mrf.grid[i].msg[UP][j];
            cost += mrf.grid[i].msg[DOWN][j];
            cost += mrf.grid[i].msg[DATA][j];

            if(cost < best)
            {
                best = cost;
                mrf.grid[i].best_assignment = j;
            }
        }
    }

    int width = mrf.width;
    int height = mrf.height;

    // Energy
    TYPE energy = 0;

    for(int y=0; y < mrf.height; y++) {
        for(int x=0; x < mrf.width; x++) {
            int cur_label = mrf.grid[y*width+x].best_assignment;

            // Data cost
            energy += mrf.grid[y*width+x].msg[DATA][cur_label];

            if(x-1 >= 0)     energy += SmoothnessCost(cur_label, mrf.grid[y*width+x-1].best_assignment);
            if(x+1 < width)  energy += SmoothnessCost(cur_label, mrf.grid[y*width+x+1].best_assignment);
            if(y-1 >= 0)     energy += SmoothnessCost(cur_label, mrf.grid[(y-1)*width+x].best_assignment);
            if(y+1 < height) energy += SmoothnessCost(cur_label, mrf.grid[(y+1)*width+x].best_assignment);
        }
    }

    return energy;
}
