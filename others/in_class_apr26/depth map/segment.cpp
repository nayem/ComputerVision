// Skeleton code for B657 A4 Part 2.
// D. Crandall
//
//
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <CImg.h>
#include <assert.h>
#include "utill.h"
#include "loopy_belief.h"
using namespace cimg_library;
using namespace std;


CImg<double> naive_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
	cout << "beta = " << beta << endl;
  // implement this in step 1...
  //  this placeholder just returns a random disparity map

  CImg<double> mean(3,1);
  mean.fill(0);
  CImg<double> variance(3,3);
  variance.fill(0);

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

  double D[2][img.width() * img.height()];

  int pixelCount = 0;
  for(int i=0; i<img.height(); i++)
    for(int j=0; j<img.width(); j++)
    {
    	for(int l = 0; l < 2; l++)
    	{
			if(isIn(fg, Point(j,i)))
			{
				if(l)
				{
					D[l][pixelCount] = 0;
				}
				else
					D[l][pixelCount] = 9999999999;
			}
			else
				if(isIn(bg, Point(j,i)))
				{
					if(l)
					{
						D[l][pixelCount] = 9999999999;
					}
					else
					{
						D[l][pixelCount] = 0;
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
						D[l][pixelCount] = -log(multivariateGaussian(pixel, mean, variance));
					}
					else
					{
						D[l][pixelCount] = beta;
					}
				}
    	}
    	pixelCount++;
    }


  CImg<double> result(img.width(), img.height());
  pixelCount = 0;
  //cout <<"working" << endl;
  for(int i=0; i<img.height(); i++)
    for(int j=0; j<img.width(); j++)
    {
      double l0 = D[0][pixelCount];
      double l1 = D[1][pixelCount];
      if(l0 < l1)
      {
    	  result(j, i, 0, 0) = 0;
      }
      else
      {
    	  result(j, i, 0, 0) = 1;
      }

      pixelCount++;
    }

  //cout << "not working" << endl;
  return result;
}

CImg<double> mrf_segment(const CImg<double> &img, const vector<Point> &fg, const vector<Point> &bg)
{
  return getMRFSegmentation(img, fg, bg);
}

// Take in an input image and a binary segmentation map. Use the segmentation map to split the 
//  input image into foreground and background portions, and then save each one as a separate image.
//
// You'll just need to modify this to additionally output a disparity map.
//
void output_segmentation(const CImg<double> &img, const CImg<double> &labels, const string &fname)
{
  // sanity checks. If one of these asserts fails, you've given this function invalid arguments!
  assert(img.height() == labels.height());
  assert(img.width() == labels.width());

  CImg<double> img_fg = img, img_bg = img;

  for(int i=0; i<labels.height(); i++)
    for(int j=0; j<labels.width(); j++)
      {
	if(labels(j,i) == 0)
	  img_fg(j,i,0,0) = img_fg(j,i,0,1) = img_fg(j,i,0,2) = 0;
	else if(labels(j,i) == 1)
	  img_bg(j,i,0,0) = img_bg(j,i,0,1) = img_bg(j,i,0,2) = 0;
	else
	  assert(0);
      }

  img_fg.get_normalize(0,255).save((fname + "_fg.png").c_str());
  img_bg.get_normalize(0,255).save((fname + "_bg.png").c_str());
}

int main(int argc, char *argv[])
{
    
  if(argc < 3)
    {
      cerr << "usage: " << argv[0] << " image_file seeds_file beta" << endl;
      return 1;
    }
    

  string input_filename1 = argv[1], input_filename2 = argv[2];
    if(argc == 4)
        beta = atof(argv[3]);
  //BP_ITERATIONS = atoi(argv[4]);
  // read in images and gt
  CImg<double> image_rgb(input_filename1.c_str());
  CImg<double> seeds_rgb(input_filename2.c_str());

  // figure out seed points 
  vector<Point> fg_pixels, bg_pixels;
  for(int i=0; i<seeds_rgb.height(); i++)
    for(int j=0; j<seeds_rgb.width(); j++)
      {
	// blue --> foreground
	if(max(seeds_rgb(j, i, 0, 0), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 2) > 100)
	  fg_pixels.push_back(Point(j, i));

	// red --> background
	if(max(seeds_rgb(j, i, 0, 2), seeds_rgb(j, i, 0, 1)) < 100 && seeds_rgb(j, i, 0, 0) > 100)
	  bg_pixels.push_back(Point(j, i));
      }

  // do naive segmentation
  CImg<double> labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
  output_segmentation(image_rgb, labels, input_filename1 + "-naive_segment_result");

  // do mrf segmentation
  labels = mrf_segment(image_rgb, fg_pixels, bg_pixels);
  output_segmentation(image_rgb, labels, input_filename1 + "-mrf_segment_result");

  return 0;
}
