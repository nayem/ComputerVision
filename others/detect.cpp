//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

#include <SImage.h>
#include <SImageIO.h>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;
SDoublePlane convolve_general(const SDoublePlane &, const SDoublePlane &);
SDoublePlane transpose(const SDoublePlane &);
// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent squared gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

// Below is a helper functions that overlays rectangles
// on an image plane for visualization purpose. 

SDoublePlane generateGaussian(){
  double kernal[3][3] ={{ 0.1070, 0.1131 ,0.1070},
			{0.1131, 0.1196, 0.1131},
			{0.1070, 0.1131, 0.1070}};
  SDoublePlane output(3,3);
  for(int i=0;i<3;i++)
   for(int j=0;j<3;j++)
     output[i][j] = kernal[i][j];
  return output;
}

double getAngle(double degree){
  int angle;
  if((degree>=22.5 && degree<67.5) || (degree<-112.5 && degree>-157.5)) 
    angle = 45;
  if((degree>=67.5 && degree<112.5) || (degree<-67.5 && degree>-112.5))
    angle = 90;
  if((degree>=112.5 && degree<157.5) || (degree<-22.5 && degree>-67.5))
    angle = 135;
  if((degree>=-22.5 && degree<22.5) || (degree<-157.5 || degree>157.5))
    angle = 0;
  return angle;
}

SDoublePlane transpose(const SDoublePlane &input){
  SDoublePlane output(input.cols(),input.rows());
  for(int i=0;i<input.rows();i++){
    for(int j=0;j<input.cols();j++){
      output[j][i] = input [i][j];
    }
  }
  return output;
}

SDoublePlane matrixMultiply(const SDoublePlane &mat1,const SDoublePlane &mat2){
  SDoublePlane output(mat1.rows(),mat2.cols());
  for (int i = 0; i < mat1.rows(); i++) {
    for (int j = 0; j < mat2.cols(); j++) {
      double sum = 0;
      for (int k = 0; k < mat2.rows(); k++) {
	sum +=  mat1[i][k]*mat2[k][j];
      }
      output[i][j] = sum;
    }
  }
  return output;
}


int indexFind(int size,int index){
  if(index<0)
    return -index;
  if(index>=size)
    return size-(index-size)-1;
  return index;
}

// Draws a rectangle on an image plane, using the specified gray level value and line width.
//
void overlay_rectangle(SDoublePlane &input, int _top, int _left, int _bottom, int _right, double graylevel, int width)
{
  for(int w=-width/2; w<=width/2; w++) {
    int top = _top+w, left = _left+w, right=_right+w, bottom=_bottom+w;

    // if any of the coordinates are out-of-bounds, truncate them 
    top = min( max( top, 0 ), input.rows()-1);
    bottom = min( max( bottom, 0 ), input.rows()-1);
    left = min( max( left, 0 ), input.cols()-1);
    right = min( max( right, 0 ), input.cols()-1);
      
    // draw top and bottom lines
    for(int j=left; j<=right; j++)
	  input[top][j] = input[bottom][j] = graylevel;
    // draw left and right lines
    for(int i=top; i<=bottom; i++)
	  input[i][left] = input[i][right] = graylevel;
  }
}

// DetectedBox class may be helpful!
//  Feel free to modify.
//
class DetectedBox {
public:
  int row, col, width, height;
  double confidence;
};

// Function that outputs the ascii detection output file
void  write_detection_txt(const string &filename, const vector<DetectedBox> &cars)
{
  ofstream ofs(filename.c_str());

  for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
    ofs << s->row << " " << s->col << " " << s->width << " " << s->height << " " << s->confidence << endl;
}

// Function that outputs a visualization of detected boxes
void  write_detection_image(const string &filename, const vector<DetectedBox> &cars, const SDoublePlane &input)
{
  SDoublePlane output_planes[3];

  for(int p=0; p<3; p++)
    {
      output_planes[p] = input;
      for(vector<DetectedBox>::const_iterator s=cars.begin(); s != cars.end(); ++s)
	overlay_rectangle(output_planes[p], s->row, s->col, s->row+s->height-1, s->col+s->width-1, p==2?255:0, 2);
    }

  SImageIO::write_png_file(filename.c_str(), output_planes[0], output_planes[1], output_planes[2]);
}



// The rest of these functions are incomplete. These are just suggestions to 
// get you started -- feel free to add extra functions, change function
// parameters, etc.

// Convolve an image with a separable convolution kernel
//
SDoublePlane convolve_separable(const SDoublePlane &input, const SDoublePlane &row_filter, const SDoublePlane &col_filter)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane temp(input.rows(), input.cols());
  SDoublePlane tRow = transpose(row_filter);
  // Convolution code here
   
  SDoublePlane res =  matrixMultiply(tRow,col_filter);
  int rlen = floor(tRow.rows()/2);
  for(int r = 0; r < input.rows(); r++){
    for(int c = 0; c < input.cols(); c++){
      double sum = 0.0;
      for(int kr = -(rlen); kr <= (rlen); kr++){
	int iR = indexFind(input.rows(), r + kr);
	sum += (input[iR][c]*tRow[kr+rlen][0]);
      }
      temp[r][c] = sum;
      //cout<<sum<<" ";
    }
  }

  int clen = floor(col_filter.cols()/2);
  for(int r = 0; r < input.rows(); r++){
    for(int c = 0; c < input.cols(); c++){
      double sum = 0.0;
      for(int kc = -(clen); kc <= clen; kc++){
	int iC = indexFind(input.cols(), c + kc);
	sum += (temp[r][iC]*col_filter[0][kc+clen]);
	//cout<<temp[r][iC]<<" "<<col_filter[0][kc+clen]<<" "<<sum<<endl;
      }
      output[r][c] = sum;
    }
  }
 
  
  //output = convolve_general(input,res);
  
  return output;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
  SDoublePlane output(input.rows(), input.cols());
  // Convolution code here
  int rlen = floor(filter.rows()/2);
  int clen = floor(filter.cols()/2);
  for(int r=0; r < input.rows(); r++){
    for(int c = 0; c < input.cols(); c++){
      double sum = 0;
      for(int kr = -(rlen); kr <= rlen; kr++){
	for(int kc = -(clen) ; kc <= clen;kc++){
	  int iR = indexFind(input.rows(),r+kr);
	  int iC = indexFind(input.cols(),c+kc);
	  sum += (input[iR][iC]*filter[kr+rlen][kc+clen]);
	  //cout<<input[iR][iC]<<" "<<filter[kr+1][kc+1]<<" "<<sum<<endl;
	}
      }
      output[r][c] = sum;
    }
  }
  
  return output;
}


// Apply a sobel operator to an image, returns the result
// 
SDoublePlane sobel_gradient_filter(const SDoublePlane &input, bool _gx)
{
  SDoublePlane output(input.rows(), input.cols());
  // Implement a sobel gradient estimation filter with 1-d filters
  const double filter1[] = {1,2,1};
  const double filter2[] = {-1,0,1};
  SDoublePlane  sFilter1(1,3,filter1);
  SDoublePlane  sFilter2(1,3,filter2);
  if(_gx)
    output = convolve_separable(input,sFilter1,sFilter2);
  else
    output = convolve_separable(input,sFilter2,sFilter1);
  return output;
}

// Apply an edge detector to an image, returns the binary edge map
// 
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
  SDoublePlane output(input.rows(), input.cols());
  SDoublePlane magnitude(input.rows(), input.cols());
  SDoublePlane direction(input.rows(), input.cols());
  // Implement an edge detector of your choice, e.g.
  // use your sobel gradient operator to compute the gradient magnitude and threshold
  SDoublePlane outputX = sobel_gradient_filter(input,false);
  SDoublePlane outputY = sobel_gradient_filter(input,true);
  double pi180 = 180/3.14159;
  for (int r = 0; r < outputX.rows(); r++) {
    for (int c = 0; c < outputX.cols(); c++) {
      magnitude[r][c] = sqrt(pow(outputY[r][c],2)+pow(outputX[r][c],2));
      if(magnitude[r][c]<thresh){
	magnitude[r][c] = 0;
      }
      direction[r][c] = getAngle(atan2(outputY[r][c],outputX[r][c]) * pi180);
      //cout<<direction[r][c]<<endl;
    }
  }

  for (int r = 1; r < outputX.rows()-1; r++) {
    for (int c = 1; c < outputX.cols()-1; c++) {
      if(direction[r][c] == 0){
	if(magnitude[r][c]>magnitude[r-1][c] && magnitude[r][c]>magnitude[r+1][c])
	  output[r][c] = 255;
	else
	  output[r][c] = 0;
      }else if(direction[r][c] == 45){
	if(magnitude[r][c]>magnitude[r+1][c-1] && magnitude[r][c]>magnitude[r-1][c+1])
	  output[r][c] = 255;
	else
	  output[r][c] = 0;
      }else if(direction[r][c] == 90){
	if(magnitude[r][c]>magnitude[r][c-1] && magnitude[r][c]>magnitude[r][c+1])
	  output[r][c] = 255;
	else
	  output[r][c] = 0;
      }else if(direction[r][c] == 135){
	if(magnitude[r][c]>magnitude[r-1][c-1] && magnitude[r][c]>magnitude[r+1][c+1])
	  output[r][c] = 255;
	else
	  output[r][c] = 0;
      }
	
      //cout<<direction[r][c]<<endl;
    }
  }

  SImageIO::write_png_file("direction.png",direction,direction,direction);
  return output;
}


//
// This main file just outputs a few test images. You'll want to change it to do 
//  something more interesting!
//
int main(int argc, char *argv[])
{
  if(!(argc == 2))
    {
      cerr << "usage: " << argv[0] << " input_image" << endl;
      return 1;
    }

  string input_filename(argv[1]);
  SDoublePlane input_image= SImageIO::read_png_file(input_filename.c_str());
  for(int i=0;i<input_image.rows();i++){
    for(int j=0;j<input_image.cols();j++){
      //input_image[i][j]/=255;
      //cout<<input_image[i][j]<<" ";
    }
  }
  // test step 2 by applying mean filters to the input image
  SDoublePlane mean_filter(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      mean_filter[i][j] = 1/9.0;
  SDoublePlane gFilter = generateGaussian();
  SDoublePlane output_image = convolve_general(input_image, mean_filter);
  SDoublePlane output_image1 =  find_edges(output_image,100);
  for(int i=0;i<output_image1.rows();i++){
    for(int j=0;j<output_image1.cols();j++){
      //cout<< output_image1[i][j]<<" ";
    }
  }
  SImageIO::write_png_file("mean_res.png",output_image1,output_image1,output_image1);
  /*
  // randomly generate some detected cars -- you'll want to replace this
  //  with your car detection code obviously!
  vector<DetectedBox> cars;
  for(int i=0; i<10; i++)
    {
      DetectedBox s;
      s.row = rand() % input_image.rows();
      s.col = rand() % input_image.cols();
      s.width = 20;
      s.height = 20;
      s.confidence = rand();
      cars.push_back(s);
    }

  write_detection_txt("detected.txt", cars);
  write_detection_image("detected.png", cars, input_image);
  */
}
