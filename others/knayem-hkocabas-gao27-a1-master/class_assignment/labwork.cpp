//
// Watermark.cpp : Add watermark to an image, or inspect if a watermark is present.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include "SImage.h"
#include "SImageIO.h"
#include "fft.h"
#include <cmath>

using namespace std;

// This code requires that input be a *square* image, and that each dimension
//  is a power of 2; i.e. that input.width() == input.height() == 2^k, where k
//  is an integer. You'll need to pad out your image (with 0's) first if it's
//  not a square image to begin with. (Padding with 0's has no effect on the FT!)
//
// Forward FFT transform: take input image, and return real and imaginary parts.
//
void fft(const SDoublePlane &input, SDoublePlane &fft_real, SDoublePlane &fft_imag)
{
  fft_real = input;
  fft_imag = SDoublePlane(input.rows(), input.cols());

  FFT_2D(1, fft_real, fft_imag);
}

// Inverse FFT transform: take real and imaginary parts of fourier transform, and return
//  real-valued image.
//
void ifft(const SDoublePlane &input_real, const SDoublePlane &input_imag, SDoublePlane &output_real)
{
  output_real = input_real;
  SDoublePlane output_imag = input_imag;

  FFT_2D(0, output_real, output_imag);
}

// Write this in Part 1.1
SDoublePlane fft_magnitude(const SDoublePlane &fft_real, const SDoublePlane &fft_imag);

// Write this in Part 1.2
SDoublePlane remove_interference(const SDoublePlane &input);

// Write this in Part 1.3 -- add watermark N to image
SDoublePlane mark_image(const SDoublePlane &input, int N);

// Write this in Part 1.3 -- check if watermark N is in image
SDoublePlane check_image(const SDoublePlane &input, int N);


int main(int argc, char **argv)
{
  try {

    if(argc < 4)
      {
	cout << "Insufficent number of arguments; correct usage:" << endl;
	cout << "    p2 problemID inputfile outputfile" << endl;
	return -1;
      }
    
    //string part = argv[1];
    string inputFile_tiger = argv[1];
    string inputFile_dog = argv[2];
    string outputFile = argv[3];
    cout << "In1: " << inputFile_tiger << "In2: " << inputFile_dog <<"  Out: " << outputFile << endl;
    
    SDoublePlane input_image_tiger = SImageIO::read_png_file(inputFile_tiger.c_str());
    SDoublePlane input_image_dog = SImageIO::read_png_file(inputFile_dog.c_str());
    
    SDoublePlane fft_real_t= SDoublePlane(input_image_tiger.rows(), input_image_tiger.cols());
    SDoublePlane fft_imag_t= SDoublePlane(input_image_tiger.rows(), input_image_tiger.cols());
    SDoublePlane fft_real_d= SDoublePlane(input_image_dog.rows(), input_image_dog.cols());
    SDoublePlane fft_imag_d= SDoublePlane(input_image_dog.rows(), input_image_dog.cols());
      
      SDoublePlane fft_real_output= SDoublePlane(input_image_dog.rows(), input_image_dog.cols());
      SDoublePlane fft_imag_output= SDoublePlane(input_image_dog.rows(), input_image_dog.cols());
      
    fft(input_image_tiger, fft_real_t, fft_imag_t);
    fft(input_image_dog, fft_real_d, fft_imag_d);
    
    SDoublePlane output_image(input_image_tiger.rows(), input_image_tiger.cols());
//     SImageIO::write_png_file(outputFile.c_str(), fft_imag_t,fft_imag_t,fft_imag_t);
      
      cout<<input_image_tiger.rows()<<endl;
      cout<<input_image_tiger.cols()<<endl;
      
    for(int i=0; i<input_image_tiger.rows(); i++){
          for(int j=0; i<input_image_tiger.cols(); j++){
              
              //fft_real_output[i][j] = (int)floor(sqrt(fft_real_t[i][j]*fft_real_t[i][j] + fft_imag_t[i][j]*fft_imag_t[i][j]));
              
              cout<<"T:"<<sqrt(fft_real_t[i][j]*fft_real_t[i][j] + fft_imag_t[i][j]*fft_imag_t[i][j])<<",";
              //cout<<"D:"<<(fft_real_d[i][j]*fft_real_d[i][j] + fft_imag_d[i][j]*fft_imag_d[i][j])<<endl;
//              if( (fft_real_t[i][j]*fft_real_t[i][j] + fft_imag_t[i][j]*fft_imag_t[i][j])< 16384.0){
//                  fft_real_output[i][j] = fft_real_t[i][j];
//                  fft_imag_output[i][j] = fft_imag_t[i][j];
//              }
//              else if( (fft_real_d[i][j]*fft_real_d[i][j] + fft_imag_d[i][j]*fft_imag_d[i][j])> 16384.0){
//                  fft_real_output[i][j] = fft_real_d[i][j];
//                  fft_imag_output[i][j] = fft_imag_d[i][j];
//              }
//              else{
//                  fft_real_output[i][j] = 0;
//                  fft_imag_output[i][j] = 0;
//              }
          }
        cout<<endl;
      }
      SImageIO::write_png_file(outputFile.c_str(), fft_real_output,fft_real_output,fft_real_output);
      cout<<"oka"<<endl;
      
//      SDoublePlane output_image(input_image_tiger.rows(), input_image_tiger.cols());
//      ifft(fft_real_output, fft_imag_output, output_image);
//      SImageIO::write_png_file("output.png", output_image,output_image,output_image);
      
      //cout << fft_real <<endl;
      
    /*
    if(part == "1.1")
      {
	// do something here!
      }
    else if(part == "1.2")
      {
	// do something here!
      }
    else if(part == "1.3")
      {
	if(argc < 6)
	  {
	    cout << "Need 6 parameters for watermark part:" << endl;
	    cout << "    p2 1.3 inputfile outputfile operation N" << endl;
	    return -1;
	  }
	string op(argv[4]);
	if(op == "add")
	  {
	    // add watermark
	  }
	else if(op == "check")
	  {
	    // check watermark
	  }
	else
	  throw string("Bad operation!");
       
	int N = atoi(argv[5]);
      }
    else
      throw string("Bad part!");
     */

  } 
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}








