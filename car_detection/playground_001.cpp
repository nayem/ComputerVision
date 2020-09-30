//
// detect.cpp : Detect cars in satellite images.
//
// Based on skeleton code by D. Crandall, Spring 2017
//
// PUT YOUR NAMES HERE
//
//

#include "SImage.h"
#include "SImageIO.h"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

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




/***************** New Functions: Start ********************/

const int THRESHOLD_MODE_NONE = 0;
const int THRESHOLD_MODE_MAXMIN = 1;
const int THRESHOLD_MODE_MEDIAN = 2;
const int THRESHOLD_MODE_VALUE = 3;

// Initialize a SDoublePlane
void initializeSDoublePlane(const SDoublePlane &input, double vec[])
{
    SDoublePlane output(input.rows(), input.cols() );
    
    for(int i=0, index=0; i < input.rows(); ++i)     // rows
        for(int j=0; j < input.cols(); ++j)          // columns
            input[i][j]=vec[index++];
}

// Transpose a SDoublePlane
SDoublePlane transpose(const SDoublePlane &input)
{
    SDoublePlane output(input.cols(), input.rows());
    
    for(int i=0; i < input.rows(); ++i)              // rows
        for(int j=0; j < input.cols(); ++j)          // columns
            output[j][i] = input[i][j];
    
    return output;
}

// Normalize a SDoublePlane
SDoublePlane normalize(const SDoublePlane &input, double scale=1.0)
{
    SDoublePlane output(input.rows(), input.cols());
    double min=1E10, max=-1E10;
    
    for(int i=0; i < input.rows(); ++i) {             // rows
        for(int j=0; j < input.cols(); ++j) {         // columns
            if (max<input[i][j])
                max = input[i][j];
            if (min>input[i][j])
                min = input[i][j];
        }
    }
    
    for(int i=0; i < input.rows(); ++i)
        for(int j=0; j < input.cols(); ++j)
            output[i][j] = scale*(input[i][j]-min)/(max-min);
    
    return output;
}

// Calculate Magnitude and Gradient using various Threshold value
SDoublePlane magnitude_gradient(const SDoublePlane &magnitude, const SDoublePlane &gx,const SDoublePlane &gy, int mode=THRESHOLD_MODE_NONE, double threshold_value = 0.0)
{
    SDoublePlane gradient(magnitude.rows(), magnitude.cols() );
    double min=1E10, max=-1E10, sum=0.0;
    
    for(int i=0; i < magnitude.rows(); ++i){              // rows
        for(int j=0; j < magnitude.cols(); ++j){          // columns
            magnitude[i][j] = sqrt(gx[i][j]*gx[i][j]+gy[i][j]*gy[i][j]);
            gradient[i][j] = atan2(gy[i][j],gx[i][j])*180/3.1416;  // radian (.. * 180 / PI; for degree)
            
            switch (mode) {
                case THRESHOLD_MODE_MAXMIN:
                case THRESHOLD_MODE_VALUE:
                    if (max<magnitude[i][j])
                        max = magnitude[i][j];
                    if (min>magnitude[i][j])
                        min = magnitude[i][j];
                    break;
                    
                case THRESHOLD_MODE_MEDIAN:
                    sum += magnitude[i][j];
                    break;
            }
        }
    }
    
    //    cout<<"Max: "<<max<<endl;
    //    cout<<"Min: "<<min<<endl;
    
    // Recompute Magnitude
    for(int i=0; mode!= THRESHOLD_MODE_NONE && i < magnitude.rows(); ++i){              // rows
        for(int j=0; j < magnitude.cols(); ++j){          // columns
            
            switch (mode) {
                case THRESHOLD_MODE_MAXMIN:
                    if ( magnitude[i][j] >= (max+min)/2.0 )
                        //magnitude[i][j] = 255.0;
                        magnitude[i][j] = magnitude[i][j];
                    else
                        magnitude[i][j] = 0.0;
                    break;
                    
                case THRESHOLD_MODE_MEDIAN:
                    if ( magnitude[i][j] >= sum/ (magnitude.rows()*magnitude.cols()) )
                        //magnitude[i][j] = 255.0;
                        magnitude[i][j] = magnitude[i][j];
                    else
                        magnitude[i][j] = 0.0;
                    break;
                    
                case THRESHOLD_MODE_VALUE:
                    if (  abs( (magnitude[i][j]-abs(min))/(max-abs(min))) >= threshold_value )
                        magnitude[i][j] = 1.0;
                    //                        magnitude[i][j] = (magnitude[i][j]-min)/(max-min);
                    else
                        magnitude[i][j] = 0.0;
                    //                        magnitude[i][j] = (magnitude[i][j]-min)/(max-min);
                    break;
            }
            
        }
    }
    
    return gradient;
}

int getAngle(double degree){
    int angle;
    if((degree>=22.5 && degree<67.5) || (degree<-112.5 && degree>=-157.5))
        angle = 45;
    if((degree>=67.5 && degree<112.5) || (degree<-67.5 && degree>=-112.5))
        angle = 90;
    if((degree>=112.5 && degree<157.5) || (degree<-22.5 && degree>=-67.5))
        angle = 135;
    if((degree>=-22.5 && degree<22.5) || (degree<-157.5 || degree>=157.5))
        angle = 0;
    return angle;
}



// Non-Maximum Suppression
SDoublePlane nonMaximumSuppression(const SDoublePlane &magnitude, const SDoublePlane &gradient, const SDoublePlane &gx,const SDoublePlane &gy){
    
    SDoublePlane output(magnitude.rows(), magnitude.cols() );
    
    double y[][2]={0,0,0,0};
    double x_est = 0;
    
    ////////////////////////////////
    //        for (int r = 1; r < magnitude.rows()-1; r++) {
    //            for (int c = 1; c < magnitude.cols()-1; c++) {
    //                if(getAngle(gradient[r][c]) == 0){
    //                    if(magnitude[r][c]>magnitude[r-1][c] && magnitude[r][c]>magnitude[r+1][c])
    //                        output[r][c] = 1;
    //                    else
    //                        output[r][c] = 0;
    //                }else if(getAngle(gradient[r][c]) == 45){
    //                    if(magnitude[r][c]>magnitude[r+1][c-1] && magnitude[r][c]>magnitude[r-1][c+1])
    //                        output[r][c] = 1;
    //                    else
    //                        output[r][c] = 0;
    //                }else if(getAngle(gradient[r][c]) == 90){
    //                    if(magnitude[r][c]>magnitude[r][c-1] && magnitude[r][c]>magnitude[r][c+1])
    //                        output[r][c] = 1;
    //                    else
    //                        output[r][c] = 0;
    //                }else if(getAngle(gradient[r][c]) == 135){
    //                    if(magnitude[r][c]>magnitude[r-1][c-1] && magnitude[r][c]>magnitude[r+1][c+1])
    //                        output[r][c] = 1;
    //                    else
    //                        output[r][c] = 0;
    //                }
    //
    //                //cout<<direction[r][c]<<endl;
    //            }
    //        }
    ///////////////////////////////
    
    for (int i=1; i< magnitude.rows()-1;++i){ // row
        for (int j=1; j< magnitude.cols()-1;++j){  //col
            if ((gradient[i][j]>=0.0 && gradient[i][j]<=45.0) || (gradient[i][j]<-135.0 && gradient[i][j]>=-180.0)){
                y[0][0] = magnitude[i][j-1];
                y[0][1] = magnitude[i-1][j-1];
                
                y[1][0] = magnitude[i][j+1];
                y[1][1] = magnitude[i+1][j+1];
                
                x_est = abs(gy[i][j]/magnitude[i][j]); //y
                if ((magnitude[i][j] >= (y[1][1]-y[1][0])*x_est+y[1][0] ) && (magnitude[i][j] >= (y[0][1]-y[0][0])*x_est+y[0][0] )) //interpolation
                    output[i][j]= magnitude[i][j];
                else
                    output[i][j]=0;
            }
            else if ((gradient[i][j]>45.0 && gradient[i][j]<=90.0) || (gradient[i][j]<-90.0 && gradient[i][j]>=-135.0)){
                y[0][0] = magnitude[i-1][j];
                y[0][1] = magnitude[i-1][j-1];
                
                y[1][0] = magnitude[i+1][j];
                y[1][1] = magnitude[i+1][j+1];
                
                x_est = abs(gx[i][j]/magnitude[i][j]);
                
                if ((magnitude[i][j] >= (y[1][1]-y[1][0])*x_est+y[1][0] ) && (magnitude[i][j] >= (y[0][1]-y[0][0])*x_est+y[0][0] )) //interpolation
                    output[i][j]= magnitude[i][j];
                else
                    output[i][j]=0;
            }
            else if ((gradient[i][j]>90.0 && gradient[i][j]<=135.0) || (gradient[i][j]<-45.0 && gradient[i][j]>=-90.0)){
                y[0][0] = magnitude[i-1][j];
                y[0][1] = magnitude[i-1][j+1];
                
                y[1][0] = magnitude[i+1][j];
                y[1][1] = magnitude[i+1][j-1];
                
                x_est = abs(gx[i][j]/magnitude[i][j]);
                
                if ((magnitude[i][j] >= (y[1][1]-y[1][0])*x_est+y[1][0] ) && (magnitude[i][j] >= (y[0][1]-y[0][0])*x_est+y[0][0] )) //interpolation
                    output[i][j]= magnitude[i][j];
                else
                    output[i][j]=0;
            }
            else if ((gradient[i][j]>135.0 && gradient[i][j]<=180.0) || (gradient[i][j]<0.0 && gradient[i][j]>=45.0)){
                y[0][0] = magnitude[i][j+1];
                y[0][1] = magnitude[i-1][j+1];
                
                y[1][0] = magnitude[i][j-1];
                y[1][1] = magnitude[i+1][j-1];
                
                x_est = abs(gx[i][j]/magnitude[i][j]);
                
                if ((magnitude[i][j] >= (y[1][1]-y[1][0])*x_est+y[1][0] ) && (magnitude[i][j] >= (y[0][1]-y[0][0])*x_est+y[0][0] )) //interpolation
                    output[i][j]= magnitude[i][j];
                else
                    output[i][j]=0;
            }
        }
    }
    return output;
    
}

void findConnectedWeakEdges(const SDoublePlane &magnitude, int row, int col);

// Double Threshold and Edge Tracking
SDoublePlane doubleThreshold(const SDoublePlane &magnitude, double upperThreshold, double lowerThreshold){
    
    double max = -1E10;
    
    for (int i=0; i< magnitude.rows();++i)
        for (int j=0; j< magnitude.cols();++j)
            if(max < magnitude[i][j])
                max =magnitude[i][j];
    
    cout<<"Max: "<<max<<endl;
    
    double high_threshold = max*upperThreshold;
    double low_threshold = high_threshold*lowerThreshold;
    
    cout<<"high_threshold: "<<high_threshold<<endl;
    cout<<"low_threshold: "<<low_threshold<<endl;
    
    int strongIndex=1, weakIndex = 1;
    
    vector<int> strongEdgesRow(magnitude.rows()*magnitude.cols(),0); // Keep track of the strong edge row index
    vector<int> strongEdgesCol(magnitude.rows()*magnitude.cols(),0); //  Keep track of the strong edge col index
    vector<int> weakEdgesRow(magnitude.rows()*magnitude.cols(),0); //  Keep track of the weak edge row index
    vector<int> weakEdgesCol(magnitude.rows()*magnitude.cols(),0); //  Keep track of the weak edge col index
    
    for (int i=1; i< magnitude.rows()-1; ++i){ // row
        for (int j=1; j< magnitude.cols()-1; ++j){  //col
            
            if (magnitude[i][j] > high_threshold){      // Strong edge
                magnitude[i][j] = 1;
                
                strongEdgesRow[strongIndex] = i;
                strongEdgesCol[strongIndex++] = j;
            }
            else if (magnitude[i][j] < low_threshold)   //No edge
                magnitude[i][j] = 0;
            else{                                       //Weak edge
                weakEdgesRow[weakIndex] = i;
                weakEdgesCol[weakIndex++] = j;
            }
        }
    }
    cout<<"strongIndex: "<<strongIndex<<endl;
    for (int limit=0; limit<10000; ++limit) {
        for (int i=0; i<strongIndex; ++i) {
            // Find the weak edges that are connected to strong edges and set them to 1
            findConnectedWeakEdges(magnitude, strongEdgesRow[i], strongEdgesCol[i]);
        }
        
        // Remove the remaining weak edges that are not actually edges and is noise instead
        for (int i=0; i<weakIndex; ++i) {
            if ( abs(magnitude[weakEdgesRow[i]][weakEdgesCol[i]] - 1)<0.001 )
                magnitude[weakEdgesRow[i]][weakEdgesCol[i]] = 0;
        }
    }
    
    return magnitude;
}

// Find weak edges that are connected to strong edges and set them to 1
void findConnectedWeakEdges(const SDoublePlane &magnitude, int row, int col){
    
    for( int i= -3; i<=3 ; ++i){
        for( int j= -3; j<=3 ; ++j){
            if ( (row+i > 0) && (col+j > 0) && (row+i < magnitude.rows()) && (col+j < magnitude.cols()) ){ // Make sure we are not out of bounds
                if ( (magnitude[row+i][col+j] > 0) && (magnitude[row+i][col+j] < 1) ){
                    magnitude[row+i][col+j] = 1;
                    findConnectedWeakEdges(magnitude, row+i, col+j);
                }
            }
        }
    }
    
}

// Copy a SDoublePlane
SDoublePlane cutSDoublePlane(const SDoublePlane &input, int row_start, int row_end, int col_start, int col_end){
    
    SDoublePlane output(row_end-row_start+1, col_end-col_start+1 );
    //    cout<<"CUT::\nrow_start: "<<row_start<<endl;
    //    cout<<"row_end: "<<row_end<<endl;
    //    cout<<"col_start: "<<col_start<<endl;
    //    cout<<"col_end: "<<col_end<<endl;
    
    for(int i=row_start,p=0 ; i <= row_end; ++i,++p)     // rows
        for(int j=col_start,q=0; j <= col_end; ++j,++q){          // columns
            //            cout<<"i:"<<i<<",j:"<<j<<",p"<<p<<",q:"<<q<<endl;
            //            cout<<"output[p][q]: "<<output[p][q]<<",input[i][j]: "<<input[i][j]<<endl;
            output[p][q]=input[i][j];
        }
    
    return output;
}

void gradientHistogram(const SDoublePlane &gradient, vector<int> &histogram, int part=0){
    
    for(int i=0 ; i < gradient.rows(); ++i){
        for(int j=0; j < gradient.cols(); ++j){
            //            cout<<gradient[i][j]<<endl;
            if ( (gradient[i][j] > -22.5) && (gradient[i][j] <= 22.5 ) ){ // 0-degree(-22.5< <=22.5)
                //                cout<<"[1]histogram 0 : "<< histogram[ (part*8 + 0)]<<",";
                histogram[ part*8 + 0]++;
                //                cout<<"[2]histogram 0 : "<< histogram[ part*8 + 0]<<endl;
            }
            if ( (gradient[i][j] > 22.5) && (gradient[i][j] <= 67.5 ) ){ // 45-degree(22.5< <=67.5)
                histogram[ part*8 +1]++;
            }
            if ( (gradient[i][j] > 67.5) && (gradient[i][j] <= 112.5 ) ){ // 90-degree(67.5< <=112.5)
                histogram[ part*8 +2]++;
            }
            if ( (gradient[i][j] > 112.5) && (gradient[i][j] <= 157.5 ) ){ // 135-degree(112.5< <=157.5)
                histogram[part*8 +3]++;
            }
            if (( (gradient[i][j] > 157.5) && (gradient[i][j] <= 180 ) ) || ( (gradient[i][j] >= -180) && (gradient[i][j] <= 157.5 ) )){ // 180/-180-degree(157.5< <=180 & -180< <=-157.5)
                histogram[part*8 +4]++;
            }
            if ( (gradient[i][j] > -157.5) && (gradient[i][j] <= -112.5 ) ){ // -135-degree(-157.5< <=-112.5)
                histogram[part*8 +5]++;
            }
            if ( (gradient[i][j] > -112.5) && (gradient[i][j] <= -67.5 ) ){ // -90-degree(-112.5< <=-67.5)
                histogram[part*8 +6]++;
                //                cout<<"[2]histogram -90 : "<< histogram[ part*8 + 6]<<endl;
            }
            if ( (gradient[i][j] > -67.5) && (gradient[i][j] <= -22.5 ) ){ // -45-degree(-67.5< <=-22.5)
                histogram[part*8 +7]++;
            }
            
        }
    }
    
    //    for (int p=0; p<4; ++p){
    //        for (int q=0; q<8; ++q)
    //            cout<<histogram[p*8+q]<<",";
    //        cout<<endl;
    //    }
    //    cout<<"********"<<endl;
    
}

bool isOldCar(vector<DetectedBox> cars, int window_row_size, int window_col_size, int row, int col){
    int flag = true;
    for (int i=0; i<cars.size(); i++){
//        if ( cars[i].row < row && (cars[i].row+cars[i].width) > row && ((cars[i].col < col && (cars[i].col+cars[i].height) > col) || (cars[i].col < col+window_col_size && (cars[i].col+cars[i].height) > col+window_col_size)) ){
//            flag = false;
//            break;
//        }
        bool val1 = cars[i].row < row && (cars[i].row+cars[i].width) > row;
        
        bool val2 = (cars[i].col < col && (cars[i].col+cars[i].height) > col) || (cars[i].col < col+window_row_size && (cars[i].col+cars[i].height) > col+window_row_size);
        
        if ( val1 && val2 ){
            flag = false;
            break;
        }
    }
    
    return flag;
}

/***************** End ****************/





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
    
    int rows = input.rows();
    int cols = input.cols();
    int kRows = col_filter.rows();  // 1 col
    int kCols = row_filter.cols();  // 1 row
    
    // Convolution code here
    int kCenterX = kCols / 2;
    int kCenterY = kRows / 2;
    
    // Apply Column Filter (1 col, multiple rows)
    for(int i=0; i < rows; ++i)              // rows
    {
        for(int j=0; j < cols; ++j)          // columns
        {
            for(int m=0; m < kRows; ++m)     // kernel rows
            {
                int mm = kRows - 1 - m;      // row index of flipped kernel
                
                // index of input signal, used for checking boundary
                int ii = i + (m - kCenterY);
                
                // ignore input samples which are out of bound
                if( ii >= 0 && ii < rows )
                    temp[i][j] += input[ii][j] * col_filter[mm][0];
            }
        }
    }
    
    // Apply Row Filter (1 row, multiple cols)
    for(int i=0; i < rows; ++i)              // rows
    {
        for(int j=0; j < cols; ++j)          // columns
        {
            for(int n=0; n < kCols; ++n)     // kernel coss
            {
                int nn = kCols - 1 - n;  // column index of flipped kernel
                
                // index of input signal, used for checking boundary
                int jj = j + (n - kCenterX);
                
                // ignore input samples which are out of bound
                if(jj >= 0 && jj < cols )
                    output[i][j] += temp[i][jj] * row_filter[0][nn];
            }
        }
    }
    
    return output;
}

// Convolve an image with a  convolution kernel
//
SDoublePlane convolve_general(const SDoublePlane &input, const SDoublePlane &filter)
{
    SDoublePlane output(input.rows(), input.cols());
    
    //2D convolution-slower version
    int rows = input.rows();
    int cols = input.cols();
    int kRows = filter.rows();
    int kCols = filter.cols();
    
    // Convolution code here
    int kCenterX = kCols / 2;
    int kCenterY = kRows / 2;
    
    for(int i=0; i < rows; ++i)              // rows
    {
        for(int j=0; j < cols; ++j)          // columns
        {
            for(int m=0; m < kRows; ++m)     // kernel rows
            {
                int mm = kRows - 1 - m;      // row index of flipped kernel
                
                for(int n=0; n < kCols; ++n) // kernel columns
                {
                    int nn = kCols - 1 - n;  // column index of flipped kernel
                    
                    // index of input signal, used for checking boundary
                    int ii = i + (m - kCenterY);
                    int jj = j + (n - kCenterX);
                    
                    // ignore input samples which are out of bound
                    if( ii >= 0 && ii < rows && jj >= 0 && jj < cols )
                        output[i][j] += input[ii][jj] * filter[mm][nn];
                }
            }
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
    SDoublePlane col_filter(3,1);
    SDoublePlane row_filter(1,3);
    
    double filter_values1[] = {1,2,1};
    initializeSDoublePlane(col_filter, filter_values1);
    
    if (_gx == true){
        double filter_values2[] = {-1,0,1};
        initializeSDoublePlane(row_filter, filter_values2);
        output = convolve_separable(input,row_filter,col_filter );
    }
    else{
        double filter_values2[] = {1,0,-1};
        initializeSDoublePlane(row_filter, filter_values2);
        output = convolve_separable(input,transpose(col_filter), transpose(row_filter) );
    }
    
    return output;
}

// Apply an edge detector to an image, returns the binary edge map
//
SDoublePlane find_edges(const SDoublePlane &input, double thresh=0)
{
    SDoublePlane output(input.rows(), input.cols());
    
    // Implement an edge detector of your choice, e.g.
    // use your sobel gradient operator to compute the gradient magnitude and threshold
    
    // Apply Gaussian filter
    SDoublePlane col_filter(7,1);
    
    double gaussian_filter_values[] = {0.006,0.061,0.242,0.383,0.242,0.061,0.006};
    initializeSDoublePlane(col_filter, gaussian_filter_values);
    
    SDoublePlane gaussianBlur_image = convolve_separable(input,transpose(col_filter),col_filter );
    SImageIO::write_png_file("detected-2-GaussianBlur.png", gaussianBlur_image,gaussianBlur_image,gaussianBlur_image);
    
    // Sobel filter: Gx, Gy
    SDoublePlane sobelGx_image = sobel_gradient_filter(gaussianBlur_image, true);
    SImageIO::write_png_file("detected-3-Gx.png", sobelGx_image,sobelGx_image,sobelGx_image);
    
    SDoublePlane sobelGy_image = sobel_gradient_filter(gaussianBlur_image, false);
    SImageIO::write_png_file("detected-4-Gy.png", sobelGy_image,sobelGy_image,sobelGy_image);
    
    // Sobel Magnitude & Gradient
    //    SDoublePlane magnitude_image(input_image.rows(), input_image.cols());
    
    SDoublePlane gradient_image = magnitude_gradient(output,sobelGx_image,sobelGy_image,THRESHOLD_MODE_VALUE,thresh);
    //    SDoublePlane gradient_image = magnitude_gradient(magnitude_image,sobelGx_image,sobelGy_image,THRESHOLD_MODE_VALUE,0.10);
    SImageIO::write_png_file("detected-5-magnitude.png", normalize(output,255),normalize(output,255),normalize(output,255) );
    //    SImageIO::write_png_file("detected-5-magnitude.png", normalize(magnitude_image,255),normalize(magnitude_image,255),normalize(magnitude_image,255) );
    SImageIO::write_png_file("detected-6-gradient.png", gradient_image,gradient_image,gradient_image);
    
    //    SDoublePlane nms_image = nonMaximumSuppression(output,gradient_image, sobelGx_image,sobelGy_image);
    //    //    SDoublePlane nms_image = nonMaximumSuppression(magnitude_image,gradient_image, sobelGx_image,sobelGy_image);
    //    SImageIO::write_png_file("detected-7-nms.png", normalize(nms_image,255),normalize(nms_image,255),normalize(nms_image,255));
    //    //
    //    SDoublePlane dT_image = doubleThreshold(nms_image,0.28,0.15);
    //    SImageIO::write_png_file("detected-8-dT.png", normalize(dT_image,255),normalize(dT_image,255),normalize(dT_image,255));
    cout<<"oka-2"<<endl;
    //    output = doubleThreshold(nms_image,0.2,0.05);
    //    SImageIO::write_png_file("detected-8-dT.png", normalize(output,255),normalize(output,255),normalize(output,255));
    
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
    
    // test step 2 by applying mean filters to the input image
    //    SDoublePlane mean_filter(3,3);
    //    for(int i=0; i<3; i++)
    //        for(int j=0; j<3; j++)
    //            mean_filter[i][j] = 1/9.0;
    //    SDoublePlane output_image = convolve_general(input_image, mean_filter);
    
    SImageIO::write_png_file("detected-1-original.png", input_image,input_image,input_image);
    cout<<"oka-1"<<endl;
    
    // test step 2 by applying mean filters to the input image
    SDoublePlane magnitude_image=find_edges(input_image, 0.13);
    
    // Train by Car image
    int Number_of_car_image = 17;
    string image_name;
    SDoublePlane part[4];
    SDoublePlane car_image, sobelGx_image_p,sobelGy_image_p,gradient_image_p;
    vector<int> car_histogram(4*8, 0);
    
    for (int i=1; i<Number_of_car_image; ++i) {
        image_name = "car_image/"+to_string(i)+".png";
        car_image= SImageIO::read_png_file( image_name.c_str());
        
        int rows = car_image.rows();
        int cols = car_image.cols();
        cout<<"rows: "<<rows<<endl;
        cout<<"cols: "<<cols<<endl;
        // 4part of a car image
        part[0] = cutSDoublePlane(car_image, 0, (int)floor(rows/2.0), 0, (int)floor(cols/2.0) );
        part[1] = cutSDoublePlane(car_image, (int)ceil(rows/2.0), rows-1, 0, (int)floor(cols/2.0) );
        part[2] = cutSDoublePlane(car_image, 0, (int)floor(rows/2), (int)ceil(cols/2), cols-1 );
        part[3] = cutSDoublePlane(car_image, (int)ceil(rows/2), rows-1, (int)ceil(cols/2), cols-1 );
        
        for (int j=0; j<4; ++j) {
            sobelGx_image_p = sobel_gradient_filter(part[j], true);
            sobelGy_image_p = sobel_gradient_filter(part[j], false);
            SDoublePlane magnitude_image_part(part[j].rows(), part[j].cols());
            
            gradient_image_p = magnitude_gradient(magnitude_image_part,sobelGx_image_p,sobelGy_image_p,THRESHOLD_MODE_VALUE,0.25);
            SImageIO::write_png_file("detected-10-magnitude.png", normalize(magnitude_image_part,255),normalize(magnitude_image_part,255),normalize(magnitude_image_part,255) );
            SImageIO::write_png_file("detected-10-gradient.png", gradient_image_p,gradient_image_p,gradient_image_p);
            
            gradientHistogram(gradient_image_p,car_histogram,j);
            //            for (int p=0; p<4; ++p){
            //                for (int q=0; q<8; ++q)
            //                    cout<<car_histogram[p*8+q]<<",";
            //                cout<<endl;
            //            }
        }
        
    }
    // Average histrogram
    for (int i=0; i<32; ++i)
        car_histogram[i]=(int)round(car_histogram[i]/16.0);
    cout<<"FINAL"<<endl;
    for (int p=0; p<4; ++p){
        for (int q=0; q<8; ++q)
            cout<<car_histogram[p*8+q]<<",";
        cout<<endl;
    }
    
    // Sliding window
    int window_row_size = 25;
    int window_col_size = 40;
    SDoublePlane window;
    SDoublePlane sobelGx_image_w,sobelGy_image_w,gradient_image_w;
    vector<int> window_histogram(4*8, 0);
    
    vector<DetectedBox> cars;
    cout<<"input_image.rows():"<<input_image.rows()<<",input_image.rows()-window_row_size:"<<input_image.rows()-window_row_size<<endl;
    cout<<"input_image.cols():"<<input_image.cols()<<",input_image.cols()-window_col_size:"<<input_image.cols()-window_col_size<<endl;
    cout<<"magnitude_image.rows():"<<magnitude_image.rows()<<endl;
    cout<<"magnitude_image.cols():"<<magnitude_image.cols()<<endl;
    
    for(int i=0; i < (magnitude_image.rows()-window_row_size); i=i+4)              // rows
    {
        for(int j=0; j < (magnitude_image.cols()-window_col_size); j=j+4){          // columns
            
            window = cutSDoublePlane(magnitude_image, i, window_row_size+i, j, window_col_size+j );
            
            // 4part of a window
            part[0] = cutSDoublePlane(window, 0, (int)floor(window_row_size/2.0), 0, (int)floor(window_col_size/2.0) );
            part[1] = cutSDoublePlane(window, (int)ceil(window_row_size/2.0), window_row_size-1, 0, (int)floor(window_col_size/2.0) );
            part[2] = cutSDoublePlane(window, 0, (int)floor(window_row_size/2), (int)ceil(window_col_size/2), window_col_size-1 );
            part[3] = cutSDoublePlane(window, (int)ceil(window_row_size/2), window_row_size-1, (int)ceil(window_col_size/2), window_col_size-1 );
            
            for (int p=0; p<4; ++p) {
                sobelGx_image_w = sobel_gradient_filter(part[p], true);
                sobelGy_image_w = sobel_gradient_filter(part[p], false);
                
                SDoublePlane magnitude_image_window(part[p].rows(), part[p].cols());
                
                gradient_image_w = magnitude_gradient(magnitude_image_window,sobelGx_image_w,sobelGy_image_w,THRESHOLD_MODE_VALUE,0.17);
                
                gradientHistogram(gradient_image_w,window_histogram,p);
                
                //                for (int p=0; p<4; ++p){
                //                    for (int q=0; q<8; ++q)
                //                        cout<<car_histogram[p*8+q]<<",";
                //                    cout<<endl;
                //                }
            }
            cout<<"i:"<<i<<",j:"<<j<<endl;
            
            for (int p=0; p<4; ++p){
                for (int q=0; q<8; ++q)
                    cout<<window_histogram[p*8+q]<<",";
                cout<<endl;
            }
            
            // calculate Eucledian distance
            double sum=0;
            for (int p=0; p<32; ++p) {
                sum += (car_histogram[p]-window_histogram[p])*(car_histogram[p]-window_histogram[p]);
            }
            sum = sqrt(sum);
            cout<<sum<<endl;
            if(sum < 80 && isOldCar(cars, window_row_size, window_col_size, i, j) ){
                DetectedBox s;
                s.row = i;
                s.col = j;
                s.width = window_row_size;
                s.height = window_col_size;
                s.confidence = 100.0/sum;
                cars.push_back(s);
                
                j +=(window_col_size-4);
                cout<< " %%%%%%%%%%%%%% " << cars.size()<< "%%%%%%%%"<<endl;
            }
            window_histogram = vector<int>(4*8,0);
            
        }
    }
    
    
    
    
    // randomly generate some detected cars -- you'll want to replace this
    //  with your car detection code obviously!
    //    vector<DetectedBox> cars;
    //    for(int i=0; i<10; i++)
    //    {
    //        DetectedBox s;
    //        s.row = rand() % input_image.rows();
    //        s.col = rand() % input_image.cols();
    //        s.width = 20;
    //        s.height = 20;
    //        s.confidence = rand();
    //        cars.push_back(s);
    //    }
    //
    write_detection_txt("detected.txt", cars);
    write_detection_image("detected.png", cars, input_image);
    
}


