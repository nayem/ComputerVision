// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "Sift.h"
#include <vector>
#include <algorithm>
#include <ctime>
#include "levmarq.h"
#include <numeric>
#include<string.h>
#include<fstream>
#include<sstream>
#include<dirent.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;


void warp_image(const CImg<double> &in, CImg<double> &out, double transform_matrix[][3])
{
    int x, y;
    int dimx = in.width(), dimy = in.height();
    for(int k = 0; k < 3; k++)
    {
        cimg_forXY(out,x,y)
        {
            double z = transform_matrix[2][0] * x + transform_matrix[2][1] * y + transform_matrix[2][2];
            double _x = (transform_matrix[0][0] * x + transform_matrix[0][1] * y + transform_matrix[0][2]) / z;
            double _y = (transform_matrix[1][0] * x + transform_matrix[1][1] * y + transform_matrix[1][2]) / z;
            
            if(_x - floor(_x) > 0 || _y - floor(_y) > 0)
            {
                //bilinearly interpolation
                //inspired from https://www.mathworks.com/matlabcentral/fileexchange/48135-image-warping
                double x_dash = floor(_x);
                double y_dash = floor(_y);
                double x_bar = x_dash + 1;
                double y_bar = y_dash + 1;
                
                double Ex_bar = x_bar - _x;
                double Ey_bar = y_bar - _y;
                double Ex_dash = _x - x_dash;
                double Ey_dash = _y - y_dash;
                if(x_dash > 0 && x_bar < dimx && y_dash > 0 && y_bar < dimy)
                {
                    double tmp = Ex_bar * Ey_bar * in(x_dash,y_dash,k);
                    tmp = tmp + Ex_dash * Ey_bar * in(x_bar,y_dash,k);
                    tmp = tmp + Ex_bar * Ey_dash * in(x_dash,y_bar,k);
                    tmp = tmp + Ex_dash * Ey_dash * in(x_bar,y_bar,k);
                    out(x,y,0,k) = tmp;
                }
            }
            else
            {
                out(x, y, 0, 0) = in(_x, _y, 0, k);
            }
        }
    }
}

void homographies(const double correspondence_input[4][2][2], CImg <double> &H)
{
    double xcentroid1;
    double ycentroid1;
    double xcentroid2;
    double ycentroid2;
    CImg <double> tr1(3,3);
    CImg <double> tr2(3,3);
    CImg <double> sc1(3,3);
    CImg <double> sc2(3,3);
    
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  FORMAT OF THE CORRESPONDENCE MATRIX IS AS FOLLOWS:
    //  correspondence_input[i][j][k]
    //  i = correspondence_index
    //  j = 1 for image1
    //  j = 2 for image2
    //  k = 1 for xcoordinate
    //  k = 2 for ycoordinate
    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //  NORMALIZE THE CORRESPONDENCES
    //  Find the centroid of the correspondences
    xcentroid1 = 0.0;
    ycentroid1 = 0.0;
    xcentroid2 = 0.0;
    ycentroid2 = 0.0;
    for (int rtemp = 0; rtemp < 4; rtemp++)
    {
        xcentroid1 += correspondence_input[rtemp][0][0];
        ycentroid1 += correspondence_input[rtemp][0][1];
        xcentroid2 += correspondence_input[rtemp][1][0];
        ycentroid2 += correspondence_input[rtemp][1][1];
    }
    
    xcentroid1 /= 4.0;
    ycentroid1 /= 4.0;
    xcentroid2 /= 4.0;
    ycentroid2 /= 4.0;
    
    //  Generating the Translation Matrices
    tr1(0,0) = 1.0;
    tr1(1,0) = 0.0;
    tr1(2,0) = -xcentroid1;
    tr1(0,1) = 0.0;
    tr1(1,1) = 1.0;
    tr1(2,1) = -ycentroid1;
    tr1(0,2) = 0.0;
    tr1(1,2) = 0.0;
    tr1(2,2) = 1;
    tr2(0,0) = 1.0;
    tr2(1,0) = 0.0;
    tr2(2,0) = -xcentroid2;
    tr2(0,1) = 0.0;
    tr2(1,1) = 1.0;
    tr2(2,1) = -ycentroid2;
    tr2(0,2) = 0.0;
    tr2(1,2) = 0.0;
    tr2(2,2) = 1;
    //cout <<"0" << endl;
    
    //  Generating Translated Correspondences
    double correspondence_input_tr[4][2][2] = {0};
    for(int i = 0; i < 4; i++)
    {
        CImg<double> temp_res(1,3);
        CImg<double> temp(1,3);
        temp(0,0) = correspondence_input[i][0][0];
        temp(0,1) = correspondence_input[i][0][1];
        temp(0,2) = 1;
        temp_res = tr1 * temp;
        correspondence_input_tr[i][0][0] = temp_res(0,0);
        correspondence_input_tr[i][0][1] = temp_res(0,1);
        
        temp(0,0) = correspondence_input[i][1][0];
        temp(0,1) = correspondence_input[i][1][1];
        temp(0,2) = 1;
        temp_res = tr2 * temp;
        correspondence_input_tr[i][1][0] = temp_res(0,0);
        correspondence_input_tr[i][1][1] = temp_res(0,1);
        
    }
    
    //cout <<"1" << endl;
    
    double rms1 = 0;
    double rms2 = 0;
    //  Computing the average RMS distance of translated correspondences
    for(int i = 0; i < 4; i++)
    {
        rms1 = rms1 + sqrt(correspondence_input_tr[i][0][0] * correspondence_input_tr[i][0][0] + correspondence_input_tr[i][0][1] * correspondence_input_tr[i][0][1]);
        rms2 = rms2 + sqrt(correspondence_input_tr[i][1][0] * correspondence_input_tr[i][1][0] + correspondence_input_tr[i][1][1] * correspondence_input_tr[i][1][1]);
    }
    
    rms1 = rms1/4.0;
    rms2 = rms2/4.0;
    
    //  Generating the scaling matrices
    double scale_factor1 = sqrt(2)/rms1;
    double scale_factor2 = sqrt(2)/rms2;
    
    sc1(0,0) = scale_factor1;
    sc1(1,0) = 0.0;
    sc1(2,0) = 0;
    sc1(0,1) = 0.0;
    sc1(1,1) = scale_factor1;
    sc1(2,1) = 0;
    sc1(0,2) = 0.0;
    sc1(1,2) = 0.0;
    sc1(2,2) = 1;
    
    sc2(0,0) = scale_factor2;
    sc2(1,0) = 0.0;
    sc2(2,0) = 0;
    sc2(0,1) = 0.0;
    sc2(1,1) = scale_factor2;
    sc2(2,1) = 0;
    sc2(0,2) = 0.0;
    sc2(1,2) = 0.0;
    sc2(2,2) = 1;
    
    //cout <<"2" << endl;
    
    //  Normalizing correspondences to average RMS distance sqrt(2)
    double correspondence_input_norm[4][2][2] = {0};
    
    for (int i = 0; i < 4; i++)
    {
        CImg<double> temp_res(1,3);
        CImg<double> temp(1,3);
        temp(0,0) = correspondence_input_tr[i][0][0];
        temp(0,1) = correspondence_input_tr[i][0][1];
        temp(0,2) = 1;
        temp_res = sc1 * temp;
        correspondence_input_norm[i][0][0] = temp_res(0,0);
        correspondence_input_norm[i][0][1] = temp_res(0,1);
        
        temp(0,0) = correspondence_input_tr[i][1][0];
        temp(0,1) = correspondence_input_tr[i][1][1];
        temp(0,2) = 1;
        temp_res = sc2 * temp;
        correspondence_input_norm[i][1][0] = temp_res(0,0);
        correspondence_input_norm[i][1][1] = temp_res(0,1);
    }
    
    //  Generating Transformation matrices and required inverse for
    //  denormalization
    //cout <<"3" << endl;
    CImg<double> T1 = sc1*tr1;
    CImg<double> T2 = sc2*tr2;
    //  T2i = inv(T2);
    
    //  Generation of the A matrix for computation of the hmatrix
    CImg<double> A(9,8);
    for (int i = 0; i < 4; i++)
    {
        double x1 = correspondence_input_norm[i][0][0];
        double y1 = correspondence_input_norm[i][0][1];
        double x1_prime = correspondence_input_norm[i][1][0];
        double y1_prime = correspondence_input_norm[i][1][1];
        A(0, i*2) = x1; A(1, i*2) = y1; A(2, i*2) = 1;A(3, i*2) = 0;A(4, i*2) = 0;A(5, i*2) = 0;A(6, i*2) = -x1_prime*x1;A(7, i*2) = -x1_prime*y1;A(8, i*2) = -x1_prime;
        A(0, i*2 + 1) = 0; A(1, i*2+1) = 0; A(2, i*2+1) = 0;A(3, i*2+1) = x1;A(4, i*2+1) = y1;A(5, i*2+1) = 1;A(6, i*2+1) = -y1_prime*x1;A(7, i*2+1) = -y1_prime*y1;A(8, i*2+1) = -y1_prime;
        
    }
    
    //A.print("test1",true);
    //  Solving the linear equation Ax = 0 using the SVD
    CImg<double> U, S, V;
    A.SVD(U, S, V, true);
    //cout <<"4" << endl;
    //cout << "V.rows = " << V.height() << " V.cols = " << V.width() << endl;
    //CImg<double> H(3,3);
    H(0,0) = -V(8,0);
    H(1,0) = -V(8,1);
    H(2,0) = -V(8,2);
    
    H(0,1) = -V(8,3);
    H(1,1) = -V(8,4);
    H(2,1) = -V(8,5);
    
    H(0,2) = -V(8,6);
    H(1,2) = -V(8,7);
    H(2,2) = -V(8,8);
    //cout <<"5" << endl;
    //H.print("test1", true);
    //  Denormalizing the hmatrix_norm
    H = T2.invert(false) * H;
    H = H*T1;
    
}

double myfun( double *par, int x, void *fdata)
{
    //%MYFUN LevMar cost calculator
    CImg<double> hmatrixtemp(3,3);
    hmatrixtemp(0,0) = par[0];
    hmatrixtemp(1,0) = par[1];
    hmatrixtemp(2,0) = par[2];
    hmatrixtemp(0,1) = par[3];
    hmatrixtemp(1,1) = par[4];
    hmatrixtemp(2,1) = par[5];
    hmatrixtemp(0,2) = par[6];
    hmatrixtemp(1,2) = par[7];
    hmatrixtemp(2,2) = par[8];
    
    CImg<double> *points = (CImg<double> *)fdata;
    CImg<double> xdata = points[x];
    CImg<double> tmp1(1,3);
    tmp1(0,0) = xdata(0,0);
    tmp1(0,1) = xdata(1,0);
    tmp1(0,2) = 1;
    CImg<double> temp1 = hmatrixtemp*tmp1;
    temp1 = temp1/temp1(2,0);
    
    CImg<double> tmp2(1,3);
    tmp2(0,0) = xdata(0,1);
    tmp2(0,1) = xdata(1,1);
    tmp2(0,2) = 1;
    CImg<double> temp2 = hmatrixtemp.invert(false) * tmp2;
    temp2 = temp2/temp2(2,0);
    
    double error1 = abs(tmp2 - temp1).sum();
    double error2 = abs(tmp1 - temp2).sum();
    
    return error1 + error2;
}

void gradient(double *g, double *par, int x, void *fdata)
{
    for(int i = 0; i < 9; i++)
    {
        g[i] =1;
    }
}

void autoHomographies( const vector<CImg<double> > correspondence_input, CImg <double> &H )
{
    //AUTOHOMEST2D This function computes the automatic homography between two
    //perspective images using RANSAC and Levenberg Marquardt optimization.
    
    int corr_req = 4;
    int corr_num = correspondence_input.size();
    int max_no_inliers = 0;
    double p = 0.99;
    int iter_no = 0;
    CImg<double> max_inlier_detail;
    CImg<double> hmatrixauto;
    
    vector<int> indexes;
    
    for(int i = 0; i < correspondence_input.size(); i++)
    {
        indexes.push_back(i);
    }
    
    while(1)
    {
        iter_no = iter_no + 1;
        //cout << "iter_no = " <<iter_no <<endl;
        CImg<double>inlier_detail(corr_num,1);
        inlier_detail.fill(0);
        
        // Randomly shuffle the points
        
        random_shuffle ( indexes.begin(), indexes.end() );
        
        //Choosing the correspondences
        
        double selected_correspondences[4][2][2] = {0};
        
        for(int i=0; i < corr_req; i++)
        {
            CImg<double> correspondence = correspondence_input[indexes[i]];
            selected_correspondences[i][0][0] = correspondence(0,0);
            selected_correspondences[i][0][1] = correspondence(1,0);
            selected_correspondences[i][1][0] = correspondence(0,1);
            selected_correspondences[i][1][1] = correspondence(1,1);
        }
        
        //Computing the temporary H matrix
        CImg <double> hmatrixtemp(3,3);
        homographies(selected_correspondences, hmatrixtemp);
        //hmatrixtemp.print("hmatrixtemp", true);
        CImg<double> hmatrixtempi = hmatrixtemp;
        hmatrixtempi.invert(false);
        //Calculating the symmetric cost and checking for inliers
        
        int no_inliers = 0;
        //cout << "1" << endl;
        for(int i = 0; i < corr_num; i++)
        {
            //CImg<double> ci = correspondence_input[i];
            
            CImg<double> tmp1(1,3);
            tmp1(0,0) = correspondence_input[i](0,0);
            tmp1(0,1) = correspondence_input[i](1,0);
            tmp1(0,2) = 1;
            CImg<double> temp1 = hmatrixtemp*tmp1;
            //hmatrixtemp.print("hmatrixtemp", true);
            temp1 = temp1/temp1(0,2);
            
            //temp1.print("temp1", true);
            
            CImg<double> tmp2(1,3);
            tmp2(0,0) = correspondence_input[i](0,1);
            tmp2(0,1) = correspondence_input[i](1,1);
            tmp2(0,2) = 1;
            CImg<double> temp2 = hmatrixtempi*tmp2;
            //temp2.print("temp2", true);
            temp2 = temp2/temp2(0,2);
            
            //temp2.print("temp2", true);
            
            double error1 = abs(tmp2 - temp1).sum();
            double error2 = abs(tmp1 - temp2).sum();
            
            double error = error1 + error2;
            
            //cout << "error = " << error << endl;
            
            if(error < 20)
            {
                no_inliers = no_inliers + 1;
                inlier_detail(i,0) = 1;
            }
            
        }
        
        //cout << "no_inliers = "<<no_inliers << endl;
        
        if(no_inliers > max_no_inliers)
        {
            max_no_inliers = no_inliers;
            max_inlier_detail = inlier_detail;
            hmatrixauto = hmatrixtemp;
        }
        
        //Dynamic update of the no of iterations, N
        
        double w = (double)max_no_inliers/(double)corr_num;
        double N = abs((log(1-p))/(log(1-w*w*w*w)));
        
        //cout <<"w = " << w << endl;
        //cout << "N = " << N << endl;
        if(N < iter_no)
            break;
        
    }
    
    //H = hmatrixauto;
    //hmatrixauto.print("hmatrixauto", true);
    
    // Isolating inlier correspondences
    vector <CImg<double> > inlier_corr;
    int k = 0;
    for(int i=0; i < corr_num; i++)
    {
        if(max_inlier_detail(i,0) == 1)
        {
            inlier_corr.push_back(correspondence_input[i]);// correspondence_input(i,:,:);
            k = k+1;
        }
    }
    
    // LevMar Optimization
    double hvector[9] = {hmatrixauto(0,0), hmatrixauto(1,0), hmatrixauto(2,0), hmatrixauto(0,1), hmatrixauto(1,1), hmatrixauto(2,1),hmatrixauto(0,2),hmatrixauto(1,2),hmatrixauto(2,2)};
    
    double *ydata = (double *) malloc(max_no_inliers * sizeof(double));
    memset(ydata, 0, max_no_inliers * sizeof(*ydata));
    
    int n_iterations;
    LMstat lmstat;
    levmarq_init(&lmstat);
    int N_PARAMS = 9;
    n_iterations = levmarq(N_PARAMS, hvector, max_no_inliers, ydata, NULL, &myfun, &gradient, &inlier_corr[0], &lmstat);
    
    // Reassigning the H matrix
    H(0,0) = hvector[0]/hvector[8];
    H(1,0) = hvector[1]/hvector[8];
    H(2,0) = hvector[2]/hvector[8];
    H(0,1) = hvector[3]/hvector[8];
    H(1,1) = hvector[4]/hvector[8];
    H(2,1) = hvector[5]/hvector[8];
    H(0,2) = hvector[6]/hvector[8];
    H(1,2) = hvector[7]/hvector[8];
    H(2,2) = hvector[8]/hvector[8];
    
    //H = hmatrixauto/hmatrixauto(2,2);
}

double calcDist(vector<float> p1, vector<float> p2)
{
    double dist = 0;
    
    for(int i = 0; i < p1.size(); i++)
    {
        dist += (p1[i]-p2[i]) * (p1[i] - p2[i]);
    }
    
    return sqrt(dist);
}

bool myfunction (int i,int j) { return (i<j); }

int getMin(vector<double> v)
{
    int minIndex = 0;
    for(int i=0;i<v.size();i++)
    {
        if(v[minIndex]>v[i])
            minIndex = i;
    }
    return minIndex;
}

int getSecondMin(vector<double> v,int minIndex)
{
    int secondMinIndex = -1;
    for(int i=0;i<v.size();i++)
    {
        if(minIndex!=i)
        {
            secondMinIndex = i;
            break;
        }
    }
    for(int i=0;i<v.size();i++)
    {
        if(minIndex==i)
            continue;
        if(v[secondMinIndex]>v[i])
            secondMinIndex = i;
    }
    return secondMinIndex;
}

vector<CImg<double> > matching( CImg<double> I1, CImg<double>I2 )
{
    //SIFT matching
    double threshold = 10.0;
    
    // Get the Key Points
    CImg<double> gray1 = I1.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts1 = Sift::compute_sift(gray1);
    
    CImg<double> gray2 = I2.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts2 = Sift::compute_sift(gray2);
    int minimum,secondMinimum;
    double minDist,secondMinDist,ratio;
    vector<CImg<double> > points;
    int count = 0;
    for(int i = 0; i < Ipts1.size(); i++)
    {
        vector<double> distances;
        for(int j = 0; j < Ipts2.size(); j++)
        {
            distances.push_back(calcDist(Ipts1[i].descriptor,Ipts2[j].descriptor));
        }
        count++;
        minimum = getMin(distances);
        minDist = distances[minimum];
        secondMinimum = getSecondMin(distances,minimum);
        secondMinDist = distances[secondMinimum];
        ratio = secondMinDist/minDist;
        //cout<<ratio<<" ";
        if(ratio > threshold) //take this pair as match
        {
            CImg<double> img(2,2);
            img(0,0) = Ipts1[i].col;
            img(1,0) = Ipts1[i].row;
            
            img(0,1) = Ipts2[minimum].col;
            img(1,1) = Ipts2[minimum].row;
            
            points.push_back(img);
        }
    }
    cout<<"Number of pairs considered: "<<count<<endl;
    cout<<"Number of pairs taken as matches: "<<points.size()<<endl;
    return points;
}

vector<CImg<double> > readPoints()
{
    string path = "files/";
    
    string   str;
    DIR *pDIR;
    struct dirent *entry;
    vector<CImg<double> > points;
    if( pDIR=opendir(path.c_str()) )
    {
        while((entry = readdir(pDIR)) != NULL)
        {
            if( strcmp(entry->d_name, ".") != 0 && strcmp(entry->d_name, "..") != 0 )
            {
                //cout << "file name = " << entry->d_name << endl;
                string filename = path + string(entry->d_name);
                ifstream inn(filename.c_str());
                double x1,x2,y1,y2;
                inn >> x1;
                inn >> x2;
                inn >> y1;
                inn >> y2;
                //cout << x1 <<" " << x2 << " " << y1 << " " << y2 << endl;
                CImg<double> im(2,2);
                im(0,0) = x1;
                im(1,0) = y1;
                im(0,1) = x2;
                im(1,1) = y2;
                
                points.push_back(im);
                inn.close();
            }
        }
        closedir(pDIR);
    }
    return points;
}


// Nayem: part2
void autoHomographies_part2( const vector<CImg<double> > correspondence_input, CImg <double> &H, vector <CImg<double> > &inlier_correct)
{
    //AUTOHOMEST2D This function computes the automatic homography between two
    //perspective images using RANSAC and Levenberg Marquardt optimization.
    
    int corr_req = 4;
    int corr_num = correspondence_input.size();
    int max_no_inliers = 0;
    double p = 0.99;
    int iter_no = 0;
    CImg<double> max_inlier_detail;
    CImg<double> hmatrixauto;
    
    vector<int> indexes;
    
    for(int i = 0; i < correspondence_input.size(); i++)
    {
        indexes.push_back(i);
    }
    
    while(1)
    {
        iter_no = iter_no + 1;
        //cout << "iter_no = " <<iter_no <<endl;
        CImg<double>inlier_detail(corr_num,1);
        inlier_detail.fill(0);
        
        // Randomly shuffle the points
        
        random_shuffle ( indexes.begin(), indexes.end() );
        
        //Choosing the correspondences
        
        double selected_correspondences[4][2][2] = {0};
        
        for(int i=0; i < corr_req; i++)
        {
            CImg<double> correspondence = correspondence_input[indexes[i]];
            selected_correspondences[i][0][0] = correspondence(0,0);
            selected_correspondences[i][0][1] = correspondence(1,0);
            selected_correspondences[i][1][0] = correspondence(0,1);
            selected_correspondences[i][1][1] = correspondence(1,1);
        }
        
        //Computing the temporary H matrix
        CImg <double> hmatrixtemp(3,3);
        homographies(selected_correspondences, hmatrixtemp);
        //hmatrixtemp.print("hmatrixtemp", true);
        CImg<double> hmatrixtempi = hmatrixtemp;
        hmatrixtempi.invert(false);
        //Calculating the symmetric cost and checking for inliers
        
        int no_inliers = 0;
        //cout << "1" << endl;
        for(int i = 0; i < corr_num; i++)
        {
            //CImg<double> ci = correspondence_input[i];
            
            CImg<double> tmp1(1,3);
            tmp1(0,0) = correspondence_input[i](0,0);
            tmp1(0,1) = correspondence_input[i](1,0);
            tmp1(0,2) = 1;
            CImg<double> temp1 = hmatrixtemp*tmp1;
            //hmatrixtemp.print("hmatrixtemp", true);
            temp1 = temp1/temp1(0,2);
            
            //temp1.print("temp1", true);
            
            CImg<double> tmp2(1,3);
            tmp2(0,0) = correspondence_input[i](0,1);
            tmp2(0,1) = correspondence_input[i](1,1);
            tmp2(0,2) = 1;
            CImg<double> temp2 = hmatrixtempi*tmp2;
            //temp2.print("temp2", true);
            temp2 = temp2/temp2(0,2);
            
            //temp2.print("temp2", true);
            
            double error1 = abs(tmp2 - temp1).sum();
            double error2 = abs(tmp1 - temp2).sum();
            
            double error = error1 + error2;
            
            //cout << "error = " << error << endl;
            
            if(error < 20)
            {
                no_inliers = no_inliers + 1;
                inlier_detail(i,0) = 1;
            }
            
        }
        
        //cout << "no_inliers = "<<no_inliers << endl;
        
        if(no_inliers > max_no_inliers)
        {
            max_no_inliers = no_inliers;
            max_inlier_detail = inlier_detail;
            hmatrixauto = hmatrixtemp;
        }
        
        //Dynamic update of the no of iterations, N
        
        double w = (double)max_no_inliers/(double)corr_num;
        double N = abs((log(1-p))/(log(1-w*w*w*w)));
        
        //cout <<"w = " << w << endl;
        //cout << "N = " << N << endl;
        if(N < iter_no)
            break;
        
    }
    
    //H = hmatrixauto;
    //hmatrixauto.print("hmatrixauto", true);
    
    // Isolating inlier correspondences
    // inlier_corr -> is needed for part2
    /*vector <CImg<double> > inlier_corr;*/
    int k = 0;
    for(int i=0; i < corr_num; i++)
    {
        if(max_inlier_detail(i,0) == 1)
        {
            inlier_correct.push_back(correspondence_input[i]);// correspondence_input(i,:,:);
            k = k+1;
//            cout<<"img1:("<<correspondence_input[i](0,0)<<","<<correspondence_input[i](1,0)<<") - ";
//            cout<<"img2:("<<correspondence_input[i](0,1)<<","<<correspondence_input[i](1,1)<<")\n";
        }
    }
    
    vector <CImg<double> > inlier_corr = inlier_correct;
    // LevMar Optimization
    double hvector[9] = {hmatrixauto(0,0), hmatrixauto(1,0), hmatrixauto(2,0), hmatrixauto(0,1), hmatrixauto(1,1), hmatrixauto(2,1),hmatrixauto(0,2),hmatrixauto(1,2),hmatrixauto(2,2)};
    
    double *ydata = (double *) malloc(max_no_inliers * sizeof(double));
    memset(ydata, 0, max_no_inliers * sizeof(*ydata));
    
    int n_iterations;
    LMstat lmstat;
    levmarq_init(&lmstat);
    int N_PARAMS = 9;
    n_iterations = levmarq(N_PARAMS, hvector, max_no_inliers, ydata, NULL, &myfun, &gradient, &inlier_corr[0], &lmstat);
    
    // Reassigning the H matrix
    H(0,0) = hvector[0]/hvector[8];
    H(1,0) = hvector[1]/hvector[8];
    H(2,0) = hvector[2]/hvector[8];
    H(0,1) = hvector[3]/hvector[8];
    H(1,1) = hvector[4]/hvector[8];
    H(2,1) = hvector[5]/hvector[8];
    H(0,2) = hvector[6]/hvector[8];
    H(1,2) = hvector[7]/hvector[8];
    H(2,2) = hvector[8]/hvector[8];
    
    //H = hmatrixauto/hmatrixauto(2,2);
}



int main(int argc, char **argv)
{
    srand ( unsigned ( std::time(0) ) );
    try {
        
        if(argc < 2)
        {
            cout << "Insufficent number of arguments; correct usage:" << endl;
            cout << "    a2 part_id ..." << endl;
            return -1;
        }
        
        string part = argv[1];
        string inputFile = argv[2];
        
        if(part == "part1")
        {
            // This is just a bit of sample code to get you started, to
            // show how to use the SIFT library.
            
            //CImg<double> input_image(inputFile.c_str());
            CImg<double> input_image1("lincoln.png");
            CImg<double> input_image2("wrap.png");
            vector<CImg<double> > cv = matching(input_image1, input_image2);
            
            CImg<double> gray1 = input_image1.get_RGBtoHSI().get_channel(2);
            vector<SiftDescriptor> descriptors = Sift::compute_sift(gray1);
            
            for(int i=0; i<descriptors.size(); i++)
            {
                /*cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
                 for(int l=0; l<128; l++)
                 cout << descriptors[i].descriptor[l] << "," ;
                 cout << ")" << endl;*/
                
                for(int j=0; j<5; j++)
                {
                    for(int k=0; k<5; k++)
                    {
                        if(j==2 || k==2)
                        {
                            for(int p=0; p<3; p++)
                            {
                                int col = descriptors[i].col+k-1;
                                int row = descriptors[i].row+j-1;
                                if(col >= 0 && col < input_image1.width() && row >= 0 && row < input_image1.height())
                                {
                                    input_image1(col, row, 0, p) = 0;
                                }
                            }
                        }
                    }
                }
            }
            input_image1.get_normalize(0,255).save("sift.png");
        }
        else if(part == "part2")
        {
            // do something here!
            string inputFile_2 = argv[3];
            
            CImg<double> input_image_1(inputFile.c_str());
            CImg<double> input_image_2(inputFile_2.c_str());
            
            CImg <double> H(3,3);
            vector<CImg<double> > cv;
            
            cv = matching(input_image_1, input_image_2);
            vector <CImg<double> > inlier_corr;
            autoHomographies_part2(cv,H, inlier_corr);
            cout<<"autoHomographies_part2 done!"<<inlier_corr.size()<<endl;
            for(int i=0; i<inlier_corr.size(); i++)
            {
                CImg<double> corresp_point =inlier_corr[i];
                for(int j=0; j<5; j++)
                    for(int k=0; k<5; k++)
                        if(j==2 || k==2)
                            for(int p=0; p<3; p++){
                                int col = corresp_point(1,0)+k-1;
                                int row = corresp_point(0,0)+j-1;
                                if(col >= 0 && col < input_image_1.width() && row >= 0 && row < input_image_1.height())
                                    input_image_1(row, col, 0, p) = 0;
                                
                                cout<<"img1:("<<row<<","<<col<<") - ";
                                col = corresp_point(1,1)+k-1;
                                row = corresp_point(0,1)+j-1;
                                if(col >= 0 && col < input_image_2.width() && row >= 0 && row < input_image_2.height())
                                    input_image_2(row, col, 0, p) = 0;
                                
                                cout<<"img2:("<<row<<","<<col<<")\n";
                            }
            }
            
            input_image_1.get_normalize(0,255).save("part2_im1.png");
            input_image_2.get_normalize(0,255).save("part2_im2.png");
            
        }
        else if(part == "part3")
        {
            double trans_matrix[3][3] = {{1.1247,-0.3147,222.9409}, {0.1088, 0.6851, -19.9247}, {0.0003, -0.0006, 1.0828}};
            CImg<double> input_image(inputFile.c_str());
            CImg<double> out(input_image.width(), input_image.height(), 1, 3);
            out.fill(0);
            
            warp_image(input_image, out, trans_matrix);
            out.get_normalize(0,255).save("wrap.png");
            
            CImg <double> H(3,3);
            vector<CImg<double> > cv;
            /*
             double C[4][2][2] = {{{24.0079,  257.6604}, {25.5874,  328.1342}},
             {{440.5009 ,  17.8450}, {710.2346,   17.9124}},
             {{573.9990,  116.5888}, {479.3390,  152.2850}},
             {{13.8227,  330.3478}, {13.8710,  357.5412}}};
             
             //homest2d(C, H);
             
             for(int i = 0; i < 4; i++)
             {
            	CImg<double> c1(2,2);
            	c1(0,0) = C[i][0][0];
            	c1(0,1) = C[i][1][0];
            	c1(1,0) = C[i][0][1];
            	c1(1,1) = C[i][1][1];
            	cv.push_back(c1);
             }
             */
            CImg<double> im1("lincoln.png");
            CImg<double> im2("wrap.png");
            cv = matching(im1, im2);
            //cv = readPoints();
            autoHomographies(cv,H);
            H.print("test", true);
            double mat[3][3];
            H.invert(false);
            for(int i = 0; i < 3; i++)
            {
                for(int j = 0; j < 3; j++)
                {
                    mat[i][j]= H(j,i);
                }
            }
            warp_image(input_image, out, mat);
            out.get_normalize(0,255).save("warp1.png");
            
        }
        else
            throw std::string("unknown part!");
        
        // feel free to add more conditions for other parts (e.g. more specific)
        //  parts, for debugging, etc.
    }
    catch(const string &err) {
        cerr << "Error: " << err << endl;
    }
}







