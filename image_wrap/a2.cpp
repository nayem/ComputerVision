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
#include <time.h>

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

// /**********************************************/
// Global Constants
const int MODE_WIN_ALL_TRIAL = 1;
const int MODE_WIN_MAX_TRIAL = 2;
// /**********************************************/


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

vector<CImg<double> > matching( CImg<double> I1, CImg<double>I2 )
{
    //SIFT matching
    double thresh = 200;
    // Get the Key Points
    CImg<double> gray1 = I1.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts1 = Sift::compute_sift(gray1);
    
    CImg<double> gray2 = I2.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts2 = Sift::compute_sift(gray2);
    vector<int> correspondingIndex;
    vector<double> correspondingDist;
    for(int i = 0; i < Ipts1.size(); i++)
    {
        vector<float> pt1 = Ipts1[i].descriptor;
        double minDist = 9999999;
        int minIndex = 0;
        for(int j = 0; j < Ipts2.size(); j++)
        {
            vector<float> pt2 = Ipts2[j].descriptor;
            double dist = calcDist(pt1,pt2);
            if(minDist > dist)
            {
                minDist = dist;
                minIndex = j;
            }
            //			cout << "dist = " <<dist <<endl;
        }
        
        //		cout << "min index = " << minIndex;
        //		cout << ", minDist = " << minDist << endl;
        correspondingDist.push_back(minDist);
        correspondingIndex.push_back(minIndex);
    }
    
    //cout << "correspondingDist.size() = " <<correspondingDist.size() <<endl;
    vector<pair<double,int> > inds;
    for(int i = 0; i < correspondingDist.size(); i++)
    {
        inds.push_back(make_pair(correspondingDist[i],i));
    }
    
    sort( inds.begin(), inds.end());
    vector<CImg<double> > points;
    
    for(int i = 0; i < inds.size(); i++)
    {
        //cout << "dist = " << inds[i].first << endl;
        if(inds[i].first <= thresh)
        {
            CImg<double> img(2,2);
            int j = inds[i].second;
            img(0,0) = Ipts1[j].col;
            img(1,0) = Ipts1[j].row;
            
            int k = correspondingIndex[j];
            img(0,1) = Ipts2[k].col;
            img(1,1) = Ipts2[k].row;
            
            points.push_back(img);
        }
    }
    cout << "points match() = " <<points.size() <<endl;
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
void autoHomographies_part2( const vector<CImg<double> > correspondence_input, CImg <double> &H, vector <CImg<double> > &inlier_corr)
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
        //        cout << "iter_no = " <<iter_no <<endl;
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
        if(N < iter_no || iter_no>10000)
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
            inlier_corr.push_back(correspondence_input[i]);// correspondence_input(i,:,:);
            k = k+1;
            //            cout<<"img1:("<<correspondence_input[i](0,0)<<","<<correspondence_input[i](1,0)<<") - ";
            //            cout<<"img2:("<<correspondence_input[i](0,1)<<","<<correspondence_input[i](1,1)<<")\n";
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


void getVectorOfCimg(CImg<double> cimg, vector<double> &vec){
    
    for (int i=0; i<cimg.size(); ++i)
        vec.push_back(cimg(i));
}

vector<float> getColumnVectorOfCimg(CImg<double> cimg, int col){
    vector<float> vec;
    for (int r=0; r<cimg.height(); ++r){
        vec.push_back(cimg(col, r));
    }
    return vec;
}

vector<float> getRowVectorOfCimg(CImg<double> cimg, int row){
    vector<float> vec;
    for (int c=0; c<cimg.width(); ++c){
        vec.push_back(cimg(c, row));
    }
    return vec;
}

CImg<double> getSummarizeDescriptor(CImg<double> x_random, vector<SiftDescriptor> Ipts, double W){
    CImg<double> v_image(Ipts.size(), 128);
    for (int i=0; i<Ipts.size(); i++) {    //col
        vector<float> pt = Ipts[i].descriptor;
        
        for (int j=0; j<pt.size(); j++) { // row
            v_image(i, j) = pt[j];
        }
    }
    //    cimg_forXY(v_image,x,y){
    //        cout<<v_image(x,y)<<",";
    //    }
    //    cout<<"[v_image]->c:"<<v_image.width()<<",r:"<<v_image.height()<<endl;
    
    CImg<double> f_image = (x_random*v_image)/W;
    //    cimg_forXY(f_image, x, y){
    //        cout<<f_image(x,y)<<",";
    //    }
    //    cout<<"[f_image]->c:"<<f_image.width()<<",r:"<<f_image.height()<<endl;
    
    return f_image;
}

vector<CImg<double> > matching_quantized( CImg<double> I1, CImg<double>I2,int K=10, double W=1.0, int NUM_TRIAL=10, int MODE = MODE_WIN_MAX_TRIAL )
{
    const int MODE_WIN_ALL_TRIAL = 1;
    const int MODE_WIN_MAX_TRIAL = 2;
    
    //SIFT matching
    double thresh = 150;
    // Get the Key Points
    CImg<double> gray1 = I1.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts1 = Sift::compute_sift(gray1);
    
    CImg<double> gray2 = I2.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts2 = Sift::compute_sift(gray2);
    
    //    vector<int> correspondingIndex;
    //    vector<double> correspondingDist;
    vector<vector<int> > correspondingIndexTrial;
    vector<vector<double> > correspondingDistTrial;
    
    vector<vector<pair<int,double> > > indicesForTrial(Ipts1.size());
    
    
    for (int trial=0; trial<NUM_TRIAL; ++trial) {
        
        // Random
        CImg<double> x_random(128, K);
        cimg_forXY(x_random, x, y){
            x_random(x,y) = ((double) rand() / RAND_MAX);
        }
        //        cout<<"Trial:"<<trial<<endl;
        //        cout<<"[x_random]->c:"<<x_random.width()<<",r:"<<x_random.height()<<endl;
        
        // Summarize SIFT Descriptor
        CImg<double> f_image1 = getSummarizeDescriptor(x_random,Ipts1, W);
        CImg<double> f_image2 = getSummarizeDescriptor(x_random,Ipts2, W);
        
        //        cout<<"[Ipts1]->c:"<<Ipts1.size()<<endl;
        //        cout<<"[f_image1]->c:"<<f_image1.width()<<",r:"<<f_image1.height()<<endl;
        
        vector<int> correspondingIndex;
        vector<double> correspondingDist;
        
        vector< vector<float> >columnImage1;
        vector< vector<float> >columnImage2;
        
        for(int i = 0; i < f_image1.width(); i++)
            columnImage1.push_back( getColumnVectorOfCimg(f_image1,i) );
        
        for(int j = 0; j < f_image2.width(); j++)
            columnImage2.push_back( getColumnVectorOfCimg(f_image2,j) );
        
        
        for(int i = 0; i < f_image1.width(); i++)
        {
            //vector<float> pt1 = getColumnVectorOfCimg(f_image1,i);
            vector<float> pt1 = columnImage1[i];
            //        cout<<"pt1:";
            //        for (int q = 0; q<pt1.size(); q++)
            //            cout<<pt1[q]<<",";
            //        cout<<endl;
            
            double minDist = 9999999;
            int minIndex = 0;
            
            for(int j = 0; j < f_image2.width(); j++)
            {
                //            vector<float> pt2 = getColumnVectorOfCimg(f_image2,j);
                vector<float> pt2 = columnImage2[j];
                //            cout<<"pt2:";
                //            for (int q = 0; q<pt2.size(); q++)
                //                cout<<pt2[q]<<",";
                //            cout<<endl;
                
                //            cout<<"i:"<<i<<", j:"<<j<<endl;
                //            cout<<"pt1.size():"<<pt1.size()<<", pt2.size():"<<pt2.size()<<endl;
                
                double dist = calcDist(pt1,pt2);
                if(minDist > dist)
                {
                    minDist = dist;
                    minIndex = j;
                }
                //            cout << "dist = " <<dist <<endl;
            }
            
            //            cout <<"[Trial:"<<trial<<"]"<<"[i:"<<i<<"]"<< "min index = " << minIndex;
            //            cout << ", minDist = " << minDist << endl;
            
            switch (MODE) {
                case MODE_WIN_ALL_TRIAL:
                    if (trial !=0){
                        if (correspondingIndexTrial[trial-1][i]!=minIndex)
                            minDist=9999999;
                        else
                            minDist = (correspondingDistTrial[trial-1][i]+minDist)/2.0;
                    }
                    correspondingDist.push_back(minDist);
                    correspondingIndex.push_back(minIndex);
                    
                    correspondingDistTrial.push_back(correspondingDist);
                    correspondingIndexTrial.push_back(correspondingIndex);
                    
                    break;
                    
                case MODE_WIN_MAX_TRIAL:
                    indicesForTrial[i].push_back( make_pair(minIndex,minDist) );
                    break;
            }
            //correspondingDist.push_back(minDist);
            //correspondingIndex.push_back(minIndex);
            
        }
        
    }
    
    vector<pair<double,int> > inds;
    vector<int> correspondingIndex;
    vector<double> correspondingDist;
    
    cout<<"MODE:"<<MODE<<endl;
    
    switch (MODE) {
        case MODE_WIN_ALL_TRIAL:
            correspondingDist = correspondingDistTrial[correspondingDistTrial.size()-1];
            correspondingIndex = correspondingIndexTrial[correspondingIndexTrial.size()-1];
            for(int i = 0; i < correspondingDist.size(); i++){
                inds.push_back(make_pair(correspondingDist[i],i));
                //                cout <<"min index = " << i;
                //                cout << ", minDist = " << correspondingDist[i] << endl;
            }
            break;
            
        case MODE_WIN_MAX_TRIAL:
            cout<<"correspondingDistTrial.size():"<<correspondingDistTrial.size()<<endl;
            for(int i = 0; i < Ipts1.size(); i++){
                sort(indicesForTrial[i].begin(), indicesForTrial[i].end());
                int modeIndex = indicesForTrial[i][0].first, mdIndx=indicesForTrial[i][0].first;
                int occurance = 1, occr=1;
                double modeDistance=0, mdDist=indicesForTrial[i][0].second;
                
                //                cout<<"i:"<<i<<",mdIndx:"<<mdIndx<<endl;
                bool isAllSame=true;
                
                for (int t=1; t<NUM_TRIAL; ++t) {
                    if (mdIndx != indicesForTrial[i][t].first) {
                        isAllSame=false;
                        if (occurance<occr) {
                            modeIndex=mdIndx;
                            occurance=occr;
                            modeDistance=mdDist/occurance;
                            
                            mdIndx = indicesForTrial[i][t].first;
                            occr=0;
                            mdDist=0;
                        }
                    }
                    occr+=1;
                    mdDist+=indicesForTrial[i][t].second;
                }
                if (isAllSame) {
                    modeIndex=mdIndx;
                    occurance=occr;
                    modeDistance=mdDist/occurance;
                }
                inds.push_back(make_pair(modeDistance,i));
                correspondingDist.push_back(modeDistance);
                correspondingIndex.push_back(modeIndex);
                
                //                cout<<"i:"<<i<<",x:"<<Ipts1[i].col<<",y:"<<Ipts1[i].row;
                //                cout<<",modeIndex:"<<modeIndex<<",x:"<<Ipts2[modeIndex].col<<",y:"<<Ipts2[modeIndex].row;
                //                cout<<",modeDistance:"<<modeDistance<<endl;
            }
            
            
            
    }
    //    vector<pair<double,int> > inds;
    //    for(int i = 0; i < correspondingDist.size(); i++)
    //    {
    //        inds.push_back(make_pair(correspondingDist[i],i));
    //    }
    
    sort( inds.begin(), inds.end());
    vector<CImg<double> > points;
    cout<<inds.size()<<endl;
    cout<<correspondingIndex.size()<<endl;
    
    for(int i = 0; i < inds.size(); i++)
    {
        //cout << "dist = " << inds[i].first << endl;
        if(inds[i].first <= thresh)
        {
            CImg<double> img(2,2);
            int j = inds[i].second;
            img(0,0) = Ipts1[j].col;
            img(1,0) = Ipts1[j].row;
            
            int k = correspondingIndex[j];
            img(0,1) = Ipts2[k].col;
            img(1,1) = Ipts2[k].row;
            
            points.push_back(img);
        }
    }
    cout<<points.size()<<endl;
    
    return points;
}
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

vector<CImg<double> > getMatching( CImg<double> I1, CImg<double>I2,double threshold)
{
    //SIFT matching
    //double threshold = 5;
    
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
            
            //cout<<Ipts1[i].row<<" "<<Ipts1[i].col<<endl;
            
            img(0,1) = Ipts2[minimum].col;
            img(1,1) = Ipts2[minimum].row;
            
            //cout<<Ipts2[minimum].row<<" "<<Ipts2[minimum].col<<endl;
            
            points.push_back(img);
        }
    }
    cout<<"Number of pairs considered: "<<count<<endl;
    cout<<"Number of pairs taken as matches: "<<points.size()<<endl;
    return points;
}
void bhm_line(CImg<double>& combinedImage,int x1,int y1,int x2,int y2)
{
    int x,y,dx,dy,dx1,dy1,px,py,xe,ye,i;
    dx=x2-x1;
    dy=y2-y1;
    dx1=fabs(dx);
    dy1=fabs(dy);
    px=2*dy1-dx1;
    py=2*dx1-dy1;
    if(dy1<=dx1)
    {
        if(dx>=0)
        {
            x=x1;
            y=y1;
            xe=x2;
        }
        else
        {
            x=x2;
            y=y2;
            xe=x1;
        }
        combinedImage(y, x, 0, 0) = 0;
        combinedImage(y, x, 0, 1) = 0;
        combinedImage(y, x, 0, 2) = 0;
        
        //putpixel(x,y,c);
        for(i=0;x<xe;i++)
        {
            x=x+1;
            if(px<0)
            {
                px=px+2*dy1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    y=y+1;
                }
                else
                {
                    y=y-1;
                }
                px=px+2*(dy1-dx1);
            }
            //delay(0);
            //putpixel(x,y,c);
            combinedImage(y, x, 0, 0) = 0;
            combinedImage(y, x, 0, 1) = 0;
            combinedImage(y, x, 0, 2) = 0;
        }
    }
    else
    {
        if(dy>=0)
        {
            x=x1;
            y=y1;
            ye=y2;
        }
        else
        {
            x=x2;
            y=y2;
            ye=y1;
        }
        //putpixel(x,y,c);
        combinedImage(y, x, 0, 0) = 0;
        combinedImage(y, x, 0, 1) = 0;
        combinedImage(y, x, 0, 2) = 0;
        
        for(i=0;y<ye;i++)
        {
            y=y+1;
            if(py<=0)
            {
                py=py+2*dx1;
            }
            else
            {
                if((dx<0 && dy<0) || (dx>0 && dy>0))
                {
                    x=x+1;
                }
                else
                {
                    x=x-1;
                }
                py=py+2*(dx1-dy1);
            }
            //delay(0);
            //putpixel(x,y,c);
            combinedImage(y, x, 0, 0) = 0;
            combinedImage(y, x, 0, 1) = 0;
            combinedImage(y, x, 0, 2) = 0;
        }
    }
}
void drawLine(CImg<double>& combinedImage,int row1,int col1,int row2,int col2)
{
    int dx = row2-row1;
    int dy = col2-col1;
    if(dx==0)
        return;
    
    bhm_line(combinedImage,row1,col1,row2,col2);
    
    /*int y;
     for(int x=row1;x<=row2;x++)
     {
     y = floor(col1 + dy * (x - row1)/ dx);
     if(y >= 0 && y < combinedImage.width() && x >= 0 && x < combinedImage.height())
     {
     combinedImage(y, x, 0, 0) = 0;
     combinedImage(y, x, 0, 1) = 0;
     combinedImage(y, x, 0, 2) = 0;
     }
     }*/
}
void markMatchedPoints(vector<CImg<double> > cv,CImg<double> input_image1,CImg<double> input_image2,double threshold)
{
    CImg<double> combinedImage(input_image1.width()+input_image2.width(),max(input_image1.height(),input_image2.height()) ,1,3);
    combinedImage.fill(0);
    cout<<input_image1.height()<<" "<<input_image1.width()<<endl;
    cout<<input_image2.height()<<" "<<input_image2.width()<<endl;
    cout<<combinedImage.height()<<" "<<combinedImage.width()<<endl;
    for(int i=0;i<input_image1.height();i++)
    {
        for(int j=0;j<input_image1.width();j++)
        {
            for(int p=0;p<3;p++)
            {
                combinedImage(j,i,0,p) = input_image1(j,i,0,p);
                combinedImage(input_image1.width()+j,i,0,p) = input_image2(j,i,0,p);
            }
        }
    }
    
    for(int i=0;i<input_image2.height();i++)
    {
        for(int j=0;j<input_image2.width();j++)
        {
            for(int p=0;p<3;p++)
            {
                combinedImage(input_image1.width()+j,i,0,p) = input_image2(j,i,0,p);
            }
        }
    }
    int row1,col1,row2,col2;
    for(int i=0; i<cv.size(); i++)
    {
        //row1 = cv[i](1,0);
        //col1 = cv[i](0,0);
        //cout<<row1<<" "<<col1<<endl;
        //cout<<cv[i](1,0)<<" "<<cv[i](0,0)<<endl;
        for(int j=0; j<5; j++)
        {
            for(int k=0; k<5; k++)
            {
                if(j==2 || k==2)
                {
                    for(int p=0; p<3; p++)
                    {
                        int col = (int)floor(cv[i](0,0))+k-1;
                        int row = (int)floor(cv[i](1,0))+j-1;
                        //cout<<"Here is: "<<row<<' '<<col<<endl;
                        if(col >= 0 && col < combinedImage.width() && row >= 0 && row < combinedImage.height())
                        {
                            combinedImage(col, row, 0, p) = 0;
                        }
                    }
                }
            }
        }
        //row2 = cv[i](1,1);
        //col2 = cv[i](0,1);
        for(int j=0; j<5; j++)
        {
            for(int k=0; k<5; k++)
            {
                if(j==2 || k==2)
                {
                    for(int p=0; p<3; p++)
                    {
                        int col = input_image1.width()+(int)floor(cv[i](0,1))+k-1;
                        int row = (int)floor(cv[i](1,1))+j-1;
                        if(col >= 0 && col < combinedImage.width() && row >= 0 && row < combinedImage.height())
                        {
                            combinedImage(col, row, 0, p) = 0;
                        }
                    }
                }
            }
        }
        row1 = (int)floor(cv[i](1,0));
        col1 = (int)floor(cv[i](0,0));
        row2 = (int)floor(cv[i](1,1));
        col2 = input_image1.width()+(int)floor(cv[i](0,1));
        /*
         int dx = row2-row1;
         int dy = col2-col1;
         if(dx==0)
         continue;
         int y;
         for(int x=row1;x<=row2;x++)
         {
         y = floor(col1 + dy * (x - row1)/ dx);
         if(y >= 0 && y < combinedImage.width() && x >= 0 && x < combinedImage.height())
         {
         combinedImage(y, x, 0, 0) = 0;
         combinedImage(y, x, 0, 1) = 0;
         combinedImage(y, x, 0, 2) = 0;
         }
         }*/
        drawLine(combinedImage,row1,col1,row2,col2);
        
        //break;
    }
    
    char name[50];
    sprintf(name,"combined%.2lf.png",threshold);
    combinedImage.get_normalize(0,255);
    combinedImage.save(name);
    
}
struct less_than_key
{
    bool operator()(const std::pair<int,char*> &left, const std::pair<int,char*> &right)
    {
        return left.first > right.first;
    }
};
vector<char*> retrieveImage(string qI,int argc,char** argv,double threshold)
{
    vector<char*> bestMatches;
    vector<char*> database;
    string fileName;
    for(int i=3;i<argc;i++)
    {
        database.push_back(argv[i]);
    }
    
    CImg<double> queryImage(qI.c_str());
    vector<pair<int,char*> > matchesFound;
    for(int i=0;i<database.size();i++)
    {
        CImg<double> testImage(database[i]);
        matchesFound.push_back(make_pair((getMatching(queryImage,testImage,threshold)).size(),database[i]));
    }
    
    std::sort(matchesFound.begin(),matchesFound.end(),less_than_key());
    //cout<<"Sorted output: ";
    //for(int i=0;(i<10)&&(i<matchesFound.size());i++)
    for(int i=0;i<matchesFound.size();i++)
    {
        cout<<(i+1)<<" "<<matchesFound[i].first<<" "<<matchesFound[i].second<<endl;
        bestMatches.push_back(matchesFound[i].second);
    }
    return bestMatches;
}

vector<CImg<double> > getMatchingQuantized( CImg<double> I1, CImg<double>I2,double threshold, int K=10, double W=1.0, int NUM_TRIAL=10, int MODE = MODE_WIN_MAX_TRIAL)
{
    //SIFT matching
    //double threshold = 5;
    
    // Get the Key Points
    CImg<double> gray1 = I1.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts1 = Sift::compute_sift(gray1);
    
    CImg<double> gray2 = I2.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> Ipts2 = Sift::compute_sift(gray2);
    
    vector<vector<pair<int,double> > > indicesForTrial(Ipts1.size());
    
    for (int trial=0; trial<NUM_TRIAL; ++trial) {
        
        // Random
        CImg<double> x_random(128, K);
        cimg_forXY(x_random, x, y){
            x_random(x,y) = ((double) rand() / RAND_MAX);
        }
        // Summarize SIFT Descriptor
        CImg<double> f_image1 = getSummarizeDescriptor(x_random,Ipts1, W);
        CImg<double> f_image2 = getSummarizeDescriptor(x_random,Ipts2, W);
        
        vector< vector<float> >columnImage1;
        vector< vector<float> >columnImage2;
        
        for(int i = 0; i < f_image1.width(); i++)
            columnImage1.push_back( getColumnVectorOfCimg(f_image1,i) );
        
        for(int j = 0; j < f_image2.width(); j++)
            columnImage2.push_back( getColumnVectorOfCimg(f_image2,j) );
        
        int minimum,secondMinimum;
        double minDist,secondMinDist,ratio;
        vector<CImg<double> > points;
        int count = 0;
        
        for(int i = 0; i < f_image1.width(); i++)
        {
            vector<float> pt1 = columnImage1[i];
            vector<double> distances;
            for(int j = 0; j < f_image2.width(); j++)
            {
                vector<float> pt2 = columnImage2[j];
                //double dist = calcDist(pt1,pt2);
                distances.push_back(calcDist(pt1,pt2));
                //distances.push_back(calcDist(Ipts1[i].descriptor,Ipts2[j].descriptor));
            }
            count++;
            minimum = getMin(distances);
            minDist = distances[minimum];
            secondMinimum = getSecondMin(distances,minimum);
            secondMinDist = distances[secondMinimum];
            ratio = secondMinDist/minDist;
            //cout<<ratio<<" ";
            
            //            if(ratio > threshold) //take this pair as match
            //            {
            //                CImg<double> img(2,2);
            //                img(0,0) = Ipts1[i].col;
            //                img(1,0) = Ipts1[i].row;
            //
            //                //cout<<Ipts1[i].row<<" "<<Ipts1[i].col<<endl;
            //
            //                img(0,1) = Ipts2[minimum].col;
            //                img(1,1) = Ipts2[minimum].row;
            //
            //                //cout<<Ipts2[minimum].row<<" "<<Ipts2[minimum].col<<endl;
            //
            //                points.push_back(img);
            //            }
            switch (MODE) {
                case MODE_WIN_ALL_TRIAL:
                    //                    if (trial !=0){
                    //                        if (correspondingIndexTrial[trial-1][i]!=minIndex)
                    //                            minDist=9999999;
                    //                        else
                    //                            minDist = (correspondingDistTrial[trial-1][i]+minDist)/2.0;
                    //                    }
                    //                    correspondingDist.push_back(minDist);
                    //                    correspondingIndex.push_back(minIndex);
                    //
                    //                    correspondingDistTrial.push_back(correspondingDist);
                    //                    correspondingIndexTrial.push_back(correspondingIndex);
                    
                    break;
                    
                case MODE_WIN_MAX_TRIAL:
                    indicesForTrial[i].push_back( make_pair(minimum,ratio) );
                    break;
            }
        }
        //        cout<<"Number of pairs considered: "<<count<<endl;
        //        cout<<"Number of pairs taken as matches: "<<points.size()<<endl;
    }
    
    vector<pair<double,int> > inds;
    vector<int> correspondingIndex;
    vector<double> correspondingDist;
    
    cout<<indicesForTrial[0].size()<<endl;
    
    switch (MODE) {
        case MODE_WIN_ALL_TRIAL:
            break;
            
        case MODE_WIN_MAX_TRIAL:
            for(int i = 0; i < Ipts1.size(); i++){
                sort(indicesForTrial[i].begin(), indicesForTrial[i].end());
                int modeIndex = indicesForTrial[i][0].first, mdIndx=indicesForTrial[i][0].first;
                int occurance = 1, occr=1;
                double modeDistance=0, mdDist=indicesForTrial[i][0].second;
                
                bool isAllSame=true;
                
                for (int t=1; t<NUM_TRIAL; ++t) {
                    if (mdIndx != indicesForTrial[i][t].first) {
                        isAllSame=false;
                        if (occurance<occr) {
                            modeIndex=mdIndx;
                            occurance=occr;
                            modeDistance=mdDist/occurance;
                            
                            mdIndx = indicesForTrial[i][t].first;
                            occr=0;
                            mdDist=0;
                        }
                    }
                    occr+=1;
                    mdDist+=indicesForTrial[i][t].second;
                }
                if (isAllSame) {
                    modeIndex=mdIndx;
                    occurance=occr;
                    modeDistance=mdDist/occurance;
                }
                inds.push_back(make_pair(modeDistance,i));
                correspondingDist.push_back(modeDistance);
                correspondingIndex.push_back(modeIndex);
            }
            
    }
    
    vector<CImg<double> > points;
    for(int i = 0; i < inds.size(); i++)
    {
        if(inds[i].first <= threshold)
        {
            CImg<double> img(2,2);
            int j = inds[i].second;
            img(0,0) = Ipts1[j].col;
            img(1,0) = Ipts1[j].row;
            
            int k = correspondingIndex[j];
            img(0,1) = Ipts2[k].col;
            img(1,1) = Ipts2[k].row;
            
            points.push_back(img);
        }
    }
    cout<<"Point():"<<points.size()<<endl;
    
    return points;
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
        
        if(part == "part1")
        {
            string inputFile = argv[2];
            //CImg<double> input_image(inputFile.c_str());
            //CImg<double> input_image1("lincoln.png");
            //CImg<double> input_image2("wrap.png");
            double threshold = 10.0;
            //vector<CImg<double> > cv = getMatching(input_image1, input_image2,threshold);
            string queryImage = inputFile;
            vector<char*> bestMatches = retrieveImage(queryImage,argc,argv,threshold);
            for(int i=0;i<bestMatches.size();i++)
            {
                cout<<(i+1)<<": "<<bestMatches[i]<<endl;
            }
            
            //markMatchedPoints(cv,input_image1,input_image2,threshold);
            
            
            //CImg<double> gray1 = input_image1.get_RGBtoHSI().get_channel(2);
            //vector<SiftDescriptor> descriptors = Sift::compute_sift(gray1);
            
            /*for(int i=0; i<descriptors.size(); i++)
            	{
             //cout<<descriptors[i].row<<' '<<descriptors[i].col<<endl;
             cout << "Descriptor #" << i << ": x=" << descriptors[i].col << " y=" << descriptors[i].row << " descriptor=(";
             for(int l=0; l<128; l++)
             cout << descriptors[i].descriptor[l] << "," ;
             cout << ")" << endl;
             
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
             //cout<<row<<' '<<col<<endl;
             if(col >= 0 && col < input_image1.width() && row >= 0 && row < input_image1.height())
             {
             input_image1(col, row, 0, p) = 0;
             }
             }
             }
             }
             }
            	}
            	input_image1.get_normalize(0,255).save("sift.png");*/
        }
        else if(part == "part2")
        {
            time_t timer;
            double seconds;
            
            // for command: ./a2 part2 query.png img_1.png
            // Do RANSAC for these 2 image only and draw inliners
            if (argc == 4) {
                string inputFile = argv[2];
                // string part = argv[1];
                // string inputFile = argv[2];
                string inputFile_2 = argv[3];
                
                CImg<double> input_image_1(inputFile.c_str());
                CImg<double> input_image_2(inputFile_2.c_str());
                
                CImg <double> H(3,3);
                vector<CImg<double> > cv;
                
                //double threshold = atof(argv[4]);
                double threshold = 10;
                
                //cv = getMatchingQuantized(input_image_1,input_image_2,threshold);
                cv = getMatching(input_image_1,input_image_2,threshold);
                // Quantized
                //                int k_good=45;
                //                double w_good = 100;
                //                int trial_good=10;
                //                cv = getMatchingQuantized(input_image_1,input_image_2,threshold,k_good, w_good, trial_good, MODE_WIN_MAX_TRIAL);
                
                vector <CImg<double> > inlier_corr;
                autoHomographies_part2(cv,H, inlier_corr);
                
                // Draw the inliners
                for(int i=0; i<inlier_corr.size(); i++)
                {
                    CImg<double> corresp_point =inlier_corr[i];
                    for(int j=0; j<5; j++)
                        for(int k=0; k<5; k++)
                            if(j==2 || k==2)
                                for(int p=0; p<3; p++){
                                    
                                    int col = corresp_point(0,0)+k-1;
                                    int row = corresp_point(1,0)+j-1;
                                    if(col >= 0 && col < input_image_1.width() && row >= 0 && row < input_image_1.height())
                                        input_image_1( col,row, 0, p) = 0;
                                    
                                    col = corresp_point(0,1)+k-1;
                                    row = corresp_point(1,1)+j-1;
                                    if(col >= 0 && col < input_image_2.width() && row >= 0 && row < input_image_2.height())
                                        input_image_2(col, row,0, p) = 0;
                                }
                }
                markMatchedPoints(cv,input_image_1,input_image_2,threshold);
                input_image_1.get_normalize(0,255).save("part2_im1.png");
                input_image_2.get_normalize(0,255).save("part2_im2.png");
                // draw figure of connecting line
            }
            // for command: ./a2 part2
            else{
                
                // 2-1: RANSAC implemented
                // With out Quantized Projection
                
                // Read file
                ifstream file("a2_part1_imageNames.txt");
                string fileDirectory = "a2-images/part1_images/";
                string line;
                
                double threshold = atof(argv[4]);
                vector<string> fileNames;
                
                while (getline(file, line))
                    fileNames.push_back(fileDirectory+line);
                
                int NUM_PRECISION=1;
                
                timer = time(NULL);
                
                for (int i=0, img_rand; i<NUM_PRECISION; ++i) {
                    img_rand = (int)floor(((double)rand()/RAND_MAX)*100.0);
                    
                    CImg<double> input_image_1(fileNames[img_rand].c_str());
                    
                    vector<pair<int, int> > score_files;
                    for (int j=0; j<fileNames.size(); ++j) {
                        if (j==img_rand)
                            continue;
                        
                        cout<<img_rand<<"$$$$$$$$$ ["<<j<<"] $$$$$$$$$"<<fileNames[img_rand]<<endl;
                        CImg<double> input_image_2(fileNames[j].c_str());
                        
                        CImg <double> H(3,3);
                        vector<CImg<double> > cv;
                        vector <CImg<double> > inlier_corr;
                        
                        //cv = getMatchingQuantized(input_image_1,input_image_2,1.04);
                        cv = getMatching(input_image_1,input_image_2,threshold);
                        autoHomographies_part2(cv,H, inlier_corr);
                        
                        score_files.push_back(make_pair(inlier_corr.size(),j));
                    }
                    
                    // Top 10 files
                    sort(score_files.begin(),score_files.end());
                    double precision = 0;
                    for (int j=0; j<10; ++j){
                        cout<<i<<"->File:"<<fileNames[score_files[j].second]<<",score:"<<score_files[j].first<<endl;
                        
                        if((int)img_rand/10 == (int)score_files[j].second/10)
                            precision+=1;
                    }
                    precision = precision/10.0;
                    cout<<"Precision:"<<precision<<endl;
                    seconds = difftime(time(NULL), timer);
                    cout<<"Execution Time for Precision:"<<(seconds)<<endl;
                    timer = time(NULL);
                    
                }
                
                //2-2
                // With Quantized Projection
                /*///////////////////////////////////////////
                 ofstream resultOutputFile;
                 resultOutputFile.open ("result_part2b.csv");
                 resultOutputFile << "Writing this to a file.\n";
                 
                 const int MAX_K = 100;
                 const int MIN_K = 30;
                 const double MAX_W = 100;
                 const double MIN_W = 1.0;
                 const int MAX_TRIAL = 20;
                 const int MIN_TRIAL = 1;
                 
                 cout<<"Part2"<<endl;
                 for(int k=MIN_K; k<MAX_K; k+=5){
                 cout<<"Part2 inside K"<<endl;
                 for(double w=MIN_W; w<MAX_W; w+=10){
                 cout<<"Part2 inside W"<<endl;
                 for(int trial=MIN_TRIAL; trial<MAX_TRIAL; trial+=2){
                 
                 cout<<"Part2 inside T"<<endl;
                 int NUM_PRECISION=1;
                 
                 timer = time(NULL);
                 
                 for (int i=0, img_rand; i<NUM_PRECISION; ++i) {
                 img_rand = (int)floor(((double)rand()/RAND_MAX)*100.0);
                 
                 CImg<double> input_image_1(fileNames[img_rand].c_str());
                 
                 vector<pair<int, int> > score_files;
                 for (int j=0; j<fileNames.size(); ++j) {
                 cout<<img_rand<<"$$$$$$$$$ ["<<j<<"] $$$$$$$$$"<<fileNames[img_rand]<<endl;
                 CImg<double> input_image_2(fileNames[j].c_str());
                 
                 CImg <double> H(3,3);
                 vector<CImg<double> > cv;
                 vector <CImg<double> > inlier_corr;
                 
                 int k_good=45;
                 double w_good = 100;
                 int trial_good=10;
                 double threshold = 10;
                 
                 // Quantized Projection
                 //matching_quantized( CImg<double> I1, CImg<double>I2,int K=10, double W=1.0, int NUM_TRIAL=10, int MODE = MODE_WIN_MAX_TRIAL )
                 //MODE_WIN_ALL_TRIAL= take the always winner inlier
                 //MODE_WIN_MAX_TRIAL= take the maximum winner inlier
                 cv = getMatchingQuantized(input_image_1,input_image_2,threshold,k_good, w_good, trial_good, MODE_WIN_MAX_TRIAL);
                 autoHomographies_part2(cv,H, inlier_corr);
                 
                 score_files.push_back(make_pair(inlier_corr.size(),j));
                 }
                 
                 // Top 10 files
                 sort(score_files.begin(),score_files.end());
                 double precision = 0;
                 for (int j=0; j<10; ++j){
                 cout<<i<<"->File:"<<fileNames[score_files[j].second]<<",score:"<<score_files[j].first<<endl;
                 
                 if((int)img_rand/10 == (int)score_files[j].second/10)
                 precision+=1;
                 }
                 precision = precision/10.0;
                 cout<<"Precision:"<<precision<<endl;
                 seconds = difftime(time(NULL), timer);
                 cout<<"Execution Time for Precision:"<<(seconds)<<endl;
                 timer = time(NULL);
                 
                 // w|k|Trial|Accuracy|Time
                 resultOutputFile <<w<<","<<k<<","<<trial<<","<<seconds<<","<<precision<<"\n";
                 }
                 } // loop-trial
                 } // loop-w
                 } // loop-k
                 
                 resultOutputFile.close();
                 //////////////////////////////////////////*/
            }
            
        }
        else if(part == "part3")
        {
            /*3-1: Warp the image using given transformation matrix.*/
            /*We inverted the matrix and inverse warping algorithm using bilinear interpolation*/
            /*
             
             double trans_matrix[3][3] = {{1.1247,-0.3147,222.9409}, {0.1088, 0.6851, -19.9247}, {0.0003, -0.0006, 1.0828}};
             CImg<double> input_image(inputFile.c_str());
             cout<<"He"<<endl;
             CImg<double> out(input_image.width(), input_image.height(), 1, 3);
             out.fill(0);
             
             warp_image(input_image, out, trans_matrix);
             out.get_normalize(0,255).save("wrap.png");
             */
            
            //            CImg <double> H(3,3);
            //            vector<CImg<double> > cv;
            //            vector <CImg<double> > inlier_corr;
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
            string inputFile = argv[2];
            double threshold = 1.05;
            CImg<double> im1(inputFile.c_str());
            
            for (int k=3; k<argc; ++k) {
                
                string inputFile_2 = argv[k];
                
                CImg<double> im2(inputFile_2.c_str());
                
                CImg <double> H(3,3);
                vector<CImg<double> > cv;
                vector <CImg<double> > inlier_corr;
                
                //cv = getMatchingQuantized(im2,im1,threshold);
                cv = getMatching(im2,im1,threshold);
                autoHomographies_part2(cv,H, inlier_corr);
                //cv = readPoints();
                //autoHomographies(cv,H);
                
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
                
                CImg<double> out2(im2.width(), im2.height(), 1, 3);
                out2.fill(0);
                //warp_image(input_image, out, mat);
                warp_image(im2, out2, mat);
                char name[20];
                sprintf(name,"wrap_image_%d.png",k);
                //                (string("warp_img_")+string.itoa(i)+string(".png")).c_str()
                out2.get_normalize(0,255).save(name);
            }
            
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








