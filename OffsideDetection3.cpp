#include <stdio.h>
#include <cv.h>
#include "highgui.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#define PI 3.14159265
#define Root2 1.4142135623
#define Root3 1.7320508075
using namespace cv;
using namespace std;

Size matSize;
int origin[1000][2000][3];
int grass[1000][2000][3];
int gaussian[1000][2000][3];
double gradient[1000][2000][3];
uchar tshGradient[1000][2000][3];
int mark[1000][2000];
int qI[2000000],qJ[2000000];
Mat myMat;
vector<Vec4i> lines;
double m[1000];
double c[1000];
double M[1000];
double C[1000];
double m1, m2, d1, d2;
int fieldLineSize=0;
const double maxDifM = 0.2,;
const double maxDifC = 5;
vector<Vec4i> fieldLines;
vector< vector<Vec4i> > groupLines;
vector<Vec4i> horizLines;
vector<Vec4i> vertLines;
Point refPt;
  
void findPlayer()
{
     int I,J,i,j,k,r,p;
     int H=matSize.height;
     int W=matSize.width;
     int area[200000];
     int numArea=0;
     
     for(i=0;i<H;i++)
     for(j=0;j<W;j++)
     mark[i][j]=0;
     
     for(i=0;i<H;i++)
     for(j=0;j<W;j++)
     {
        if(tshGradient[i][j][0]==0&&mark[i][j]==0)
        {
           r=0;p=0;
           qI[r]=i;qJ[r]=j;r++;
           
           while(p<r)
           {
              I=qI[p];J=qJ[p];p++;
              if(J+1<W){if(tshGradient[I][J+1][0]==0&&mark[I][J+1]==0){qI[r]=I;qJ[r]=J+1;mark[qI[r]][qJ[r]]=1;r++;}}
              if(J-1>=0){if(tshGradient[I][J-1][0]==0&&mark[I][J-1]==0){qI[r]=I;qJ[r]=J-1;mark[qI[r]][qJ[r]]=1;r++;}}
              if(I+1<H){if(tshGradient[I+1][J][0]==0&&mark[I+1][J]==0){qI[r]=I+1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
              if(I-1>=0){if(tshGradient[I-1][J][0]==0&&mark[I-1][J]==0){qI[r]=I-1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
           }
           if(50<r && r<150)
           {
              //printf("%d ",r);
              for(k=0;k<r;k++)
              {tshGradient[qI[k]][qJ[k]][0]=128;tshGradient[qI[k]][qJ[k]][1]=128;tshGradient[qI[k]][qJ[k]][2]=128;}
           }
           area[numArea++]=r;
        }
     }    
     printf("\n");
}   

bool myfunction (Vec4i x,Vec4i y) { 
     double distX = sqrt( pow((double)x[0]-x[2],2) + pow((double)x[1]-x[3],2));
     double distY = sqrt( pow((double)y[0]-y[2],2) + pow((double)y[1]-y[3],2));
     return (distX > distY); 
} 

bool isInSameGroup(Vec4i x, vector<Vec4i> groupLines) {
    double m,c;
    const double trsh = 20.0;
    for(int i=0; i<groupLines.size(); i++){
        Vec4i L = groupLines[i];
        m = (double)(L[3]-L[1])/(L[2]-L[0]);
        c = (double) L[1]-m*L[0];
        if( !( x[1] > m*x[0] + c - trsh && x[1] < m*x[0] + c + trsh 
            && x[3] > m*x[2] + c - trsh && x[3] < m*x[2] + c + trsh ) ) return false;
    }
    return true;
}

void findFieldLines() {
    //find slope of each lines
    int i;
    for(i = 0; i < lines.size(); i++ ) {
        Vec4i l = lines[i];
        m[i] = (double)(l[3]-l[1])/(l[2]-l[0]);
        c[i] = l[1]-m[i]*l[0];
    }
        
    sort(lines.begin(), lines.end(), myfunction);
    
    //find real field line    
    for(i = 0; i < lines.size(); i++ ) {
        Vec4i L = lines[i];
        bool dupLine=false;
        for( int j = 0; j < groupLines.size(); j++) {
            if (isInSameGroup(L, groupLines[j])){
                groupLines[j].push_back(L);
                dupLine = true;
                break;
            }
        }
        if(!dupLine){
            vector<Vec4i> v;
            v.push_back(L);
            groupLines.push_back(v);
        }
    }
        
    // bestfit
    for(int g = 0; g < groupLines.size(); g++ ) {
        vector<Point2i> points;
        double d2, maxd2=-1.0;
        for( int j = 0; j < groupLines[g].size(); j++) {
            Vec4i L = groupLines[g][j];
            points.push_back(Point(L[0],L[1]));
            points.push_back(Point(L[2],L[3]));
        }
        // find the farest point
        Point2i start = Point(points[0].x, points[0].y);
        Point2i end = Point(points[1].x, points[1].y);
        for( int i = 0; i < points.size(); i++) 
            for( int j = i+1; j < points.size(); j++) {
                d2 = pow(points[i].x-points[j].x,2.0) + pow(points[i].y-points[j].y,2.0);
                if( d2 - maxd2 > 0.001 ) {
                    start = points[i];
                    end = points[j];
                    maxd2 = d2;
                }
            }
            
        // find best fit line
        Vec4f bestFitLine;
        fitLine(points, bestFitLine, CV_DIST_L12, 0, 0.01, 0.01);
        Vec4f L = bestFitLine;
        if( L[1]<0.001 && L[1]>-0.001 || L[0]<0.001 && L[0]>-0.001) {
            groupLines.erase(groupLines.begin()+g);
            g--;
            continue;
        }
        //cout<<i<<": "<<L[0]<<" "<<L[1]<<" "<<L[2]<<" "<<L[3]<<endl;
        line(myMat,Point(L[2],L[3]),Point(L[2]+L[0]*1000,L[3]+L[1]*1000),Scalar(0,0,255),1,CV_AA);
        line(myMat,Point(L[2],L[3]),Point(L[2]-L[0]*1000,L[3]-L[1]*1000),Scalar(0,0,255),1,CV_AA);
        
        M[g] = -L[1]/L[0]; // m = -dy/dx;
        C[g] = L[3] + M[g]*L[2]; // c = y + mx
        int x1 = start.x;
        int x2 = end.x;
        int y1 = (int)(- M[g]*x1 + C[g]);
        int y2 = (int)(- M[g]*x2 + C[g]);
        fieldLines.push_back(Vec4i(x1,y1,x2,y2));
        //circle(myMat, Point(x1,y1), 3, Scalar(255,0,0), -1);
        //circle(myMat, Point(x2,y2), 3, Scalar(255,0,0), -1);
    }
}
/*
void mapTo2d(){
    // homography
        vector<Point2f> pts_src, pts_dst;
        const uchar RIGHT = 2, CENTER = 1, LEFT = 0;
        // TRY TRY TRY
        //Mat dst_bw_image (height, width, CV_8UC1);
        
        uchar field_side = LEFT;
        
        for(int i=0; i<fieldLines.size(); i++){
            Vec4i L = fieldLines[i];
            //cout<<i<<": "<<L[0]<<" "<<L[1]<<" "<<L[2]<<" "<<L[3]
            //<<" "<<sqrt( pow((double)L[0]-L[2],2) + pow((double)L[1]-L[3],2))<<endl;
            circle(myMat, Point(L[0],L[1]), 3, Scalar(255,0,0), -1);
            circle(myMat, Point(L[2],L[3]), 3, Scalar(255,0,0), -1);
            //line(myMat,Point(L[2],L[3]),Point(L[2]+L[0]*100,L[3]+L[1]*100),Scalar(0,0,255),1,CV_AA);
        }
        if ( field_side==LEFT ) {
            pts_src.push_back(Point(570,180));
            pts_src.push_back(Point(width-1,200));
            pts_src.push_back(Point(0,240));
            pts_src.push_back(Point(616,230));
            
            int d;
            pts_dst.push_back(Point(0,0));
            pts_dst.push_back(Point(width-1,0));
            pts_dst.push_back(Point(0,height-1));
            pts_dst.push_back(Point(200,200));
            
            
        }
        
        
        
        // pts_src and pts_dst are vectors of points in source 
        // and destination images. They are of type vector<Point2f>. 
        // We need at least 4 corresponding points. 
        
        Mat h = findHomography(pts_src, pts_dst);
 
        // The calculated homography can be used to warp 
        // the source image to destination. im_src and im_dst are
        // of type Mat. Size is the size (width,height) of im_dst. 
        
        warpPerspective(bw_image, dst_bw_image, h, matSize);     
}
    */
    
void findRealFieldLines() {
    double slopeTrsh = 0.5;
    int x,y;
    
    // CASES
    if (fieldLines.size() == 2) {
        x = (C[1] - C[0]) / (M[1] - M[0]);
        y = -M[0]*refPt.x + C[0];
        if(M[0]<0  ) {
            vertLines.push_back(fieldLines[0]);
            horizLines.push_back(fieldLines[1]);
        }
        else {
            vertLines.push_back(fieldLines[1]);
            horizLines.push_back(fieldLines[0]);
        }
        circle(myMat, Point(x,y), 3, Scalar(255,0,0), -1);
        refPt.x = x; refPt.y = y;
    }
}

int main(int argc, char** argv)
{
    IplImage* frame=0;
    IplImage* frame_bw=0;     
    IplImage* frame_gray=0;
    CvCapture* capture=0;
    char key;  
    const int reduce = 1;
    int i,j,k,width, height, widthStep, channels;
   
    char const* avifile = "D:\\Google Drive\\SeniorProject\\6.mp4";
    capture = cvCaptureFromAVI(avifile);
    if(!capture)throw "Error when reading steam_avi";
    cvNamedWindow( "small", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
          
    
    for(int framePointer=0 ; framePointer<100 ; ++framePointer)
    {
        frame = cvQueryFrame( capture );
        if(!frame)break;
        CvSize size = cvSize(frame->width / reduce, frame->height / reduce);
        IplImage* tmpsize = cvCreateImage(size, frame->depth, frame->nChannels);
        cvResize(frame, tmpsize, CV_INTER_LINEAR);


        height = tmpsize->height;
        width = tmpsize->width;
        widthStep = tmpsize->widthStep;
        channels = tmpsize->nChannels;
        myMat = cvarrToMat(tmpsize);
		matSize = myMat.size();
        	
		
		for(i=0; i<height ;i++)
		for(j=0; j<width ;j++)
		for(k=0; k<3 ;k++)
		origin[i][j][k]=myMat.at<cv::Vec3b>(i,j)[k];
		
        
        for(i=2; i<height-2 ;i++)
        for(j=2; j<width-2 ;j++)
        for(k=0; k<3 ;k++) 
        //gaussian[i][j][k] = origin[i][j][k];        
        gaussian[i][j][k] = ( 41*origin[i][j][k]
                            + 26*(origin[i-1][j][k]+origin[i+1][j][k]+origin[i][j-1][k]+origin[i][j+1][k] )
                            + 16*(origin[i-1][j-1][k]+origin[i+1][j-1][k]+origin[i-1][j+1][k]+origin[i+1][j+1][k])
                            + 7*(origin[i-2][j][k]+origin[i+2][j][k]+origin[i][j-2][k]+origin[i][j+2][k])
                            + 1*(origin[i-2][j-2][k]+origin[i+2][j-2][k]+origin[i-2][j+2][k]+origin[i+2][j+2][k])
                            + 4*(origin[i-2][j-1][k]+origin[i-2][j+1][k]+origin[i+2][j-1][k]+origin[i+2][j+1][k]+origin[i-1][j-2][k]+origin[i+1][j-2][k]+origin[i-1][j+2][k]+origin[i+1][j+2][k])
                            ) /273;       
        

        for(i=0; i<height ;i++) 
        for(j=0; j<width ;j++) 
        if(gaussian[i][j][1]>gaussian[i][j][2]+10 && gaussian[i][j][1]>gaussian[i][j][0]+10)
        {grass[i][j][0]=gaussian[i][j][0];grass[i][j][1]=gaussian[i][j][1];grass[i][j][2]=gaussian[i][j][2];}
        else {grass[i][j][0]=0;grass[i][j][1]=0;grass[i][j][2]=0;}
         
			
        for(i=1; i<height ;i++)
        for(j=1; j<width ;j++)
        for(k=0; k<3 ;k++)
        gradient[i][j][k] = sqrt ( (double) ( 
                               pow ( 
                               (double) (grass[i-1][j+1][k]+grass[i+1][j+1][k]+grass[i][j+1][k]*2-grass[i-1][j-1][k]-grass[i+1][j-1][k]-grass[i][j-1][k]*2) 
                               , 2.0)
                               +  pow ( 
                               (double) (grass[i+1][j+1][k]+grass[i+1][j-1][k]+grass[i+1][j][k]*2-grass[i-1][j-1][k]-grass[i-1][j+1][k]-grass[i-1][j][k]*2)
                               , 2.0) 
                               ));
        	
		
        const int tshGD =30;
        for(i=0; i<height ;i++)
        for(j=0; j<width ;j++)
        if(gradient[i][j][0]>tshGD && gradient[i][j][1]>tshGD && gradient[i][j][2]>tshGD) 
        {tshGradient[i][j][0]=255;tshGradient[i][j][1]=255;tshGradient[i][j][2]=255;}
        //if(gradient[i][j][0]>tshGD) tshGradient[i][j][0]=255;
        else {tshGradient[i][j][0]=0;tshGradient[i][j][1]=0;tshGradient[i][j][2]=0;}
        

        Mat bw_image (height, width, CV_8UC1);
		for(i=0; i<height ;i++)
		for(j=0; j<width ;j++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][0];    
       
        HoughLinesP(bw_image, lines, 5, CV_PI/180, 200, 200, 10 );
        
        findFieldLines();
        
        findPlayer();
        
        findRealFieldLines();
 

        imshow("small", myMat);
        imshow("bw", bw_image);
        //imshow("2d", dst_bw_image);
        
        lines.clear();
        groupLines.clear();
        fieldLines.clear();
        horizLines.clear();
        vertLines.clear();
        
        key = cvWaitKey(33);
        if(key==27)break;
    }

   
    cvWaitKey(0); 
    cvDestroyWindow("bw");
    cvDestroyWindow( "small" );
    //cvDestroyWindow( "1C" );
    //cvDestroyWindow( "2d" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}


