#include <stdio.h>
#include <cv.h>
#include "highgui.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.14159265
#define Root2 1.4142135623
#define Root3 1.7320508075
using namespace cv;

Size matSize;
int origin[1000][2000][3];
int grass[1000][2000][3];
int gaussian[1000][2000][3];
double gradient[1000][2000][3];
int tshGradient[1000][2000][3];
int mark[1000][2000];
int qI[2000000],qJ[2000000];

  
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
    
    
int main(int argc, char** argv)
{
    IplImage* frame=0;
    IplImage* frame_bw=0;     
    IplImage* frame_gray=0;
    CvCapture* capture=0;
    char key;  
    const int reduce = 2;
    int i,j,k,width, height, widthStep, channels;
   
    char const* avifile = "D:\\Google Drive\\SeniorProject\\2.mp4";
    capture = cvCaptureFromAVI(avifile);
    if(!capture)throw "Error when reading steam_avi";
    cvNamedWindow( "small", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
          
    
    for(int framePointer=0 ; framePointer<1000 ; ++framePointer)
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
        Mat myMat (tmpsize);
		matSize = myMat.size();
        	
		
		for(i=0; i<matSize.height ;i++)
		for(j=0; j<matSize.width ;j++)
		for(k=0; k<3 ;k++)
		origin[i][j][k]=myMat.at<cv::Vec3b>(i,j)[k];
		
        
        for(i=2; i<matSize.height-2 ;i++)
        for(j=2; j<matSize.width-2 ;j++)
        for(k=0; k<3 ;k++) 
        //gaussian[i][j][k] = origin[i][j][k];        
        gaussian[i][j][k] = ( 41*origin[i][j][k]
                            + 26*(origin[i-1][j][k]+origin[i+1][j][k]+origin[i][j-1][k]+origin[i][j+1][k] )
                            + 16*(origin[i-1][j-1][k]+origin[i+1][j-1][k]+origin[i-1][j+1][k]+origin[i+1][j+1][k])
                            + 7*(origin[i-2][j][k]+origin[i+2][j][k]+origin[i][j-2][k]+origin[i][j+2][k])
                            + 1*(origin[i-2][j-2][k]+origin[i+2][j-2][k]+origin[i-2][j+2][k]+origin[i+2][j+2][k])
                            + 4*(origin[i-2][j-1][k]+origin[i-2][j+1][k]+origin[i+2][j-1][k]+origin[i+2][j+1][k]+origin[i-1][j-2][k]+origin[i+1][j-2][k]+origin[i-1][j+2][k]+origin[i+1][j+2][k])
                            ) /273;       
        

        for(i=0; i<matSize.height ;i++) 
        for(j=0; j<matSize.width ;j++) 
        if(gaussian[i][j][1]>gaussian[i][j][2]+10 && gaussian[i][j][1]>gaussian[i][j][0]+10)
        {grass[i][j][0]=gaussian[i][j][0];grass[i][j][1]=gaussian[i][j][1];grass[i][j][2]=gaussian[i][j][2];}
        else {grass[i][j][0]=0;grass[i][j][1]=0;grass[i][j][2]=0;}

			
        for(i=1; i<matSize.height-1 ;i++)
        for(j=1; j<matSize.width-1 ;j++)
        for(k=0; k<3 ;k++)
        gradient[i][j][k] = sqrt ( (double) ( 
                               pow ( 
                               (double) (grass[i-1][j+1][k]+grass[i+1][j+1][k]+grass[i][j+1][k]*2-grass[i-1][j-1][k]-grass[i+1][j-1][k]-grass[i][j-1][k]*2) 
                               , 2.0)
                               +  pow ( 
                               (double) (grass[i+1][j+1][k]+grass[i+1][j-1][k]+grass[i+1][j][k]*2-grass[i-1][j-1][k]-grass[i-1][j+1][k]-grass[i-1][j][k]*2)
                               , 2.0) 
                               ));
                               
			
        const int tshGD = 30;
        for(i=0; i<matSize.height ;i++)
        for(j=0; j<matSize.width ;j++)
        if(gradient[i][j][0]>tshGD && gradient[i][j][1]>tshGD && gradient[i][j][2]>tshGD) 
        {tshGradient[i][j][0]=255;tshGradient[i][j][1]=255;tshGradient[i][j][2]=255;}
        else {tshGradient[i][j][0]=0;tshGradient[i][j][1]=0;tshGradient[i][j][2]=0;}


        Mat bw_image (height, width, CV_8UC1);
		matSize = bw_image.size();
		for(i=0; i<matSize.height ;i++)
		for(j=0; j<matSize.width ;j++)
		for(k=0; k<3 ;k++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][k];    
       
        vector<Vec4i> lines;
        HoughLinesP(bw_image, lines, 5, CV_PI/180, 50, 100, 5 );
        
///////////////////////////////////////////////////////////////////////////////        
        //find slope of each lines
        double m[1000];
        double c[1000];
        for( i = 0; i < lines.size(); i++ ) {
             Vec4i l = lines[i];
             m[i] = (double)(l[1]-l[3])/(l[2]-l[0]);
             c[i] = l[1]-m[i]*l[0];
        }
        
        //find real field line
        double M[1000];
        double C[1000];
        int fieldLineSize=0;
        
        for( i = 0; i < lines.size(); i++ ) 
        {
             Vec4i l1 = lines[i];
             int repeatC=0,repeatM=0;
             for( int j = 0; j < fieldLineSize; j++) 
             {
                  if(m[i]-M[j]<0.3 && m[i]-M[j]>-0.3)
                  {
                     if(c[i]-C[j]<10 && c[i]-C[j]>-10)repeatC=1;
                     else repeatM=1;
                     break;        
                  }
             }
                  if(repeatC==1){}
             else if(repeatM==1){}
             else {M[fieldLineSize]=m[i];C[fieldLineSize++]=c[i];}
        }
///////////////////////////////////////////////////////////////////////////////////////////

        //draw lines
        for(i=0; i<lines.size() ;i++) 
        {
            Vec4i l = lines[i];
            line( myMat, Point(l[0], l[1]), Point(l[2], l[3]), Scalar(0,0,255), 1, CV_AA);
            //std::cout << m[i] << " ";
        }
        //std::cout << "\n";
        
        findPlayer();
 
        for(i=0; i<matSize.height ;i++)
        for(j=0; j<matSize.width ;j++)
        for(k=0; k<3 ;k++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][k];
        
        
        imshow("small", myMat);
        imshow("bw", bw_image);
        
        key = cvWaitKey(33);
        if(key==27)break;
    }

   
    cvWaitKey(0); 
    cvDestroyWindow("bw");
    cvDestroyWindow( "small" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}
