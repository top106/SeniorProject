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
int width, height, widthStep, channels;
int origin[1000][2000][3];
int grass[1000][2000][3];
int gaussian[1000][2000][3];
double gradient[1000][2000][3];
uchar tshGradient[1000][2000][3];
int qI[2000000],qJ[2000000],pI[100][10000],pJ[100][10000];
int mark[1000][2000],mark2[1000][2000];
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
Vec4i refVertLine;
Vec4i prevVertLine;

int abs(int x,int y){if(x>y)return x-y;else return y-x;} 
  
void findPlayer()
{
     int A,I,J,i,j,k,r,p;
     int R[1000],G[1000],B[1000],X[1000],Y[1000],team[1000],area[100],Harea[100];
     int H=matSize.height,W=matSize.width;
     int numArea=0,MinH=H,max;
     
     for(i=0;i<H;i++)
     for(j=0;j<W;j++)
     {mark[i][j]=0;mark2[i][j]=0;}
     
     for(i=0;i<H;i++)
     for(j=0;j<W;j++)
     {
        if(tshGradient[i][j][0]==0&&mark[i][j]==0)
        {
           r=0;p=0;max=0;
           mark[i][j]=0; 
           qI[r]=i;qJ[r]=j;r++;            
           
           while(p<r)
           {
               I=qI[p];J=qJ[p];p++;
               if(I-i>max)max=I-i;
               if(J+1<W){if(tshGradient[I][J+1][0]==0&&mark[I][J+1]==0){qI[r]=I;qJ[r]=J+1;mark[qI[r]][qJ[r]]=1;r++;}}
               if(J-1>=0){if(tshGradient[I][J-1][0]==0&&mark[I][J-1]==0){qI[r]=I;qJ[r]=J-1;mark[qI[r]][qJ[r]]=1;r++;}}
               if(I+1<H){if(tshGradient[I+1][J][0]==0&&mark[I+1][J]==0){qI[r]=I+1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
               if(I-1>=0){if(tshGradient[I-1][J][0]==0&&mark[I-1][J]==0){qI[r]=I-1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
           }
           
           if(10<r && r<2000)
           {
               Harea[numArea]=max;
               area[numArea]=r;
               X[numArea]=j;
               Y[numArea]=i;
               
               for(k=0;k<r;k++)
               {pI[numArea][k]=qI[k];pJ[numArea][k]=qJ[k];}
               numArea++;     
               
               r=0;p=0;
               qI[r]=i-1;qJ[r]=j;r++;
               mark2[i-1][j]=0;
               while(p<r)
               {
                  I=qI[p];J=qJ[p];p++;
                  tshGradient[I][J][0]=1;tshGradient[I][J][1]=1;tshGradient[I][J][2]=1;
                  if(J+1<W){if(tshGradient[I][J+1][0]==255&&mark2[I][J+1]==0){qI[r]=I;qJ[r]=J+1;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(J-1>=0){if(tshGradient[I][J-1][0]==255&&mark2[I][J-1]==0){qI[r]=I;qJ[r]=J-1;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(I+1<H){if(tshGradient[I+1][J][0]==255&&mark2[I+1][J]==0){qI[r]=I+1;qJ[r]=J;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(I-1>=0){if(tshGradient[I-1][J][0]==255&&mark2[I-1][J]==0){qI[r]=I-1;qJ[r]=J;mark2[qI[r]][qJ[r]]=1;r++;}}
               }
               //for(k=0;k<r;k++){tshGradient[qI[k]][qJ[k]][0]=120;tshGradient[qI[k]][qJ[k]][1]=120;tshGradient[qI[k]][qJ[k]][2]=120;}             
           }
        }   
     }
     
     
     /*
     max=0;
     for(i=0;i<numArea;i++)
     {
         int count=0;
         for(j=0;j<numArea;j++)
         if(abs(Harea[i]-Harea[j])<30)count++;
         if(count>max){max=count;I=i;}
     }
     
     int min=Harea[I];
     for(i=0;i<numArea;i++)
     if(abs(Harea[i]-Harea[I])<30&&Harea[i]<min)min=Harea[i];
     */
     
     //for(i=0;i<numArea;i++)
     //if(abs(Harea[i]-Harea[I])<30)
     //for(k=0;k<area[i];k++){tshGradient[pI[i][k]][pJ[i][k]][0]=120;tshGradient[pI[i][k]][pJ[i][k]][1]=120;tshGradient[pI[i][k]][pJ[i][k]][2]=120;}
     
/*     
     for(i=0;i<numArea;i++)
     //if(abs(Harea[i]-Harea[I])<30)
     {
         R[i]=0;G[i]=0;B[i]=0;A=0;
         
         for(k=0;k<area[i];k++)
         //if(pI[i][k]-Y[i]<min)
         {B[i]+=origin[pI[i][k]][pJ[i][k]][0];G[i]+=origin[pI[i][k]][pJ[i][k]][1];R[i]+=origin[pI[i][k]][pJ[i][k]][2];A++;}
                     
         R[i]/=A;G[i]/=A;B[i]/=A;
         printf("%d: (%d,%d) [%d,%d,%d]   ",framePointer,X[i],Y[i],R[i],G[i],B[i]);
     }

     
     int tsh=40;
     int numManySameTeam=0;
     int sameTeam[1000];
     int manySameTeam[1000],MST;
     for(i=0;i<numArea;i++)
     {
        int count=0;
        for(j=0;j<numArea;j++)
        if(abs(R[i]-R[j])+abs(G[i]-G[j])+abs(B[i]-B[j])<tsh)count++;
        sameTeam[i]=count;
        //printf("\\%d %d %d %d\\",count,R[i],G[i],B[i]);
        
        if(count>1)manySameTeam[numManySameTeam++]=i;
     }
     
     max=0;
     for(i=0;i<numManySameTeam;i++)
     if(sameTeam[manySameTeam[i]]>max){max=sameTeam[manySameTeam[i]];I=manySameTeam[i];} 
     
     max=0;
     for(i=0;i<numManySameTeam;i++)
     {
        MST = manySameTeam[i];
        if(sameTeam[MST]>max && MST!=I && abs(R[MST]-R[I])+abs(G[MST]-G[I])+abs(B[MST]-B[I])>=tsh)
        {max=sameTeam[MST];J=MST;}
     }
     //printf("--%d--%d,%d,%d--%d--%d,%d,%d\n",sameTeam[I],R[I],G[I],B[I],sameTeam[J],R[J],G[J],B[J]);
     
     for(i=0;i<numManySameTeam;i++)
     {
        MST = manySameTeam[i];
        if(abs(R[MST]-R[I])+abs(G[MST]-G[I])+abs(B[MST]-B[I])<tsh)team[MST]=I;
        else if(abs(R[MST]-R[J])+abs(G[MST]-G[J])+abs(B[MST]-B[J])<tsh)team[MST]=J;                      
     }
     
     for(i=0;i<numManySameTeam;i++)
     {
        MST = manySameTeam[i];
        if(team[MST]==I)
           for(j=X[MST];j<=X[MST]+50;j++)
           {tshGradient[Y[MST]][j][0]=120;tshGradient[Y[MST]][j][1]=120;tshGradient[Y[MST]][j][2]=120;
           tshGradient[Y[MST]-1][j][0]=120;tshGradient[Y[MST]-1][j][1]=120;tshGradient[Y[MST]-1][j][2]=120;
           tshGradient[Y[MST]+1][j][0]=120;tshGradient[Y[MST]+1][j][1]=120;tshGradient[Y[MST]+1][j][2]=120;}
        else if(team[MST]==J)
           for(j=Y[MST];j<=Y[MST]+50;j++)
           {tshGradient[j][X[MST]][0]=120;tshGradient[j][X[MST]][1]=120;tshGradient[j][X[MST]][2]=120;
           tshGradient[j][X[MST]-1][0]=120;tshGradient[j][X[MST]-1][1]=120;tshGradient[j][X[MST]-1][2]=120;
           tshGradient[j][X[MST]+1][0]=120;tshGradient[j][X[MST]+1][1]=120;tshGradient[j][X[MST]+1][2]=120;}
     }
     //printf(".");
*/
}



void findBall()
{
     int iBall,jBall;
     int i,j,k,l,r=7,score;
     int H=matSize.height;
     int W=matSize.width;
     double BB,GG,RR,DD=0;
     
     for(i=r;i<H-r;i++)
     for(j=r;j<W-r;j++)
     {
        score=0;
        /*BB=0;GG=0;RR=0;DD=0;
        for(k=i-r;k<=i+r;k++)
        for(l=j-r;l<=j+r;l++)
        if((k-i)*(k-i)+(l-j)*(l-j)<=r*r)
        {BB+=origin[k][l][0];GG+=origin[k][l][1];RR+=origin[k][l][2];DD++;}*/
        
        for(k=i-r;k<=i+r;k++)
        for(l=j-r;l<=j+r;l++)
        if((k-i)*(k-i)+(l-j)*(l-j)<=r*r){if(tshGradient[k][l][0]==255)score++;}
        else{if(tshGradient[k][l][0]==0||tshGradient[k][l][0]==1)score++;}
           
        if(score*100>r*r*4*90)
        {
           iBall = i;
           jBall = j;
           //printf("%d %d",i,j);
           for(k=i-r;k<=i+r;k++)
           for(l=j-r;l<=j+r;l++)
           if((k-i)*(k-i)+(l-j)*(l-j)<=r*r)
           {tshGradient[k][l][0]=120;tshGradient[k][l][1]=120;tshGradient[k][l][2]=120;}
        }
     } 
} 
double round(double d)
{
  return floor(d + 0.5);
}

bool myfunction (Vec4i x,Vec4i y) { 
     double distX = sqrt( pow((double)x[0]-x[2],2) + pow((double)x[1]-x[3],2));
     double distY = sqrt( pow((double)y[0]-y[2],2) + pow((double)y[1]-y[3],2));
     return (distX > distY); 
} 


void findVerticalLine() {
    sort(lines.begin(), lines.end(), myfunction);
    int x1,y1,x2,y2;
    double m;
    Vec4i l;
    // find longest vertical line 
    for(int i = 0; i < lines.size(); i++) {
        l = lines[i];
        x1=l[0]; y1=l[1]; x2=l[2]; y2=l[3];
        m = (double)(y2-y1)/(x2-x1);
        if(x2==x1) continue;
        if(m > 1 || m < -1) {
            refVertLine = l;
            prevVertLine = l;
            break;
        }
    }
    x1=refVertLine[0];
    y1=refVertLine[1];
    x2=refVertLine[2];
    y2=refVertLine[3];
    if(x1==0 && y1==0 && x2==0 && y2==0) refVertLine = prevVertLine;
    // draw
    x1 = refVertLine[0];
    y1 = refVertLine[1];
    x2 = refVertLine[2];
    y2 = refVertLine[3];
    if(x1==0 && y1==0 && x2==0 && y2==0) return; // prevVertLine is also 0 0 0 0
    line(myMat,Point(x1,y1),Point(x2,y2),Scalar(0,0,255),1,CV_AA);
    circle(myMat, Point(x1,y1), 3, Scalar(255,0,0), -1);
    circle(myMat, Point(x2,y2), 3, Scalar(255,0,0), -1);
    
    //print
    m = (y2-y1)/(x2-x1); 
    cout<<m<<" "<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<"\n\n";
    
}


void clearVec4i (Vec4i &v){
     v[0]=0; v[1]=0; v[2]=0; v[3]=0;
}


int main(int argc, char** argv)
{
    IplImage* frame=0;
    IplImage* frame_bw=0;     
    IplImage* frame_gray=0;
    CvCapture* capture=0;
    char key;
    const int reduce = 1;
    int i,j,k;
    
    char const* avifile = "D:\\Google Drive\\SeniorProject\\Roma vs Fiorentina - cut2.avi";
    capture = cvCaptureFromAVI(avifile);
    if(!capture)throw "Error when reading steam_avi";
    cvNamedWindow( "small", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
          
    
    for(int framePointer=0 ; framePointer<200 ; ++framePointer)
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
        else {tshGradient[i][j][0]=0;tshGradient[i][j][1]=0;tshGradient[i][j][2]=0;}
        

        Mat bw_image (height, width, CV_8UC1);
		for(i=0; i<height ;i++)
		for(j=0; j<width ;j++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][0];
       
        HoughLinesP(bw_image, lines, 5, CV_PI/180, 200, 200, 10 );
        
        findVerticalLine();
        
        findPlayer();
        
        findBall(); 
        
        for(i=0; i<matSize.height ;i++)
        for(j=0; j<matSize.width ;j++)
        for(k=0; k<3 ;k++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][k];
        imshow("small", myMat);
        imshow("bw", bw_image);
        
        lines.clear();
        clearVec4i(refVertLine);
        
        key = cvWaitKey(33);
        if(key==27)break;
    }

   
    cvWaitKey(0); 
    cvDestroyWindow("bw");
    cvDestroyWindow( "small" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}


