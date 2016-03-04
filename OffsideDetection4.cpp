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
vector<Vec4i> vertLines;
vector<Vec4i> horizLines;
Vec4i uppermostHoriz;
Vec4i aHorizLine;
Vec4i leftmostLine;
Vec4i rightmostLine;
Vec4i aVertLine;
//Point refPt;

double round(double d)
{
  return floor(d + 0.5);
}

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
        double x = round(L[2]), y = round(L[3]);
        //line(myMat,Point(x,y),Point(x+L[0]*1000,y+L[1]*1000),Scalar(0,0,255),1,CV_AA);
        //line(myMat,Point(x,y),Point(x-L[0]*1000,y-L[1]*1000),Scalar(0,0,255),1,CV_AA);
        
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


void groupParallelLines() { // and delete curve line of the penalty box
    const double trsh = 0.1;
    double  m;
    vector<double> mAvg;
    vector<int> groupSize;
    vector< vector<int> > index;
    bool hasGroup;
    for(int i=0; i<fieldLines.size(); i++) {
        hasGroup=false;
        if(M[i]>1) m=1.0+1.0/M[i];
        else m=M[i];
        for(int j=0; j<groupSize.size(); j++) {
            if( m<mAvg[j]+trsh && m>mAvg[j]-trsh) {
                hasGroup=true;
                mAvg[j]=(mAvg[j]*groupSize[j]+m)/(groupSize[j]+1);
                groupSize[j]++;
                index[j].push_back(i);
                break;
            }
        }
        if(!hasGroup){
            mAvg.push_back(m);
            groupSize.push_back(1);
            vector<int> v;
            v.push_back(i);
            index.push_back(v);
        }
    }
    
    if(mAvg.size() > 3) {
        cout<<"group of lines > 3, YOU MAD!?\n\n";
        return;
    }
    if(mAvg.size() == 1) {
        cout<<"detect only 1 group of lines, YOU MAD!?\n\n";
        return;
    }
    if(mAvg.size() == 0) {
        cout<<"cannot detect any line, YOU MAD!?\n\n";
        return;
    }
    // max+ max- OR min- min+
    //double max_pos=0.0, min_pos=2000000.0, max_neg=-2000000.0, min_neg=0.0, d;
    
    if(mAvg.size()==3) {
        vector<int> m_pos, m_neg; // index of group
        for(int i=0; i<mAvg.size(); i++){
            if(mAvg[i]>0) m_pos.push_back(i);
            else m_neg.push_back(i);
        }
        int i;
        if(m_pos.size() == 2) {
            if(mAvg[m_pos[0]] > mAvg[m_pos[1]]) 
                for(i=0; i<index[m_pos[0]].size(); i++) 
                    vertLines.push_back(fieldLines[index[m_pos[0]][i]]);
            else 
                for(i=0; i<index[m_pos[1]].size(); i++)
                    horizLines.push_back(fieldLines[index[m_pos[1]][i]]);
            for(i=0; i<index[m_neg[0]].size(); i++) 
                horizLines.push_back(fieldLines[index[m_neg[0]][i]]);
        }
        else if(m_neg.size() == 2){
            if(mAvg[m_neg[0]] < mAvg[m_neg[1]]) 
                for(i=0; i<index[m_neg[0]].size(); i++) 
                    vertLines.push_back(fieldLines[index[m_neg[0]][i]]);
            else 
                for(i=0; i<index[m_neg[1]].size(); i++)
                    horizLines.push_back(fieldLines[index[m_neg[1]][i]]);
            for(i=0; i<index[m_pos[0]].size(); i++) 
                horizLines.push_back(fieldLines[index[m_pos[0]][i]]);
        }
        else {
             cout<<"YOU MAD!?\n\n";
             return;
        }
    }
    
    if(mAvg.size()==2) { // assume the curve line is harder to detect than normal lines
        int i;
        if( mAvg[0]>0 && fabs(mAvg[0])>fabs(mAvg[1]) || mAvg[0]<0 && fabs(mAvg[0])>fabs(mAvg[1])) {
            for(i=0; i<index[0].size(); i++) 
                vertLines.push_back(fieldLines[index[0][i]]);
            for(i=0; i<index[1].size(); i++) 
                horizLines.push_back(fieldLines[index[1][i]]);
        }
        else if( mAvg[1]>0 && fabs(mAvg[1])>fabs(mAvg[0]) || mAvg[1]<0 && fabs(mAvg[1])>fabs(mAvg[0])) {
            for(i=0; i<index[1].size(); i++) 
                vertLines.push_back(fieldLines[index[1][i]]);
            for(i=0; i<index[0].size(); i++) 
                horizLines.push_back(fieldLines[index[0][i]]);
        }
        else {
             cout<<"both are vertical or horizontal lines, YOU MAD!?\n\n";
             return;
        }
    }
    // print
    int x1,y1,x2,y2; 
    Vec4i L;
    cout<<"vertical:\n";
    for(int i=0; i<vertLines.size(); i++) {
        L = vertLines[i];
        x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
        cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        line(myMat,Point(x1,y1),Point(x2,y2),Scalar(0,0,255),1,CV_AA);
        circle(myMat, Point(x1,y1), 3, Scalar(255,0,0), -1);
        circle(myMat, Point(x2,y2), 3, Scalar(255,0,0), -1);
    }
    cout<<"horizontal:\n";
    for(int i=0; i<horizLines.size(); i++) {
        L = horizLines[i];
        x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
        cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        line(myMat,Point(x1,y1),Point(x2,y2),Scalar(0,0,255),1,CV_AA);
        circle(myMat, Point(x1,y1), 3, Scalar(255,0,0), -1);
        circle(myMat, Point(x2,y2), 3, Scalar(255,0,0), -1);
    }
    
    // find grass point
    /*Point leftuppermost, rightuppermost;
    bool found=false;
    int i,j;
    for(i=0; i<width ;i++) {
        for(j=0; j<height ;j++) {
            if(grass[j][i][0]!=0 && grass[j][i][1]!=0 && grass[j][i][2]!=0) {
                leftuppermost = Point(i,j);
                found=true;
                break;
            }
        }
        if(found) break;
    }
    found=false;
    for(i=width-1; i>=0 ;i--) {
        for(j=0; j<height ;j++) {
            if(grass[j][i][0]!=0 && grass[j][i][1]!=0 && grass[j][i][2]!=0) {
                rightuppermost = Point(i,j);
                found=true;
                break;
            }
        }
        if(found) break;
    }*/
    //cout<<"left uppermost point: "<<leftuppermost.x<<" "<<leftuppermost.y<<endl
    //<<"right uppermost point: "<<rightuppermost.x<<" "<<rightuppermost.y<<endl;
    
    // y avg is smallest
    L = vertLines[0];
    x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
    double mVert = (y1-y2)/(x1-x2);
    L = horizLines[0];
    x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
    double mHoriz = (y1-y2)/(x1-x2);
    double uppermostY=20000.0,y;
    
    if(mVert>0 && mHoriz<0) { // left side of the field
        for(int i=0; i<vertLines.size(); i++) {
            L = vertLines[0];
            x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
            y=(y1+y2)/2;
            if(y<uppermostY){
                uppermostY=y;
                // assign leftuppermostLine with L
            }
        }
    }
    else if(mVert<0 && mHoriz>0) { // right side of the field
        
        
        
        //c=(double)y1-m*x1;
        //y1=m*x+c;
    }
    else { //you mad?
        cout<<"YOU MAD!?"<<endl;
        return;
    }
    for(int i=0; i<horizLines.size(); i++) {
        L = horizLines[i];
        x1=L[0]; y1=L[1]; x2=L[2]; y2=L[3];
        cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        line(myMat,Point(x1,y1),Point(x2,y2),Scalar(0,0,255),1,CV_AA);
        circle(myMat, Point(x1,y1), 3, Scalar(255,0,0), -1);
        circle(myMat, Point(x2,y2), 3, Scalar(255,0,0), -1);
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
    
// Finds the intersection of two lines, or returns false.
// The lines are defined by (o1, p1) and (o2, p2).
bool intersection(Point2f o1, Point2f p1, Point2f o2, Point2f p2, Point2f &r)
{
    Point2f x = o2 - o1;
    Point2f d1 = p1 - o1;
    Point2f d2 = p2 - o2;

    float cross = d1.x*d2.y - d1.y*d2.x;
    if (abs(cross) < /*EPS*/1e-8)
        return false;

    double t1 = (x.x * d2.y - x.y * d2.x)/cross;
    d1.x*=t1; d1.y*=t1;
    r = o1 + d1;
    return true;
}
    
void findRealFieldLines() {
    const double slopeTrsh = 0.5;
    double x,y;
    
    // CASES
    /*if (fieldLines.size() == 2) {
        x = round((C[1] - C[0]) / (M[1] - M[0]));
        y = round(-M[0]*refPt.x + C[0]);
        //if(M[i] > 1) M1 = 1 + 1/M[i];
        Vec4i L1 = fieldLines[0];
        Vec4i L2 = fieldLines[1];
        //cout<<"1: "<<L1[0]<<" "<<L1[1]<<" "<<L1[2]<<" "<<L1[3]<<endl;
        //cout<<"2: "<<L1[0]<<" "<<L1[1]<<" "<<L1[2]<<" "<<L1[3]<<endl;
        
        Point2f r;
        bool intersects = intersection(Point(L1[0],L1[1]),Point(L1[2],L1[3]),Point(L2[0],L2[1]),Point(L2[2],L2[3]),r);
        if (intersects) {
                       
                       
        
        if( M[i]<0.5 && M[i]>-0.5 ) {
            if( L[1] < height/2 && L[3] < height/2 && (L[0] > width*2/3 || L[2] > width*2/3) )
            upTouchLine = L;
                
        }
        else if( M[i] < -2 ) {
            if( L[1] < height/2 && L[3] < height/2 && L[0] > width*2/3 && L[2] < 20 )
            cout<<"";
                    //upTouchLine = L;
                  //vertLines.push_back(fieldLines[1]);
                  //horizLines.push_back(fieldLines[0]);
        }
        circle(myMat, Point(L[0],L[1]), 3, Scalar(255,0,0), -1);
        circle(myMat, Point(L[2],L[3]), 3, Scalar(255,0,0), -1);
        circle(myMat, Point(x,y), 3, Scalar(255,0,0), -1);
        refPt.x = x; refPt.y = y;
    }
    cout<<upTouchLine[0]<<" "<<upTouchLine[1]<<" "<<upTouchLine[2]<<" "<<upTouchLine[3]<<endl
        <<downTouchLine[0]<<" "<<downTouchLine[1]<<" "<<downTouchLine[2]<<" "<<downTouchLine[3]<<endl
        <<leftGoalLine[0]<<" "<<leftGoalLine[1]<<" "<<leftGoalLine[2]<<" "<<leftGoalLine[3]<<endl
        <<rightGoalLine[0]<<" "<<rightGoalLine[1]<<" "<<rightGoalLine[2]<<" "<<rightGoalLine[3]<<endl
        <<halfwayLine[0]<<" "<<halfwayLine[1]<<" "<<halfwayLine[2]<<" "<<halfwayLine[3]<<endl;
    cout << "------------------------------------------------\n";*/
}

void clearVec4i (Vec4i &v){
     v[0]=0; v[1]=0; v[2]=0; v[3]=0;
}

void swapPoint (Vec4i &v){
     int tmp;
     tmp = v[0]; v[0] = v[2]; v[2] = tmp;
     tmp = v[1]; v[1] = v[3]; v[3] = tmp;
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
        
        groupParallelLines();
 

        imshow("small", myMat);
        imshow("bw", bw_image);
        //imshow("2d", dst_bw_image);
        
        lines.clear();
        groupLines.clear();
        fieldLines.clear();
        vertLines.clear();
        horizLines.clear();
        clearVec4i(uppermostHoriz);
        clearVec4i(aHorizLine);
        clearVec4i(leftmostLine);
        clearVec4i(rightmostLine);
        clearVec4i(aVertLine);
        
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


