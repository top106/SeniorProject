 #include <stdio.h>
#include <cv.h>
#include "highgui.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <string>
#define PI 3.14159265
#define Root2 1.4142135623
#define Root3 1.7320508075
using namespace cv;
using namespace std;

// least line's slope that will be vertical
const double vertSlopeTrsh = 0.31;
// least pixels from margins that will be considered margins
const int marginTrsh = 15;

const int ballRTrsh = 3;
const int grassTrsh = 120;
const int tshGD =40;
const double fontScale = 0.5;
const int FRAMES_WANTED = 200;
// offence directions
const int NA = -1;
const int LEFT = 0;
const int RIGHT = 1;
// colors
const Scalar BLUE = Scalar(255,0,0);
const Scalar RED = Scalar(0,0,255);
const Scalar PINK = Scalar(153,153,255);
const Scalar PURPLE = Scalar(102,0,102);
const Scalar BLACK = Scalar(0,0,0);
const Scalar GRAY = Scalar(192,192,192);
const Scalar WHITE = Scalar(255,255,255);
const Scalar YELLOW = Scalar(0,255,255);
const Scalar LIGHTBLUE = Scalar(255,255,0);


Size matSize;
int width, height, widthStep, channels, framePointer;
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
vector<Vec4i> fieldLines;
vector< vector<Vec4i> > groupLines;
vector<Vec4i> vertLines;
vector<Vec4i> prevVertLines;
Point2f intPoint, prevIntPoint;
int offenseDirection = NA;
bool parallel = false;
int yLeftMargin = -1, yRightMargin = -1;
int sameTeam[1000],manySameTeam[1000],MST;
int numArea,numManySameTeam,R[10000],G[10000],B[10000],X[10000],Y[10000],team[10000],area[10000],Harea[10000],Warea[10000];     
int thisR,lastR,iBall,jBall,IBall=-1,JBall=-1,lastI;
int lastDFIndex;
vector<int> offsidePlyrIndexes; 

// returns a slope of a Vec4i line
double slope(Vec4i L) {
    return (double) (L[1]-L[3])/(L[2]-L[0]+0.01);
}
  
void findOffsidePlayer()
{
     int i,j,k,I,J,M,Count[10000];
     double min=-1.0,max1=1,max2=-1,min1=-1,min2=1;
     int max=-1,Max1=-1,Max2=-1,Min1=-1,Min2=-1;
     //uncomment for anartovic
     //if(lastR>=3) 
     {for(i=0;i<numArea;i++)
     if((abs(X[i]-JBall)+abs(Y[i]-IBall)<min||min==-1)&&team[i]!=-1){min=abs(X[i]-JBall)+abs(Y[i]-IBall);I=i;}}
     //else{I=lastI;}
     
     //lastI=I;
     
     for(i=0;i<numArea;i++)Count[i]=0;
     for(i=0;i<numArea;i++)if(team[i]!=team[I]&&team[i]!=-1)Count[team[i]]++;
     for(i=0;i<numArea;i++)if(Count[i]>max){max=Count[i];J=i;}

     for(i=0;i<numArea;i++)
     {
        if(team[i]==team[J])
        {
           double m = ((double)(intPoint.y)+(double)(Y[i]))/((double)(intPoint.x)-(double)(X[i]));
           char output[50];
           snprintf(output,50,"%lf",m);
           putText(myMat, output, Point(X[i]+10, Y[i]), FONT_HERSHEY_SIMPLEX, fontScale, GRAY, 2 );
           line(myMat,Point(intPoint.x,-intPoint.y),Point(X[i],Y[i]),PINK,1,CV_AA);
           if(offenseDirection==1)
           {
              if(m<0){if(m>max1||max1==1){max1=m;Max1=i;}}
              else{if(m>max2){max2=m;Max2=i;}}
           }
           else if(offenseDirection==0)
           {
              if(m>0){if(m<min1||min1==-1){min1=m;Min1=i;}}
              else{if(m<min2){min2=m;Min2=i;}}
           }         
        }                      
     }
     
          if(offenseDirection==1){if(Max1==-1)M=Max2;else M=Max1;}
     else if(offenseDirection==0){if(Min1==-1)M=Min2;else M=Min1;}
     
     lastDFIndex = M;
     
     for(j=0;j<area[M];j++)
     for(k=0;k<3;k++)
     tshGradient[pI[M][j]][pJ[M][j]][k]=120; 
     
     //printf("%d %d\n",X[J],Y[J]);
     
     for(i=0;i<numArea;i++)
     {
        if(team[i]==team[I])
        {
           double m = ((double)(intPoint.y)+(double)(Y[i]))/((double)(intPoint.x)-(double)(X[i]));
           
           if(offenseDirection==1)
           {
              if((Max1==-1&&(m<0||m>max2))||(Max1!=-1&&m<0&&m>max1)) 
              {
                  for(j=0;j<area[i];j++)
                  for(k=0;k<3;k++)
                  tshGradient[pI[i][j]][pJ[i][j]][k]=120;
                  offsidePlyrIndexes.push_back(i);
              }
           }
           else if(offenseDirection==0)
           {
              if((Min1==-1&&(m>0||m<min2))||(Min1!=-1&&m>0&&m<min1)) 
              {
                  for(j=0;j<area[i];j++)
                  for(k=0;k<3;k++)
                  tshGradient[pI[i][j]][pJ[i][j]][k]=120;
                  offsidePlyrIndexes.push_back(i); 
              }
           }
        }
     }
}


void findPlayer()
{
     int A,I,J,i,j,k,l,r,p,z;
     int H=matSize.height,W=matSize.width;
     int min,max,maxJ,minJ,count;
     
     numArea=0;numManySameTeam=0;
     
     for(i=0;i<H;i++)
     for(j=0;j<W;j++)
     {mark[i][j]=0;mark2[i][j]=0;}
     
     for(i=H-1;i>=0;i--)
     for(j=0;j<W;j++)
     {
        if(tshGradient[i][j][0]==0&&mark[i][j]==0)
        {
           r=0;p=0;max=0;
           minJ=W;maxJ=0;
           mark[i][j]=0; 
           qI[r]=i;qJ[r]=j;r++;            
           
           while(p<r)
           {
               I=qI[p];J=qJ[p];p++;
               if(i-I>max)max=i-I;
               if(J>maxJ)maxJ=J;
               if(J<minJ)minJ=J;
               
               if(J+1<W){if(tshGradient[I][J+1][0]==0&&mark[I][J+1]==0){qI[r]=I;qJ[r]=J+1;mark[qI[r]][qJ[r]]=1;r++;}}
               if(J-1>=0){if(tshGradient[I][J-1][0]==0&&mark[I][J-1]==0){qI[r]=I;qJ[r]=J-1;mark[qI[r]][qJ[r]]=1;r++;}}
               if(I+1<H){if(tshGradient[I+1][J][0]==0&&mark[I+1][J]==0){qI[r]=I+1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
               if(I-1>=0){if(tshGradient[I-1][J][0]==0&&mark[I-1][J]==0){qI[r]=I-1;qJ[r]=J;mark[qI[r]][qJ[r]]=1;r++;}}
           }
           
           if(20<r && r<5000)
           {
               if(36<r && r<1500 && max>6 && maxJ-minJ>6)
               {
                   Harea[numArea]=max;
                   Warea[numArea]=maxJ-minJ;
                   area[numArea]=r;
                   X[numArea]=j;
                   Y[numArea]=i;
                   
                   for(k=0;k<r;k++)
                   {pI[numArea][k]=qI[k];pJ[numArea][k]=qJ[k];}
                   
                   R[numArea]=0;G[numArea]=0;B[numArea]=0;A=0;
                   for(k=0;k<r;k++)
                   {
                      B[numArea]+=origin[pI[numArea][k]][pJ[numArea][k]][0];
                      G[numArea]+=origin[pI[numArea][k]][pJ[numArea][k]][1];
                      R[numArea]+=origin[pI[numArea][k]][pJ[numArea][k]][2];
                      A++;
                   }
                   R[numArea]/=A;G[numArea]/=A;B[numArea]/=A;
                   numArea++;
               }     
               
               r=0;p=0;
               qI[r]=i+1;qJ[r]=j;r++;
               mark2[i-1][j]=0;
               while(p<r)
               {
                  I=qI[p];J=qJ[p];p++;
                  if(J+1<W){if(tshGradient[I][J+1][0]==255&&mark2[I][J+1]==0){qI[r]=I;qJ[r]=J+1;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(J-1>=0){if(tshGradient[I][J-1][0]==255&&mark2[I][J-1]==0){qI[r]=I;qJ[r]=J-1;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(I+1<H){if(tshGradient[I+1][J][0]==255&&mark2[I+1][J]==0){qI[r]=I+1;qJ[r]=J;mark2[qI[r]][qJ[r]]=1;r++;}}
                  if(I-1>=0){if(tshGradient[I-1][J][0]==255&&mark2[I-1][J]==0){qI[r]=I-1;qJ[r]=J;mark2[qI[r]][qJ[r]]=1;r++;}}
               }
               
               int value=253;
               if(r>500)value=252;
               for(l=0;l<r;l++)
               for(k=0;k<3;k++)
               tshGradient[qI[l]][qJ[l]][k]=value;     
           }
        }   
     }

     max=0;
     for(i=0;i<numArea;i++)
     {
        count=0;
        for(j=0;j<numArea;j++)
        if(abs(area[i]-area[j])<50 && abs(R[i]-R[j])+abs(G[i]-G[j])+abs(B[i]-B[j])<50)count++;
        if(count>max){max=count;I=i;}
     }

     
     int tsh=70;
     for(i=0;i<numArea;i++)
     {
        count=0;
        for(j=0;j<numArea;j++)
        if(abs(R[i]-R[j])+abs(G[i]-G[j])+abs(B[i]-B[j])<tsh)count++;
        sameTeam[i]=count;
        
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
     
     for(i=0;i<numArea;i++)team[i]=-1;
     
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
        rectangle(myMat, 
                  Point(X[MST]-Warea[MST],Y[MST]-Harea[MST]), 
                  Point(X[MST]+Warea[MST],Y[MST]+Harea[MST]/2), 
                  WHITE,2,CV_AA,0);
                  
                  
        if(team[MST]==J)          
        rectangle(myMat, 
                  Point(X[MST]-Warea[MST],Y[MST]-Harea[MST]), 
                  Point(X[MST]+Warea[MST],Y[MST]+Harea[MST]/2), 
                  BLACK,2,CV_AA,0);
                  
        if(team[MST]==I)
           for(j=X[MST];j<=X[MST]+50;j++)
           {tshGradient[Y[MST]][j][0]=120;tshGradient[Y[MST]][j][1]=120;tshGradient[Y[MST]][j][2]=120;
           tshGradient[Y[MST]-1][j][0]=120;tshGradient[Y[MST]-1][j][1]=120;tshGradient[Y[MST]-1][j][2]=120;
           tshGradient[Y[MST]+1][j][0]=120;tshGradient[Y[MST]+1][j][1]=120;tshGradient[Y[MST]+1][j][2]=120;}
        else if(team[MST]==J)
           for(j=Y[MST];j>=Y[MST]-50;j--)
           {tshGradient[j][X[MST]][0]=120;tshGradient[j][X[MST]][1]=120;tshGradient[j][X[MST]][2]=120;
           tshGradient[j][X[MST]-1][0]=120;tshGradient[j][X[MST]-1][1]=120;tshGradient[j][X[MST]-1][2]=120;
           tshGradient[j][X[MST]+1][0]=120;tshGradient[j][X[MST]+1][1]=120;tshGradient[j][X[MST]+1][2]=120;}
     }   
}

void findBall()
{
     int i,j,k,l,r,score,numCircle;
     int H=matSize.height;
     int W=matSize.width;
     
     for(r=0;r<=20;r++)
     {
         numCircle=0;
         for(i=r;i<H/2-r;i++)
         for(j=r;j<W-r;j++)
         {
            score=0;        
            for(k=i-r;k<=i+r;k++)
            for(l=j-r;l<=j+r;l++)
            if((k-i)*(k-i)+(l-j)*(l-j)<=r*r){if(tshGradient[k][l][0]==255)score++;}
            //else{if(tshGradient[k][l][0]==0||tshGradient[k][l][0]==1)score++;}
               
            if((double)(score*100)>(double)(r*r*PI*90) && (abs(i-IBall)+abs(j-JBall)<50||IBall==-1))
            {
               ++numCircle;
               iBall = i;
               jBall = j;      
            }
         }
         if(numCircle==0){IBall=iBall;JBall=jBall;lastR=thisR;thisR=r-1;r=21;}
     }
     
     // paint gray color on the ball in bw image
     for(k=IBall-thisR;k<=IBall+thisR;k++)
     for(l=JBall-thisR;l<=JBall+thisR;l++)
     if((k-IBall)*(k-IBall)+(l-JBall)*(l-JBall)<=thisR*thisR)
     {tshGradient[k][l][0]=120;tshGradient[k][l][1]=120;tshGradient[k][l][2]=120;}

     if(thisR<ballRTrsh)findOffsidePlayer();
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

bool sortBySlope (Vec4i x,Vec4i y) { 
     return fabs(slope(x)) > fabs(slope(y)); 
} 

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

// first OR second point is in range
bool isInSameGroup(Vec4i x, Vec4i L) 
{
    double m,m0,m1,c;
    x[1]*=-1; x[3]*=-1; L[1]*=-1; L[3]*=-1;
    double x1=x[0], y1=x[1], x2=x[2], y2=x[3];
    
    double t = 40.0;
    const double slopeT = 0.3;
    
    m = (double)(L[3]-L[1])/(L[2]-L[0]+0.1);
    c = L[1] - m*L[0];
    m0 = (double)(y2-y1)/(x2-x1+0.1);
    m1 = m;
    if( fabs(m0) > 1 && fabs(m) > 1 ) { m0 = 1/m0; m1 = 1/m1; }
    if( m>1 ) t*=m;
    if( fabs(m0-m1) < slopeT
        && ( y1 > m * x1 + c - t && y1 < m * x1 + c + t
        && y2 > m * x2 + c - t && y2 < m * x2 + c + t
/*
value = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0)
if value > 0, p2 is on the left side of the line.
if value = 0, p2 is on the same line.
if value < 0, p2 is on the right side of the line.
*/
        || (L[2] - L[0])*(y1 - L[1]) - (x1 - L[0] + t)*(L[3] - L[1]) > 0 
        && (L[2] - L[0])*(y1 - L[1]) - (x1 - L[0] - t)*(L[3] - L[1]) < 0 
        || (L[2] - L[0])*(y2 - L[1]) - (x2 - L[0] + t)*(L[3] - L[1]) > 0 
        && (L[2] - L[0])*(y2 - L[1]) - (x2 - L[0] - t)*(L[3] - L[1]) < 0 )
        ) return true;
    return false;
}



void findVerticalLine() 
{
    sort(lines.begin(), lines.end(), myfunction);
    
    //find real field line
    for(int i = 0; i < lines.size(); i++) {
        Vec4i L = lines[i];
        bool sameGroup=false;
        for(int j = 0; j < groupLines.size(); j++) {
            if (isInSameGroup(L, groupLines[j][0])) {
                groupLines[j].push_back(L);
                sameGroup = true;
                break;
            }
        }
        if(!sameGroup) {
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
        for( int i = 0; i < points.size(); i++) {
            for( int j = i+1; j < points.size(); j++) {
                d2 = pow(points[i].x-points[j].x,2.0) + pow(points[i].y-points[j].y,2.0);
                if( d2 - maxd2 > 0.001 ) {
                    start = points[i];
                    end = points[j];
                    maxd2 = d2;
                }
            }
        }
        int x1 = start.x;
        int x2 = end.x;
        int y1 = start.y;
        int y2 = end.y;
        fieldLines.push_back(Vec4i(x1,y1,x2,y2));
    }
    
    for(int i = 0; i < fieldLines.size(); i++) {
        Vec4i L1 = fieldLines[i];
        Vec4i l = L1;
        if( l[0] < marginTrsh && l[2] < marginTrsh 
            || l[0] > width-marginTrsh && l[2] > width-marginTrsh 
            || l[1] < marginTrsh && l[3] < marginTrsh 
            || l[1] > height-marginTrsh && l[3] > height-marginTrsh ) {
            fieldLines.erase(fieldLines.begin()+i);
            i--;
            continue;
        }
        for(int j = i+1; j < fieldLines.size(); j++) {
            Vec4i L2 = fieldLines[j];
            if (isInSameGroup(L1, L2)) {
                fieldLines.erase(fieldLines.begin()+j);
                j--;
            }
        }
    }
    
    int x1,y1,x2,y2,c=0;
    Point2f o1,p1,o2,p2;
    double m;
    Vec4i l;

    // find longest vertical line 
    //sort(fieldLines.begin(), fieldLines.end(), sortBySlope);
    for(int i = 0; i < fieldLines.size(); i++) {
        l = fieldLines[i];            
        m = slope(l);
        if(m > vertSlopeTrsh || m < -vertSlopeTrsh) {
            vertLines.push_back(l);
            c++;
        }
        if(c==2) {
            prevVertLines = vertLines;
            break;
        }
    }

    if(vertLines.size()<2)
        if(prevVertLines.size()==2) vertLines = prevVertLines;
        else return;
        
    if( fabs(slope(vertLines[0])-slope(vertLines[1])) > 0.0001 ) {
        o1 = Point2f(vertLines[0][0],-vertLines[0][1]);
        p1 = Point2f(vertLines[0][2],-vertLines[0][3]);
        o2 = Point2f(vertLines[1][0],-vertLines[1][1]);
        p2 = Point2f(vertLines[1][2],-vertLines[1][3]);
        intersection(o1, p1, o2, p2, intPoint);
        if( //intPoint.x>0 && 
        intPoint.y>0 ) prevIntPoint = intPoint;
        else intPoint = prevIntPoint;
    } else parallel = true;

}

bool isNull(Vec4i v) {
    return v[0]==0 && v[1]==0 && v[2]==0 && v[3]==0;
}

void findOffenseDirection() {
    if( offenseDirection != NA ) return;
    
    for (int i=0; i<marginTrsh; i++) {
        for (int j=0; j<height; j++) {
            if ( tshGradient[j][i][0] == 255 ) {
                yLeftMargin = j;
                break;
            }
        }
        if( yLeftMargin > -1 ) break;
    }
    
    for (int i=width-1; i>width-marginTrsh; i--) {
        for (int j=0; j<height; j++) {
            if ( tshGradient[j][i][0] == 255 ) {
                yRightMargin = j;
                break;
            }
        }
        if( yRightMargin > -1 ) break;
    }
    
    if( yLeftMargin < yRightMargin )
        offenseDirection = RIGHT;
    else if( yLeftMargin > yRightMargin )
        offenseDirection = LEFT;
    
    
}

void printAndPaint() {
    if( offenseDirection != NA ) {
        if( offenseDirection == RIGHT ) cout<<"offense direction = RIGHT";
        else cout<<"offense direction = LEFT";
        cout<<endl;
    }
    printf("ball: %d %d %d\n",IBall,JBall,thisR);
    
    Vec4i l; int x1,x2,y1,y2;
    
    for(int i = 0; i < fieldLines.size(); i++) {
        l = fieldLines[i];
        x1=l[0]; y1=l[1]; x2=l[2]; y2=l[3];
        line(myMat,Point(x1,y1),Point(x2,y2),RED,1,CV_AA);
        circle(myMat, Point(x1,y1), 3, BLUE, -1);
        circle(myMat, Point(x2,y2), 3, BLUE, -1);
    }

    if(vertLines.size()<2 && prevVertLines.size()<2) 
        cout<<"find only "<<vertLines.size()<<" line\n";
    else { 
        cout<<"2 vertical lines: \n";
        for(int i = 0; i < vertLines.size(); i++) {
            l = vertLines[i];
            x1=l[0]; y1=l[1]; x2=l[2]; y2=l[3];
            line(myMat,Point(x1,y1),Point(x2,y2),PURPLE,1,CV_AA);
            circle(myMat, Point(x1,y1), 3, BLUE, -1);
            circle(myMat, Point(x2,y2), 3, BLUE, -1);
            cout<<x1<<" "<<y1<<" "<<x2<<" "<<y2<<endl;
        }
    }
    
    if(parallel) cout<<"(parallel) slopes: "<<slope(vertLines[0])<<" "<<slope(vertLines[1])<<endl;
    else cout<<"intersection point: "<<intPoint.x<<" "<<intPoint.y<<endl;
    if( yLeftMargin == yRightMargin )  cout<<"seriously!? yLeft = yRight margin!\n";
    //cout<<"-----------------------------------------------------\n";
    
    if(thisR>=ballRTrsh)
        putText( myMat, "ball", Point(JBall,IBall), FONT_HERSHEY_SIMPLEX, fontScale, BLACK , 2 );
    else
        putText( myMat, "ball not found", Point(JBall,IBall), FONT_HERSHEY_SIMPLEX, fontScale, YELLOW , 2 );
    
    int i=lastDFIndex;
    double m = ((double)(intPoint.y)+(double)(Y[i]))/((double)(intPoint.x)-(double)(X[i]));
    char output[50];
    snprintf(output,50,"%lf",m);
    putText(myMat, output, Point(X[i]+10, Y[i]), FONT_HERSHEY_SIMPLEX, fontScale, BLACK, 2 );
    putText(myMat, "last DF", Point(X[i]-50, Y[i]), FONT_HERSHEY_SIMPLEX, fontScale, WHITE, 2 );
    line(myMat,Point(intPoint.x,-intPoint.y),Point(X[i],Y[i]),PINK,1,CV_AA);
    
    for(int j=0; j<offsidePlyrIndexes.size(); j++) {
        i = offsidePlyrIndexes[j];
        double m = ((double)(intPoint.y)+(double)(Y[i]))/((double)(intPoint.x)-(double)(X[i]));
        snprintf(output,50,"%lf",m);
        putText(myMat, output, Point(X[i]+10, Y[i]), FONT_HERSHEY_SIMPLEX, fontScale, BLACK, 2 );
        putText(myMat, "offside", Point(X[i]-50, Y[i]), FONT_HERSHEY_SIMPLEX, fontScale, WHITE, 2 );
        line(myMat,Point(intPoint.x,-intPoint.y),Point(X[i],Y[i]),LIGHTBLUE,1,CV_AA);
    }
    
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
    
    // for collecting results
    int team[2][3]={0},lastDF[3]={0},offsidePlyr[3]={0},tp=0,fp=0,fn=0,countedFrame=0;
    int ball[3]={0},intPt=0,b=0,skip=false,same=false;
    int tp1=0,fp1=0,fn1=0,tp2=0,fp2=0,fn2=0,tp3=0,fp3=0,fn3=0,tp4=0,fp4=0,fn4=0,tp5=0,fp5=0,fn5=0,b1=0,b2=0;
    
    char const* avifile = "D:\\Google Drive\\SeniorProject\\video\\Lionel Messi is too intelligent for the offside trap-cropped.avi";
    //Lionel Messi is too intelligent for the offside trap-cropped.avi
    //Eibar vs Barcelona - cropped.mp4"
    //Barcelona canceled goal vs Manchester City (WASN'T OFFSIDE) - cropped.mp4
    //Arnautovic OFFSIDE GOAL!! vs Liverpool - Semifinal - Capital one Cup 2016 - cut.mp4";
     //xavi;
    
    capture = cvCaptureFromAVI(avifile);
    if(!capture)throw "Error when reading steam_avi";
    cvNamedWindow( "color", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
          
    
    for(framePointer=0 ; framePointer<FRAMES_WANTED ; ++framePointer)
    {
        frame = cvQueryFrame( capture );
        if(!frame)break;
        
        //cout<<framePointer<<endl;
        if(framePointer < 24) continue;
        
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
        if(gaussian[i][j][1]>=grassTrsh && gaussian[i][j][1]>gaussian[i][j][2] && gaussian[i][j][2]>gaussian[i][j][0])
        //{grass[i][j][0]=255;grass[i][j][1]=255;grass[i][j][2]=255;}
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
        	
		
        for(i=0; i<height ;i++)
        for(j=0; j<width ;j++)
        if(gradient[i][j][0]>tshGD && gradient[i][j][1]>tshGD && gradient[i][j][2]>tshGD) 
        {tshGradient[i][j][0]=255;tshGradient[i][j][1]=255;tshGradient[i][j][2]=255;}
        else {tshGradient[i][j][0]=0;tshGradient[i][j][1]=0;tshGradient[i][j][2]=0;}
        

        Mat bw_image (height, width, CV_8UC1);
		for(i=0; i<height ;i++)
		for(j=0; j<width ;j++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][0];
        
        HoughLinesP(bw_image, lines, 5, CV_PI/180, 100, 100, 10 );
        
        findVerticalLine();
        
        findOffenseDirection();
        
        findPlayer();
        
        findBall(); 
        
        printAndPaint();
        
        for(i=0; i<matSize.height ;i++)
        for(j=0; j<matSize.width ;j++)
        for(k=0; k<3 ;k++)
        bw_image.at<uchar>(i,j) = tshGradient[i][j][k];
        imshow("color", myMat);
        imshow("bw", bw_image);
        
        lines.clear();
        fieldLines.clear();
        vertLines.clear();
        groupLines.clear();
        parallel = false;
        offsidePlyrIndexes.clear(); 
        
        // collect results
        //skip = false; 
        same = false; b=0; tp=0;fp=0;fn=0;
        cout<<"::::::::::::::::::RESULTS::::::::::::::::::"<<endl;
        if(framePointer==0) skip = true;
        if(framePointer==1) skip = false;
        if(framePointer>90 && !skip) { 
                            cout<<"*** skip? 1/0 ***: "; cin>>skip; 
                            }
        
        if(!skip) {
            countedFrame++;
            cout<<"same? (1/0): "; cin>>same;
            if(same) {
                team[0][0]+=tp1; team[0][1]+=fp1; team[0][2]+=fn1;
                team[1][0]+=tp2; team[1][1]+=fp2; team[1][2]+=fn2;
                intPt+=b1;
                ball[0]+=tp5; ball[1]+=tp5; ball[2]+=tp5;
                if(!b2) {
                    lastDF[0]+=tp3; lastDF[1]+=fp3; lastDF[2]+=fn3;
                    offsidePlyr[0]+=tp4; offsidePlyr[1]+=fp4; offsidePlyr[2]+=fn4;
                }
            }
            else {
                cout<<"team 1 (Tp Fp Fn): "; cin>>tp>>fp>>fn; 
                team[0][0]+=tp; team[0][1]+=fp; team[0][2]+=fn;
                tp1=tp; fp1=fp; fn1=fn;
                cout<<"team 2 (Tp Fp Fn): "; cin>>tp>>fp>>fn;
                team[1][0]+=tp; team[1][1]+=fp; team[1][2]+=fn;
                tp2=tp; fp2=fp; fn2=fn;
                cout<<"ball (Tp Fp Fn): "; cin>>tp>>fp>>fn;
                ball[0]+=tp; ball[1]+=fp; ball[2]+=fn;
                tp5=tp; fp5=fp; fn5=fn;
                cout<<"intersection point? (1/0): "; cin>>b;
                intPt+=b; b1=b;
                cout<<"ball found? (1/0): "; cin>>b; b2=b;
                if(!b) {
                    cout<<"last defender (Tp Fp Fn): "; cin>>tp>>fp>>fn;
                    lastDF[0]+=tp; lastDF[1]+=fp; lastDF[2]+=fn;
                    tp3=tp; fp3=fp; fn3=fn;
                    cout<<"offside player(s) (Tp Fp Fn): "; cin>>tp>>fp>>fn;
                    offsidePlyr[0]+=tp; offsidePlyr[1]+=fp; offsidePlyr[2]+=fn;
                    tp4=tp; fp4=fp; fn4=fn;
                }
            }
        }
        
        key = cvWaitKey(63);
        if(key==27)break;
    }
    
    // show results
    cout<<"\n\n\n=================ALL RESULTS=================\n\n";
    cout<<"counted frames: "<<countedFrame<<"\n\n";
    
    cout.precision(3);
    
    tp=team[0][0]; fp=team[0][1]; fn=team[0][2];
    cout<<"team 1 (Tp Fp Fn): "<<tp<<" "<<fp<<" "<<fn<<endl;
    cout<<"P = "<<fixed<<(double)tp/(tp+fp)<<", R = "<<fixed<<(double)tp/(tp+fn)<<"\n\n";
    
    tp=team[1][0]; fp=team[1][1]; fn=team[1][2];
    cout<<"team 2 (Tp Fp Fn): "<<tp<<" "<<fp<<" "<<fn<<endl;
    cout<<"P = "<<fixed<<(double)tp/(tp+fp)<<", R = "<<fixed<<(double)tp/(tp+fn)<<"\n\n";
    
    tp=ball[0]; fp=ball[1]; fn=ball[2];
    cout<<"ball (Tp Fp Fn): "<<tp<<" "<<fp<<" "<<fn<<endl;
    cout<<"P = "<<fixed<<(double)tp/(tp+fp)<<", R = "<<fixed<<(double)tp/(tp+fn)<<"\n\n";
    
    cout<<"intersection point (T F): "<<intPt<<" "<<countedFrame-intPt<<endl;
    cout<<"accuracy = "<<fixed<<(double)intPt*100/countedFrame<<"%"<<"\n\n";
    
    tp=lastDF[0]; fp=lastDF[1]; fn=lastDF[2];
    cout<<"last defender (Tp Fp Fn): "<<tp<<" "<<fp<<" "<<fn<<endl;
    cout<<"P = "<<fixed<<(double)tp/(tp+fp)<<", R = "<<fixed<<(double)tp/(tp+fn)<<"\n\n";
    
    tp=offsidePlyr[0]; fp=offsidePlyr[1]; fn=offsidePlyr[2];
    cout<<"offside player(s) (Tp Fp Fn): "<<tp<<" "<<fp<<" "<<fn<<endl;
    cout<<"P = "<<fixed<<(double)tp/(tp+fp)<<", R = "<<fixed<<(double)tp/(tp+fn)<<"\n\n";
   
   
    cvWaitKey(0); 
    cvDestroyWindow("bw");
    cvDestroyWindow( "color" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}


