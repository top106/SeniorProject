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
uchar h[365][10000];
        
int main(int argc, char** argv)
{
    CvCapture* capture=0;
    IplImage* frame=0;
    IplImage* frame_gray=0;
    IplImage* frame_bw=0;
    uchar b, g, r;
    int i,j,k,f;
    
    FILE *fp;
    char const* textfile = "D:\\Study\\Project\\Code\\test.txt";

    fp = fopen(textfile, "w+");
    if(fp == NULL){
        printf("\nNão encontrei arquivo\n");
        //system("pause");
        exit(EXIT_FAILURE);
    }
    char const* avifile = "D:\\Study\\Project\\Code\\2.mp4";
    capture = cvCaptureFromAVI(avifile); // read AVI video    
    if( !capture )
        throw "Error when reading steam_avi";

    char c;
    cvNamedWindow( "small", CV_WINDOW_AUTOSIZE );
    //cvNamedWindow( "gray", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
    //cvNamedWindow( "w", 1 );
    
    const int reduce = 2;
    
    int width, height, widthStep, channels;
    for(f=0 ; f<100 ; ++f)
    {
        frame = cvQueryFrame( capture );
        if(!frame)break;
        
        //cvShowImage("w", frame);
        
        CvSize size = cvSize(frame->width / reduce, frame->height / reduce);
        IplImage* tmpsize = cvCreateImage(size, frame->depth, frame->nChannels);
        cvResize(frame, tmpsize, CV_INTER_LINEAR);
        
//        frame_gray = cvCreateImage(cvGetSize(tmpsize),IPL_DEPTH_8U,1);
//        cvCvtColor(tmpsize, frame_gray, CV_BGR2GRAY);
//        cvShowImage("gray", frame_gray);
//        
//        frame_bw = cvCreateImage(cvGetSize(frame_gray),IPL_DEPTH_8U,1);
//        cvThreshold(frame_gray, frame_bw, 10, 255, CV_THRESH_BINARY_INV | CV_THRESH_OTSU);
//        cvShowImage("bw", frame_bw);


        
        height = tmpsize->height;
        width = tmpsize->width;
        widthStep = tmpsize->widthStep;
        channels = tmpsize->nChannels;
        
        cvShowImage("small", tmpsize);

        /*for(i = 0; i < height; i++){
            for(j = 0; j < width; j++){
                //b = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 0]; // b
                //g = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 1]; // g
                //r = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 2]; // r
                //B G R
                //fprintf(fp, "%i %i %i\n", b, g, r);
            }
        }*/
        Mat A (tmpsize);
	
		Size S = A.size();
        uchar p[height][width][3];
	
		
			for(i=0; i<S.height ;i++)
			for(j=0; j<S.width ;j++)
			for(k=0; k<3 ;k++)
			p[i][j][k]=A.at<cv::Vec3b>(i,j)[k];
		
        
        uchar d[height][width][3];
        uchar g[height][width][3];
        uchar D[height][width][3];
        uchar t[height][width][3];

        
        for(i=2; i<S.height-2 ;i++) {
			for(j=2; j<S.width-2 ;j++) {
                for(k=0; k<3 ;k++) {
                    p[i][j][k] = ( 41*p[i][j][k]
					     + 26*( p[i-1][j][k] + p[i+1][j][k] + p[i][j-1][k] + p[i][j+1][k] )
					     + 16*(p[i-1][j-1][k]+p[i+1][j-1][k]+p[i-1][j+1][k]+p[i+1][j+1][k])
					     + 7*(p[i-2][j][k]+p[i+2][j][k]+p[i][j-2][k]+p[i][j+2][k])
					     + 1*(p[i-2][j-2][k]+p[i+2][j-2][k]+p[i-2][j+2][k]+p[i+2][j+2][k])
					     + 4*( p[i-2][j-1][k]
                               +p[i-2][j+1][k]
                               +p[i+2][j-1][k]
                               +p[i+2][j+1][k]
                               +p[i-1][j-2][k]
                               +p[i+1][j-2][k]
                               +p[i-1][j+2][k]
                               +p[i+1][j+2][k] )
					     ) / 273;
                }
            }
        }
        
        for(i=0; i<S.height ;i++) 
        for(j=0; j<S.width ;j++) 
        if(p[i][j][1]>p[i][j][2]&&p[i][j][2]>p[i][j][0])
        {g[i][j][0]=p[i][j][0];g[i][j][1]=p[i][j][1];g[i][j][2]=p[i][j][2];}
        else {g[i][j][0]=0;g[i][j][1]=0;g[i][j][2]=0;}
			
        for(i=1; i<S.height-1 ;i++)
			for(j=1; j<S.width-1 ;j++)
                for(k=0; k<3 ;k++)
                    d[i][j][k] = (uchar) sqrt ( (double) ( 
                               pow ( 
                               (double) (p[i-1][j+1][k]+p[i+1][j+1][k]+p[i][j+1][k]*2-p[i-1][j-1][k]-p[i+1][j-1][k]-p[i][j-1][k]*2) 
                               , 2.0)
                               +  pow ( 
                               (double) (p[i+1][j+1][k]+p[i+1][j-1][k]+p[i+1][j][k]*2-p[i-1][j-1][k]-p[i-1][j+1][k]-p[i-1][j][k]*2)
                               , 2.0) 
					           )
                               );
/////////////////////////////////////////////
			
        const int tsh = 30;
        
        for(i=0; i<S.height ;i++)
			for(j=0; j<S.width ;j++)
                if( d[i][j][0]>tsh && d[i][j][1]>tsh && d[i][j][2]>tsh ) {
                    d[i][j][0]=255;
                    d[i][j][1]=255;
                    d[i][j][2]=255;
                }
			    else {
                    d[i][j][0]=0;
                    d[i][j][1]=0;
                    d[i][j][2]=0;
  
                }
/*   
        int a,e,I,J,r,R;        
        
        for(i=0; i<S.height ;i++)
        for(j=0; j<S.width ;j++)
        t[i][j][0]=0;
        
        for(i=0; i<S.height ;i++)
        for(j=0; j<S.width ;j++)
        {
            for(r=0; ;r++)
            {
               e = 0;
               if(d[i][j-r][0]==255)e++;
               if(d[i][j+r][0]==255)e++;
               if(d[i-r][j][0]==255)e++;
               if(d[i+r][j][0]==255)e++;
               
               if(d[i-(int)(r/Root2)][j-(int)(r/Root2)][0]==255)e++;
               if(d[i-(int)(r/Root2)][j+(int)(r/Root2)][0]==255)e++;
               if(d[i+(int)(r/Root2)][j-(int)(r/Root2)][0]==255)e++;
               if(d[i+(int)(r/Root2)][j+(int)(r/Root2)][0]==255)e++;
               
               if(d[i-(int)(r/2)][j-(int)(r*Root3/2)][0]==255)e++;
               if(d[i-(int)(r/2)][j+(int)(r*Root3/2)][0]==255)e++;
               if(d[i+(int)(r/2)][j-(int)(r*Root3/2)][0]==255)e++;
               if(d[i+(int)(r/2)][j+(int)(r*Root3/2)][0]==255)e++;
               
               if(d[i-(int)(r*Root3/2)][j-(int)(r/2)][0]==255)e++;
               if(d[i-(int)(r*Root3/2)][j+(int)(r/2)][0]==255)e++;
               if(d[i+(int)(r*Root3/2)][j-(int)(r/2)][0]==255)e++;
               if(d[i+(int)(r*Root3/2)][j+(int)(r/2)][0]==255)e++;
               
               //printf("%d %d %d %d %d\n",f,I,J,r,e);
               if(e<=10)break;
            }
            if(r>=3)
            for(R=0;R<=r;R++)
            {
               t[i][j-R][0]=255;
               t[i][j+R][0]=255;
               t[i-R][j][0]=255;
               t[i+R][j][0]=255;
               
               t[i-(int)(R/Root2)][j-(int)(R/Root2)][0]=255;
               t[i-(int)(R/Root2)][j+(int)(R/Root2)][0]=255;
               t[i+(int)(R/Root2)][j-(int)(R/Root2)][0]=255;
               t[i+(int)(R/Root2)][j+(int)(R/Root2)][0]=255;
               
               t[i-(int)(r/2)][j-(int)(r*Root3/2)][0]=255;
               t[i-(int)(r/2)][j+(int)(r*Root3/2)][0]=255;
               t[i+(int)(r/2)][j-(int)(r*Root3/2)][0]=255;
               t[i+(int)(r/2)][j+(int)(r*Root3/2)][0]=255;
               
               t[i-(int)(r*Root3/2)][j-(int)(r/2)][0]=255;
               t[i-(int)(r*Root3/2)][j+(int)(r/2)][0]=255;
               t[i+(int)(r*Root3/2)][j-(int)(r/2)][0]=255;
               t[i+(int)(r*Root3/2)][j+(int)(r/2)][0]=255;
            }
            //printf("%d %d %d %d\n",f,i,j,r);
        }
*/
////////////////////////////////////////
        int max=0,l;                                
        for(i=0; i<S.height ;i++)
			for(j=0; j<S.width ;j++)
			  if(d[i][j][0]==255)
			  {
				for(k=0; k<=360 ;k++)
				{
					l=(int) (j*cos(k*PI/180)+i*sin(k*PI/180)+S.height+S.width);
					h[k][l]++;
					if(l>max)max=l;
				}
			  }
			
			for(k=0; k<=360 ;k++)
			for(l=0; l<=max ;l++)
			  if(h[k][l]>S.height/2)
			  {
                //printf("%d %d\n",k,l);break;
				for(i=0; i<S.height ;i++)
				for(j=0; j<S.width ;j++)
				if(d[i][j][0]==255 && l==(int) (j*cos(k*PI/180)+i*sin(k*PI/180)+S.height+S.width))
				{D[i][j][0]=255;D[i][j][1]=255;D[i][j][2]=255;}
			    c++;
    			break;
			  }	
////////////////////////////////////////
        //fprintf(fp, "\n------NEXT FRAME------\n\n");
        Mat bw_image (height, width, CV_8UC1);
		S = bw_image.size();
		
		
			for(i=0; i<S.height ;i++)
			for(j=0; j<S.width ;j++)
			for(k=0; k<3 ;k++)
			bw_image.at<uchar>(i,j) = g[i][j][k];
        
        

        imshow("bw", bw_image);

        c = cvWaitKey(33);
        if( c == 27 ) break;
    }
    
    /*printf("Enter frame number: ");
    scanf("%d", &k);
    for(int i = 0; i < tmpsizeH; i++){
        for(int j = 0; j < tmpsizeW; j++){
            printf("B=%i, G=%i, R=%i\n", B[k][i][j], G[k][i][j], R[k][i][j]);
        }
        printf("\n");
    }*/
/*
    printf("press enter to continue: ");
    scanf("%c", &c);
     Seek to the beginning of the file 
    fseek(fp, SEEK_SET, 0);

    Read and display data 
    const int MAX_LINE = 100;
    char buffer[MAX_LINE];
    for(int i=0; i<5; i++){
        fgets(buffer, MAX_LINE, fp);
        printf("%d %s", i+1, buffer);
    }*/
    
    fclose(fp);
    cvWaitKey(0); // key press to close window
    //cvDestroyWindow("w");
    cvDestroyWindow("bw");
    //cvDestroyWindow( "gray" );
    cvDestroyWindow( "small" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}
