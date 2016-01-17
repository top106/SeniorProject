#include <stdio.h>
#include <cv.h>
#include "highgui.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

//using namespace cv

int main(int argc, char** argv)
{
    CvCapture* capture=0;
    IplImage* frame=0;
    IplImage* frame_gray=0;
    IplImage* frame_bw=0;
    uchar b, g, r;

    
    FILE *fp;
    char const* textfile = "D:\\Google Drive\\SeniorProject\\test.txt";

    fp = fopen(textfile, "w+");
    if(fp == NULL){
        printf("\nNão encontrei arquivo\n");
        //system("pause");
        exit(EXIT_FAILURE);
    }
    char const* avifile = "D:\\Google Drive\\SeniorProject\\sample video.avi";
    capture = cvCaptureFromAVI(avifile); // read AVI video    
    if( !capture )
        throw "Error when reading steam_avi";

    char c;
    cvNamedWindow( "small", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "gray", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "bw", CV_WINDOW_AUTOSIZE );
    cvNamedWindow( "w", 1 );
    const int reduce = 1;
    int width, height, widthStep, channels;
    int k=0;
    for( ; k<10 ; )
    {
        frame = cvQueryFrame( capture );
        if(!frame)
            break;
        cvShowImage("w", frame);
        
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

        for(int i = 0; i < height; i++){
            for(int j = 0; j < width; j++){
                //b = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 0]; // b
                //g = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 1]; // g
                //r = ((uchar *)(tmpsize->imageData + i*widthStep))[j*channels + 2]; // r
                //B G R
                //fprintf(fp, "%i %i %i\n", b, g, r);
            }
        }
        
        for(i=2; i<S.height-2 ;i++)
			for(j=2; j<S.width-2 ;j++)
			for(k=0; k<3 ;k++)
			p[i][j][k]=( 15*p[i][j][k]
					     + 12*( p[i-1][j][k] + p[i+1][j][k] + p[i][j-1][k] + p[i][j+1][k] )
					     + 9*(p[i-1][j-1][k]+p[i+1][j-1][k]+p[i-1][j+1][k]+p[i+1][j+1][k])
					     + 5*(p[i-2][j][k]+p[i+2][j][k]+p[i][j-2][k]+p[i][j+2][k])
					     + 2*(p[i-2][j-2][k]+p[i+2][j-2][k]+p[i-2][j+2][k]+p[i+2][j+2][k])
					     + 4*(p[i-2][j-1][k]+p[i-2][j+1][k]+p[i+2][j-1][k]+p[i+2][j+1][k]+p[i-1][j-2][k]+p[i+1][j-2][k]+p[i-1][j+2][k]+p[i+1][j+2][k])
					     ) / 159;
			
        for(i=1; i<S.height-1 ;i++)
			for(j=1; j<S.width-1 ;j++)
                for(k=0; k<3 ;k++)
                    d[i][j][k] = Math.sqrt(Math.pow(p[i-1][j+1][k]+p[i+1][j+1][k]+p[i][j+1][k]*2-p[i-1][j-1][k]-p[i+1][j-1][k]-p[i][j-1][k]*2,2)+
					     Math.pow(p[i+1][j+1][k]+p[i+1][j-1][k]+p[i+1][j][k]*2-p[i-1][j-1][k]-p[i-1][j+1][k]-p[i-1][j][k]*2,2));
/////////////////////////////////////////////
			
			for(i=0; i<S.height ;i++)
			for(j=0; j<S.width ;j++)
			if(d[i][j][0]>10&&d[i][j][1]>10&&d[i][j][2]>10){d[i][j][0]=255;d[i][j][1]=255;d[i][j][2]=255;}
			else {d[i][j][0]=0;d[i][j][1]=0;d[i][j][2]=0;}
        
        ++k;

        fprintf(fp, "\n------NEXT FRAME------\n\n");

        cvShowImage("small", frame);

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
    cvDestroyWindow("w");
    cvDestroyWindow("bw");
    cvDestroyWindow( "gray" );
    cvDestroyWindow( "small" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}
