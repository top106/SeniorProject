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
    uchar b, g, r;
    uchar ***B = (uchar***)malloc(10*sizeof(uchar**));
    uchar ***G = (uchar***)malloc(10*sizeof(uchar**));
    uchar ***R = (uchar***)malloc(10*sizeof(uchar**));

    
    FILE *fp;
    char const* textfile = "D:\\Google Drive\\SeniorProject\\test.txt";
    //fp = fopen(textfile, "w");
    /// Open file for both reading and writing
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
    cvNamedWindow( "w", 1 );
    
    int k=0, tmpsizeW, tmpsizeH;
    for( ; k<10 ; )
    {
        frame = cvQueryFrame( capture );
        if(!frame)
            break;
        cvShowImage("w", frame);
        
        int reduce = 20;
        CvSize size = cvSize(frame->width / reduce, frame->height / reduce);
        IplImage* tmpsize = cvCreateImage(size, frame->depth, frame->nChannels);
        tmpsizeW = tmpsize->width; tmpsizeH = tmpsize->height;
        cvResize(frame, tmpsize, CV_INTER_LINEAR);

        B[k] = (uchar**)malloc(tmpsize->height*sizeof(uchar*));
        G[k] = (uchar**)malloc(tmpsize->height*sizeof(uchar*));
        R[k] = (uchar**)malloc(tmpsize->height*sizeof(uchar*));


        for(int i = 0; i < tmpsize->height; i++){
            B[k][i] = (uchar*)malloc(tmpsize->width*sizeof(uchar));
            G[k][i] = (uchar*)malloc(tmpsize->width*sizeof(uchar));
            R[k][i] = (uchar*)malloc(tmpsize->width*sizeof(uchar));

            for(int j = 0; j < tmpsize->width; j++){
                b = ((uchar *)(tmpsize->imageData + i*tmpsize->widthStep))[j*tmpsize->nChannels + 0]; // b
                g = ((uchar *)(tmpsize->imageData + i*tmpsize->widthStep))[j*tmpsize->nChannels + 1]; // g
                r = ((uchar *)(tmpsize->imageData + i*tmpsize->widthStep))[j*tmpsize->nChannels + 2]; // r
                B[k][i][j] = b;
                G[k][i][j] = g;
                R[k][i][j] = r;

                fprintf(fp, "B=%i, G=%i, R=%i\n", b, g, r);
            }
        }
        ++k;

        fprintf(fp, "\n------NEXT FRAME------\n\n");

        cvShowImage("small", tmpsize);

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

    printf("press enter to continue: ");
    scanf("%c", &c);
    /* Seek to the beginning of the file */
    fseek(fp, SEEK_SET, 0);

    /* Read and display data */
    const int MAX_LINE = 100;
    char buffer[MAX_LINE];
    for(int i=0; i<5; i++){
        fgets(buffer, MAX_LINE, fp);
        printf("%d %s", i+1, buffer);
    }
    
    fclose(fp);
    cvWaitKey(0); // key press to close window
    cvDestroyWindow("w");
    cvDestroyWindow( "small" );
    cvReleaseImage(&frame);
    cvReleaseCapture( &capture );
}
