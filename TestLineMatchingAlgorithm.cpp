/*IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.

 By downloading, copying, installing or using the software you agree to this license.
 If you do not agree to this license, do not download, install,
 copy or use the software.


                          License Agreement
               For Open Source Computer Vision Library

Copyright (C) 2011-2012, Lilian Zhang, all rights reserved.
Copyright (C) 2013, Manuele Tamburrano, Stefano Fabri, all rights reserved.
Third party copyrights are property of their respective owners.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

  * Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.

  * Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.

  * The name of the copyright holders may not be used to endorse or promote products
    derived from this software without specific prior written permission.

This software is provided by the copyright holders and contributors "as is" and
any express or implied warranties, including, but not limited to, the implied
warranties of merchantability and fitness for a particular purpose are disclaimed.
In no event shall the Intel Corporation or contributors be liable for any direct,
indirect, incidental, special, exemplary, or consequential damages
(including, but not limited to, procurement of substitute goods or services;
loss of use, data, or profits; or business interruption) however caused
and on any theory of liability, whether in contract, strict liability,
or tort (including negligence or otherwise) arising in any way out of
the use of this software, even if advised of the possibility of such damage.
*/


#include <math.h>
#include <time.h>
#include <fstream>


#include "LineDescriptor.hh"


using namespace std;


void usage(int argc, char** argv){
	cout<<"Usage: "<<argv[0]<<"  image1.png"<<"  image2.png"<<endl;
}

unsigned int factorial(unsigned int x)
{
    unsigned int value = 1;

    for(unsigned int i = 2; i <= x; i++)
    {
        value = value * i;
    }

    return value;
}

void calcEasyCombination(unsigned int n)
{   
    int num_combinations = factorial(n)/(factorial(n-2)*2);
    
    unsigned int results [num_combinations][2];
    unsigned int a, b, index = 0;
    for(int i = 0; i<n-1; i++)
    {
        a = i;
        for(int j = i+1; j<n; j++)
        {
            b = j;
            results[index][0]=a;
            results[index][1]=b;
            std::cout<<index<<std::endl;
            index++;
        }
    }
    
    std::cout<<"unsigned int combinations [num_combinations][2] = {";
    for(int i = 0; i<num_combinations; i++)
    {
        if(i != num_combinations -1)
            std::cout<<"{"<<results[i][0]<<","<<results[i][1]<<"},";
        else
            std::cout<<"{"<<results[i][0]<<","<<results[i][1]<<"}};"<<std::endl;
    }

    //unsigned int combinations [36][2] = {{0,1},{0,2},{0,3},{0,4},{0,5},{0,6},{0,7},{0,8},{1,2},{1,3},{1,4},{1,5},{1,6},{1,7},{1,8},{2,3},{2,4},{2,5},{2,6},{2,7},{2,8},{3,4},{3,5},{3,6},{3,7},{3,8},{4,5},{4,6},{4,7},{4,8},{5,6},{5,7},{5,8},{6,7},{6,8},{7,8}};
    
}


int main(int argc, char** argv)
{
	int ret = -1;
	if(argc<3){
		usage(argc,argv);
		return ret;
	}
	
	//calcEasyCombination(9);
	
  //load first image from file
	std::string imageName1(argv[1]);
	cv::Mat leftImage;
    leftImage = imread(imageName1, cv::IMREAD_GRAYSCALE);   // Read the file
    if(! leftImage.data )                              // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }

  //load second image from file
	std::string imageName2(argv[2]);
	
    cv::Mat rightImage;
    rightImage = imread(imageName2, cv::IMREAD_GRAYSCALE);   // Read the file
    if(! rightImage.data )                              // Check for invalid input
    {
        cout <<  "Could not open or find the image" << std::endl ;
        return -1;
    }

	unsigned int imageWidth  = leftImage.cols;
	unsigned int imageHeight = leftImage.rows;

	srand((unsigned)time(0));
	int lowest=100, highest=255;
	int range=(highest-lowest)+1;
	unsigned int r, g, b; //the color of lines

	//initial variables
	cv::Mat leftColorImage(leftImage.size(), CV_8UC3);
	cv::Mat rightColorImage(rightImage.size(), CV_8UC3);
	
	cvtColor( leftImage, leftColorImage,  cv::COLOR_GRAY2RGB );
	cvtColor( rightImage, rightColorImage, cv::COLOR_GRAY2RGB );


  ///////////####################################################################
  ///////////####################################################################
	//extract lines, compute their descriptors and match lines
	LineDescriptor lineDesc;
//	PairwiseLineMatching lineMatch;

	ScaleLines   linesInLeft;
	ScaleLines   linesInRight;
	std::vector<unsigned int> matchResult;


//	for(int h = 0; h<10; h++)
//	{
////	    cv::imshow("PRIMA",leftImage);
//	    //cv::Mat m = rightImage.clone();
//	    lineDesc.GetLineDescriptor(rightImage,linesInLeft);
////	    cv::imshow("DOPO",leftImage);
////	    cv::waitKey();
//	}
//	return 1;
//	

	cv::Mat binaryDescriptors;
	
	std::cout<<"leftImage.cols: "<<leftImage.cols<<"leftImage.rows: "<<leftImage.rows<<std::endl;
	
	lineDesc.GetLineDescriptor(leftImage,linesInLeft);
	lineDesc.GetLineDescriptor(rightImage,linesInRight);
	lineDesc.GetLineBinaryDescriptor(binaryDescriptors, linesInLeft);
	//TODO remove BIAS dependecies in PairwiseMatching
//	lineMatch.LineMatching(linesInLeft,linesInRight,matchResult);

  ///////////####################################################################
  ///////////####################################################################
	//draw  extracted lines into images
	cv::Point startPoint;
	cv::Point endPoint;
	cv::Point point;

	for(unsigned int i=0; i<linesInLeft.size(); i++){
		r = lowest+int(rand()%range);
		g = lowest+int(rand()%range);
		b = lowest+int(rand()%range);
		startPoint = cv::Point(int(linesInLeft[i][0].startPointX),int(linesInLeft[i][0].startPointY));
		endPoint   = cv::Point(int(linesInLeft[i][0].endPointX),  int(linesInLeft[i][0].endPointY));
		cv::line( leftColorImage,startPoint,endPoint,cv::Scalar(r,g,b));
//		char buf[10];
//		sprintf( buf,   "%d ",  i);
//		if(i%2){
//			point = cvPoint(round(0.75*startPoint.x+0.25*endPoint.x),round(0.75*startPoint.y+0.25*endPoint.y));
//			cvPutText(cvLeftColorImage,buf,point,&font,CV_RGB(r,g,b));
//		}else{
//			point = cvPoint(round(0.25*startPoint.x+0.75*endPoint.x),round(0.25*startPoint.y+0.75*endPoint.y));
//			cvPutText(cvLeftColorImage,buf,point,&font,CV_RGB(r,g,b));
//		}
	}
	for(unsigned int i=0; i<linesInRight.size(); i++){
		r = lowest+int(rand()%range);
		g = lowest+int(rand()%range);
		b = lowest+int(rand()%range);
		startPoint = cv::Point(int(linesInRight[i][0].startPointX),int(linesInRight[i][0].startPointY));
		endPoint   = cv::Point(int(linesInRight[i][0].endPointX),  int(linesInRight[i][0].endPointY));
		cv::line( rightColorImage,startPoint,endPoint,cv::Scalar(r,g,b));
//		char buf[10];
//		sprintf( buf,   "%d ",  i);
//		if(i%2){
//			point = cvPoint(round(0.75*startPoint.x+0.25*endPoint.x),round(0.75*startPoint.y+0.25*endPoint.y));
//			cvPutText(cvRightColorImage,buf,point,&font,CV_RGB(r,g,b));
//		}else{
//			point = cvPoint(round(0.25*startPoint.x+0.75*endPoint.x),round(0.25*startPoint.y+0.75*endPoint.y));
//			cvPutText(cvRightColorImage,buf,point,&font,CV_RGB(r,g,b));
//		}
	}
	imwrite("LinesInImage1.png",leftColorImage);
	imwrite("LinesInImage2.png",rightColorImage);
	//TODO enable after BIAS dependecies in PairwiseMatching have been removed
//  ///////////####################################################################
//  ///////////####################################################################
//  //store the matching results of the first and second images into a single image
//	double ww1,ww2;
//	int lineIDLeft;
//	int lineIDRight;
//	cvCvtColor( cvLeftImage, cvLeftColorImage,  CV_GRAY2RGB );
//	cvCvtColor( cvRightImage,cvRightColorImage, CV_GRAY2RGB );
//	int lowest1=0, highest1=255;
//	int range1=(highest1-lowest1)+1;
//	std::vector<unsigned int> r1(matchResult.size()/2), g1(matchResult.size()/2), b1(matchResult.size()/2); //the color of lines
//	for(unsigned int pair=0; pair<matchResult.size()/2;pair++){
//		r1[pair] = lowest1+int(rand()%range1);
//		g1[pair] = lowest1+int(rand()%range1);
//		b1[pair] = 255 - r1[pair];
//		ww1 = 0.2*(rand()%5);
//		ww2 = 1 - ww1;
//		char buf[10];
//		sprintf( buf,   "%d ",  pair);
//		lineIDLeft = matchResult[2*pair];
//		lineIDRight= matchResult[2*pair+1];
//		startPoint = cvPoint(int(linesInLeft[lineIDLeft][0].startPointX),int(linesInLeft[lineIDLeft][0].startPointY));
//		endPoint   = cvPoint(int(linesInLeft[lineIDLeft][0].endPointX),  int(linesInLeft[lineIDLeft][0].endPointY));
//		cvLine( cvLeftColorImage,startPoint,endPoint,CV_RGB(r1[pair],g1[pair],b1[pair]),4, CV_AA);
//		startPoint = cvPoint(int(linesInRight[lineIDRight][0].startPointX),int(linesInRight[lineIDRight][0].startPointY));
//		endPoint   = cvPoint(int(linesInRight[lineIDRight][0].endPointX),  int(linesInRight[lineIDRight][0].endPointY));
//		cvLine( cvRightColorImage,startPoint,endPoint,CV_RGB(r1[pair],g1[pair],b1[pair]),4, CV_AA);
//	}
//
//	IplImage *cvResultColorImage1 = cvCreateImage(cvSize(imageWidth*2,imageHeight),IPL_DEPTH_8U, 3);
//	IplImage *cvResultColorImage2 = cvCreateImage(cvSize(imageWidth*2,imageHeight),IPL_DEPTH_8U, 3);
//	IplImage *cvResultColorImage = cvCreateImage(cvSize(imageWidth*2,imageHeight),IPL_DEPTH_8U, 3);
//	cvSetImageROI(cvResultColorImage1, cvRect(0, 0, imageWidth-1, imageHeight-1));
//	cvResize(cvLeftColorImage, cvResultColorImage1);
//	cvResetImageROI(cvResultColorImage1);
//	cvSetImageROI(cvResultColorImage1, cvRect(imageWidth, 0, imageWidth*2-1, imageHeight-1));
//	cvResize(cvRightColorImage, cvResultColorImage1);
//	cvResetImageROI(cvResultColorImage1);
//	cvCopy(cvResultColorImage1,cvResultColorImage2);
//	for(unsigned int pair=0; pair<matchResult.size()/2;pair++){
//		lineIDLeft = matchResult[2*pair];
//		lineIDRight= matchResult[2*pair+1];
//		startPoint = cvPoint(int(linesInLeft[lineIDLeft][0].startPointX),int(linesInLeft[lineIDLeft][0].startPointY));
//		endPoint   = cvPoint(int(linesInRight[lineIDRight][0].startPointX+imageWidth),int(linesInRight[lineIDRight][0].startPointY));
//		cvLine( cvResultColorImage2,startPoint,endPoint,CV_RGB(r1[pair],g1[pair],b1[pair]),2, CV_AA);
//	}
//	cvAddWeighted( cvResultColorImage1, 0.5, cvResultColorImage2, 0.5, 0.0, cvResultColorImage);
//
//	cvSaveImage("LBDSG.png",cvResultColorImage);
//	cvReleaseImage(&cvResultColorImage);
//	cvReleaseImage(&cvResultColorImage1);
//	cvReleaseImage(&cvResultColorImage2);
//	cvReleaseImage(&cvLeftImage);
//	cvReleaseImage(&cvRightImage);
//	cvReleaseImage(&cvLeftColorImage);
//	cvReleaseImage(&cvRightColorImage);
//	cout<<"number of total matches = "<<matchResult.size()/2<<endl;
  ///////////####################################################################
  ///////////####################################################################
}
