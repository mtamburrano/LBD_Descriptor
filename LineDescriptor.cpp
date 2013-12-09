/*
 * LineDescriptor.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: lz
 */

#include "LineDescriptor.hh"
#include <cv.h>
#include <highgui.h>
#include "opencv2/features2d/features2d.hpp"
#include <Base/Debug/TimeMeasure.hh>
#include <time.h>
#define SalienceScale 0.9//0.9

//#define DEBUGLinesInOctaveImages

using namespace std;
using namespace BIAS;

LineDescriptor::LineDescriptor()
{
	srand ( time(NULL) );
//	cout<<"Call LineDescriptor constructor function"<<endl;
	ksize_       = 5;
	numOfOctave_ = 5;//5
	edLineVec_.resize(numOfOctave_);
	for(unsigned int i=0; i<numOfOctave_; i++){
		edLineVec_[i] = new EDLineDetector;
	}
	numOfBand_   = 9;//9 is a good value.
	widthOfBand_ = 7;//widthOfBand_%3 must equal to 0; 7 is a good value.
	gaussCoefL_.resize(widthOfBand_*3);
	double u = (widthOfBand_*3-1)/2;
	double sigma = (widthOfBand_*2+1)/2;// (widthOfBand_*2+1)/2;
	double invsigma2 = -1/(2*sigma*sigma);
	double dis;
	for(int i=0; i<widthOfBand_*3; i++){
		dis = i-u;
		gaussCoefL_[i] = exp(dis*dis*invsigma2);
//		cout<<"gaussCoefL_="<<gaussCoefL_[i]<<endl;
	}
	gaussCoefG_.resize(numOfBand_*widthOfBand_);
	u = (numOfBand_*widthOfBand_-1)/2;
	sigma = u;
	invsigma2 = -1/(2*sigma*sigma);
	for(int i=0; i<numOfBand_*widthOfBand_; i++){
		dis = i-u;
		gaussCoefG_[i] = exp(dis*dis*invsigma2);
//		cout<<"gaussCoefG_="<<gaussCoefG_[i]<<endl;
	}
//	cout<<"LineDescriptor object is constructed"<<endl;
	LowestThreshold = 0.3;//2 is used to show recall ratio;  0.2 is used to show scale space results, 0.35 is used when verify geometric constraints.
	NNDRThreshold   = 0.6;
}

LineDescriptor::LineDescriptor(unsigned int numOfBand, unsigned int widthOfBand){
	srand ( time(NULL) );
//	cout<<"Call LineDescriptor constructor function"<<endl;
	ksize_       = 5;
	numOfOctave_ = 5;//5
	edLineVec_.resize(numOfOctave_);
	for(unsigned int i=0; i<numOfOctave_; i++){
		edLineVec_[i] = new EDLineDetector;
	}
	numOfBand_   = numOfBand;
	widthOfBand_ = widthOfBand;
	gaussCoefL_.resize(widthOfBand_*3);
	double u = (widthOfBand_*3-1)/2;
	double sigma = (widthOfBand_*2+1)/2;// (widthOfBand_*2+1)/2;
	double invsigma2 = -1/(2*sigma*sigma);
	double dis;
	for(int i=0; i<widthOfBand_*3; i++){
		dis = i-u;
		gaussCoefL_[i] = exp(dis*dis*invsigma2);
//		cout<<"gaussCoefL_="<<gaussCoefL_[i]<<endl;
	}
	gaussCoefG_.resize(numOfBand_*widthOfBand_);
	u = (numOfBand_*widthOfBand_-1)/2;
	sigma = u;
	invsigma2 = -1/(2*sigma*sigma);
	for(int i=0; i<numOfBand_*widthOfBand_; i++){
		dis = i-u;
		gaussCoefG_[i] = exp(dis*dis*invsigma2);
//		cout<<"gaussCoefG_="<<gaussCoefG_[i]<<endl;
	}
//	cout<<"LineDescriptor object is constructed"<<endl;
	LowestThreshold = 0.35;//0.35;
	NNDRThreshold   = 0.2;//0.6
}

LineDescriptor::~LineDescriptor(){
	for(unsigned int i=0; i<numOfOctave_; i++){
		if(edLineVec_[i] !=NULL){
			delete edLineVec_[i];
		}
	}
}

/*Line detection method: element in keyLines[i] includes a set of lines which is the same line
 * detected in different octave images.
 */
int LineDescriptor::OctaveKeyLines(BIAS::Image<unsigned char> & image, ScaleLines &keyLines)
{
	unsigned int numOfFinalLine = 0;
#ifdef DEBUGLinesInOctaveImages
		BIAS::TimeMeasure timer;
#endif
	BIAS::Image<float> img;
	BIAS::ImageConvert::ConvertST(image,img, BIAS::ImageBase::ST_float);//copy the input image
	float preSigma2 = 0;//orignal image is not blurred, has zero sigma;
	float curSigma2 = 1.0;//[sqrt(2)]^0=1;
	float factor = sqrt(2);//the down sample factor between connective two octave images
	BIAS::Gauss<float, float> theGauss;
	//BIAS::Rescale<float, float> theScale;
	//theScale.SetBorderHandling(BIAS::FilterBase<float,float>::TBH_full);
	for(unsigned int octaveCount = 0; octaveCount<numOfOctave_; octaveCount++){
#ifdef DEBUGLinesInOctaveImages
		timer.Start();
#endif
		BIAS::Image<float> blur;
		/* Form each level by adding incremental blur from previous level.
		 * curSigma = [sqrt(2)]^octaveCount;
		 * increaseSigma^2 = curSigma^2 - preSigma^2 */
		float increaseSigma = sqrt(curSigma2-preSigma2);
		theGauss.SetSigma(increaseSigma);
		switch(ksize_){
		case 3: theGauss.Filter3x3Grey(img, blur); break;
		case 5: theGauss.Filter5x5Grey(img, blur); break;
		case 7: theGauss.Filter7x7Grey(img, blur); break;
		case 9: theGauss.Filter9x9Grey(img, blur); break;
		case 11: theGauss.Filter11x11Grey(img, blur); break;
		default: theGauss.Filter5x5Grey(img, blur); break;
		}

		//detect line for each octave image;
		BIAS::Image<unsigned char> blurUChar;
		BIAS::ImageConvert::ConvertST(blur,blurUChar,BIAS::ImageBase::ST_unsignedchar);
		if(!edLineVec_[octaveCount]->EDline(blurUChar,true)){
			return -1;
		}
		numOfFinalLine += edLineVec_[octaveCount]->lines_.numOfLines;
#ifdef DEBUGLinesInOctaveImages
		timer.Stop();
		cout<<"Octave "<<octaveCount<<", line num = "<<edLineVec_[octaveCount]->lines_.numOfLines <<", Time to detect line "<<endl;
		timer.Print();
		timer.Reset();
		////////////////////////////////////

		//the following code is to show the detected lines for each octave image
		IplImage* cvImg = cvCreateImage(cvSize(blur.GetWidth(),blur.GetHeight()),IPL_DEPTH_32F, 1);
		IplImage* cvColorImg = cvCreateImage(cvSize(blur.GetWidth(),blur.GetHeight()),IPL_DEPTH_32F, 3);
		float *newImageData = (float*) cvImg->imageDataOrigin;
		memcpy(newImageData, blur.GetImageData(), cvImg->imageSize);
		cvCvtColor( cvImg, cvColorImg,  CV_GRAY2RGB );
		CvScalar s;
		s.val[0] = 0; s.val[1] = 250;  s.val[2] = 0;
		CvPoint point;
		CvFont  font;
		cvInitFont(&font,CV_FONT_HERSHEY_DUPLEX ,1.0,1.0,0,1);
		unsigned int *pLineXCors = edLineVec_[octaveCount]->lines_.xCors.data();
		unsigned int *pLineYCors = edLineVec_[octaveCount]->lines_.yCors.data();
		unsigned int *pLineSID   = edLineVec_[octaveCount]->lines_.sId.data();
		unsigned int offsetInLineArray;
		for(unsigned int i=0; i<edLineVec_[octaveCount]->lines_.numOfLines; i++){
			for(offsetInLineArray = pLineSID[i]; offsetInLineArray < pLineSID[i+1]; offsetInLineArray++){
				cvSet2D(cvColorImg,pLineYCors[offsetInLineArray],pLineXCors[offsetInLineArray],s);
			}
			//			char buf[10];
			//			sprintf( buf,   "%d ",  i);
			//			point = cvPoint(round(0.5*edLineVec_[octaveCount].lineEndpoints_[i][0]+0.5*edLineVec_[octaveCount].lineEndpoints_[i][2]),
			//					round(0.5*edLineVec_[octaveCount].lineEndpoints_[i][1]+0.5*edLineVec_[octaveCount].lineEndpoints_[i][3]));
			//			cvPutText(cvColorImg,buf,point,&font,CV_RGB(255,0,0));
		}
		ostringstream imagename;
		imagename<<"cvImage"<<octaveCount<<".png";
		cvSaveImage(imagename.str().c_str(),cvColorImg);
		cvReleaseImage(&cvImg);
		cvReleaseImage(&cvColorImg);
//		if(octaveCount==0){//store the dxImage and dyImage to disk for debugging.
//			BIAS::Matrix<short> dxImg(blur.GetHeight(),blur.GetWidth()), dyImg(blur.GetHeight(),blur.GetWidth());
//			short *pdxImg = edLineVec_[octaveCount]->dxImg_->data.s;
//			short *pdyImg = edLineVec_[octaveCount]->dyImg_->data.s;
//			short *pdxMat = dxImg.GetData();
//			short *pdyMat = dyImg.GetData();
//			unsigned int size = blur.GetHeight()*blur.GetWidth()*sizeof(short);
//			memcpy(pdxMat,pdxImg,size);
//			memcpy(pdyMat,pdyImg,size);
//			dxImg.Save("dxImage.txt");
//			dyImg.Save("dyImage.txt");
//		}
#endif
		////////////////////////////////////
		//down sample the current octave image to get the next octave image
		img.Init((int)((float) blur.GetWidth() / factor), (int)((float) blur.GetHeight() / factor),1);
		sample(blur.GetImageData(),img.GetImageData(), factor, blur.GetWidth(),  blur.GetHeight());
		preSigma2 = curSigma2;
		curSigma2 = curSigma2*2;
	}
#ifdef DEBUGLinesInOctaveImages
	timer.Start();
#endif
	/*lines which correspond to the same line in the octave images will be stored in the same element of ScaleLines.*/
	std::vector<OctaveLine> octaveLines(numOfFinalLine);//store the lines in OctaveLine structure
	numOfFinalLine = 0;//store the number of finally accepted lines in ScaleLines
  unsigned int lineIDInScaleLineVec = 0;
	float dx, dy;
	for(unsigned int lineCurId=0;lineCurId<edLineVec_[0]->lines_.numOfLines;lineCurId++){//add all line detected in the original image
		octaveLines[numOfFinalLine].octaveCount    = 0;
		octaveLines[numOfFinalLine].lineIDInOctave = lineCurId;
		octaveLines[numOfFinalLine].lineIDInScaleLineVec = lineIDInScaleLineVec;
		dx = fabs(edLineVec_[0]->lineEndpoints_[lineCurId][0]-edLineVec_[0]->lineEndpoints_[lineCurId][2]);//x1-x2
		dy = fabs(edLineVec_[0]->lineEndpoints_[lineCurId][1]-edLineVec_[0]->lineEndpoints_[lineCurId][3]);//y1-y2
		octaveLines[numOfFinalLine].lineLength = sqrt(dx*dx+dy*dy);
		numOfFinalLine++;
		lineIDInScaleLineVec++;
	}

	float *scale = new float[numOfOctave_];
	scale[0] = 1;
	for(unsigned int octaveCount = 1; octaveCount<numOfOctave_; octaveCount++ ){
		scale[octaveCount] = factor * scale[octaveCount-1];
	}


	float rho1, rho2, tempValue;
	float direction, near, length;
	unsigned int octaveID, lineIDInOctave;
	/*more than one octave image, organize lines in scale space.
	 *lines corresponding to the same line in octave images should have the same index in the ScaleLineVec */
	if(numOfOctave_>1){
		float twoPI = 2*M_PI;
		unsigned int closeLineID;
		float endPointDis,minEndPointDis,minLocalDis,maxLocalDis;
		float lp0,lp1, lp2, lp3, np0,np1, np2, np3;
		for(unsigned int octaveCount = 1; octaveCount<numOfOctave_; octaveCount++){
			/*for each line in current octave image, find their corresponding lines in the octaveLines,
			 *give them the same value of lineIDInScaleLineVec*/
			for(unsigned int lineCurId=0;lineCurId<edLineVec_[octaveCount]->lines_.numOfLines;lineCurId++){
				rho1 = scale[octaveCount] *  fabs(edLineVec_[octaveCount]->lineEquations_[lineCurId][2]);
				/*nearThreshold depends on the distance of the image coordinate origin to current line.
				 *so nearThreshold = rho1 * nearThresholdRatio, where nearThresholdRatio = 1-cos(10*pi/180) = 0.0152*/
				tempValue = rho1 * 0.0152;
				float nearThreshold = (tempValue>6)?(tempValue):6;
				nearThreshold = (nearThreshold<12)?nearThreshold:12;
				dx = fabs(edLineVec_[octaveCount]->lineEndpoints_[lineCurId][0]-edLineVec_[octaveCount]->lineEndpoints_[lineCurId][2]);//x1-x2
				dy = fabs(edLineVec_[octaveCount]->lineEndpoints_[lineCurId][1]-edLineVec_[octaveCount]->lineEndpoints_[lineCurId][3]);//y1-y2
				length = scale[octaveCount] * sqrt(dx*dx+dy*dy);
				minEndPointDis = 12;
				for(unsigned int lineNextId=0; lineNextId<numOfFinalLine;lineNextId++){
					octaveID = octaveLines[lineNextId].octaveCount;
					if(octaveID==octaveCount){//lines in the same layer of octave image should not be compared.
						break;
					}
					lineIDInOctave = octaveLines[lineNextId].lineIDInOctave;
					/*first check whether current line and next line are parallel.
					 *If line1:a1*x+b1*y+c1=0 and line2:a2*x+b2*y+c2=0 are parallel, then
					 *-a1/b1=-a2/b2, i.e., a1b2=b1a2.
					 *we define parallel=fabs(a1b2-b1a2)
					 *note that, in EDLine class, we have normalized the line equations to make a1^2+ b1^2 = a2^2+ b2^2 = 1*/
					direction = fabs(edLineVec_[octaveCount]->lineDirection_[lineCurId] -
							edLineVec_[octaveID]->lineDirection_[lineIDInOctave]);
					if(direction>0.1745&&(twoPI - direction>0.1745)){
						continue;//the angle between two lines are larger than 10degrees(i.e. 10*pi/180=0.1745), they are not close to parallel.
					}
					/*now check whether current line and next line are near to each other.
					 *If line1:a1*x+b1*y+c1=0 and line2:a2*x+b2*y+c2=0 are near in image, then
					 *rho1 = |a1*0+b1*0+c1|/sqrt(a1^2+b1^2) and rho2 = |a2*0+b2*0+c2|/sqrt(a2^2+b2^2) should close.
					 *In our case, rho1 = |c1| and rho2 = |c2|, because sqrt(a1^2+b1^2) = sqrt(a2^2+b2^2) = 1;
					 *note that, lines are in different octave images, so we define near =  fabs(scale*rho1 - rho2) or
					 *where scale is the scale factor between to octave images*/
					rho2 = scale[octaveID] * fabs(edLineVec_[octaveID]->lineEquations_[lineIDInOctave][2]);
					near = fabs(rho1 - rho2);
					if(near>nearThreshold){
						continue;//two line are not near in the image
					}
					/*now check the end points distance between two lines, the scale of  distance is in the original image size.
					 * find the minimal and maximal end points distance*/
					lp0 = scale[octaveCount] *edLineVec_[octaveCount]->lineEndpoints_[lineCurId][0];
					lp1 = scale[octaveCount] *edLineVec_[octaveCount]->lineEndpoints_[lineCurId][1];
					lp2 = scale[octaveCount] *edLineVec_[octaveCount]->lineEndpoints_[lineCurId][2];
					lp3 = scale[octaveCount] *edLineVec_[octaveCount]->lineEndpoints_[lineCurId][3];
					np0 = scale[octaveID] * edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][0];
					np1 = scale[octaveID] * edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][1];
					np2 = scale[octaveID] * edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][2];
					np3 = scale[octaveID] * edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][3];
					//L1(0,1)<->L2(0,1)
					dx = lp0 - np0;
					dy = lp1 - np1;
					endPointDis = sqrt(dx*dx + dy*dy);
					minLocalDis = endPointDis;
					maxLocalDis = endPointDis;
					//L1(2,3)<->L2(2,3)
					dx = lp2 - np2;
					dy = lp3 - np3;
					endPointDis = sqrt(dx*dx + dy*dy);
					minLocalDis = (endPointDis<minLocalDis)?endPointDis:minLocalDis;
					maxLocalDis = (endPointDis>maxLocalDis)?endPointDis:maxLocalDis;
					//L1(0,1)<->L2(2,3)
					dx = lp0 - np2;
					dy = lp1 - np3;
					endPointDis = sqrt(dx*dx + dy*dy);
					minLocalDis = (endPointDis<minLocalDis)?endPointDis:minLocalDis;
					maxLocalDis = (endPointDis>maxLocalDis)?endPointDis:maxLocalDis;
					//L1(2,3)<->L2(0,1)
					dx = lp2 - np0;
					dy = lp3 - np1;
					endPointDis = sqrt(dx*dx + dy*dy);
					minLocalDis = (endPointDis<minLocalDis)?endPointDis:minLocalDis;
					maxLocalDis = (endPointDis>maxLocalDis)?endPointDis:maxLocalDis;

					if((maxLocalDis<0.8*(length+octaveLines[lineNextId].lineLength))&&(minLocalDis<minEndPointDis)){//keep the closest line
						minEndPointDis = minLocalDis;
						closeLineID = lineNextId;
					}
				}
				//add current line into octaveLines
				if(minEndPointDis<12){
					octaveLines[numOfFinalLine].lineIDInScaleLineVec = octaveLines[closeLineID].lineIDInScaleLineVec;
				}else{
					octaveLines[numOfFinalLine].lineIDInScaleLineVec = lineIDInScaleLineVec;
					lineIDInScaleLineVec++;
				}
				octaveLines[numOfFinalLine].octaveCount    = octaveCount;
				octaveLines[numOfFinalLine].lineIDInOctave = lineCurId;
				octaveLines[numOfFinalLine].lineLength     = length;
				numOfFinalLine++;
			}
		}//end for(unsigned int octaveCount = 1; octaveCount<numOfOctave_; octaveCount++)
	}//end if(numOfOctave_>1)

	////////////////////////////////////
	//Reorganize the detected lines into keyLines
	keyLines.clear();
	keyLines.resize(lineIDInScaleLineVec);
  unsigned int tempID;
	float s1,e1,s2,e2;
	bool shouldChange;
	OctaveSingleLine singleLine;
	for(unsigned int  lineID = 0;lineID < numOfFinalLine; lineID++){
		lineIDInOctave = octaveLines[lineID].lineIDInOctave;
		octaveID       = octaveLines[lineID].octaveCount;
		direction      = edLineVec_[octaveID]->lineDirection_[lineIDInOctave];
		singleLine.octaveCount = octaveID;
		singleLine.direction = direction;
		singleLine.lineLength = octaveLines[lineID].lineLength;
		singleLine.salience  = edLineVec_[octaveID]->lineSalience_[lineIDInOctave];
		singleLine.numOfPixels = edLineVec_[octaveID]->lines_.sId[lineIDInOctave+1]-
		                         edLineVec_[octaveID]->lines_.sId[lineIDInOctave];
		//decide the start point and end point
		shouldChange = false;
		s1 = edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][0];//sx
		s2 = edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][1];//sy
		e1 = edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][2];//ex
		e2 = edLineVec_[octaveID]->lineEndpoints_[lineIDInOctave][3];//ey
		dx = e1 - s1;//ex-sx
		dy = e2 - s2;//ey-sy
		if(direction>=-0.75*M_PI&&direction<-0.25*M_PI){
			if(dy>0){shouldChange = true;}
		}
		if(direction>=-0.25*M_PI&&direction<0.25*M_PI){
			if(dx<0){shouldChange = true;}
		}
		if(direction>=0.25*M_PI&&direction<0.75*M_PI){
			if(dy<0){shouldChange = true;}
		}
		if((direction>=0.75*M_PI&&direction<M_PI)||(direction>=-M_PI&&direction<-0.75*M_PI)){
			if(dx>0){shouldChange = true;}
		}
		tempValue = scale[octaveID];
		if(shouldChange){
			singleLine.sPointInOctaveX = e1;
			singleLine.sPointInOctaveY = e2;
			singleLine.ePointInOctaveX = s1;
			singleLine.ePointInOctaveY = s2;
			singleLine.startPointX = tempValue * e1;
			singleLine.startPointY = tempValue * e2;
			singleLine.endPointX   = tempValue * s1;
			singleLine.endPointY   = tempValue * s2;
		}else{
			singleLine.sPointInOctaveX = s1;
			singleLine.sPointInOctaveY = s2;
			singleLine.ePointInOctaveX = e1;
			singleLine.ePointInOctaveY = e2;
			singleLine.startPointX = tempValue * s1;
			singleLine.startPointY = tempValue * s2;
			singleLine.endPointX   = tempValue * e1;
			singleLine.endPointY   = tempValue * e2;
		}
		tempID = octaveLines[lineID].lineIDInScaleLineVec;
		keyLines[tempID].push_back(singleLine);
	}

	////////////////////////////////////
#ifdef DEBUGLinesInOctaveImages
	//the following code is to show lines in ScaleLines
	timer.Stop();
	cout<<"number of final lines = "<<numOfFinalLine<<", number of lines in ScaleLines vector="<<lineIDInScaleLineVec
	<<", Time to find extreme lines is "<<endl;
	timer.Print();
	timer.Reset();
	//generate the color;
	short *r = new short[lineIDInScaleLineVec];
	short *g = new short[lineIDInScaleLineVec];
	short *b = new short[lineIDInScaleLineVec];
	int lowest=10, highest=255;
	int range=(highest-lowest)+1;
	for(int i=0; i <lineIDInScaleLineVec; i++){
		r[i] = lowest+int(rand()%range);
		g[i] = lowest+int(rand()%range);
		b[i] = lowest+int(rand()%range);
	}

	unsigned int offsetInLineArray;
	unsigned int lineID = 0;
	CvPoint startPoint;
	CvPoint endPoint;
	CvFont  font;
	cvInitFont(&font,CV_FONT_HERSHEY_DUPLEX,2.0,1.6,0.2,2,8);
	CvMat* cvImage = cvCreateMat(image.GetHeight(), image.GetWidth(),CV_8UC1);
	CvMat* cvColorImage = cvCreateMat(image.GetHeight(), image.GetWidth(),CV_8UC3);
	memcpy(cvImage->data.ptr, image.GetImageData(), image.GetHeight()*image.GetWidth()*sizeof(unsigned char));
	cvCvtColor( cvImage, cvColorImage,  CV_GRAY2RGB );
	for(unsigned int octaveCount = 0; octaveCount<numOfOctave_; octaveCount++){
		unsigned int *pLineXCors = edLineVec_[octaveCount]->lines_.xCors.data();
		unsigned int *pLineYCors = edLineVec_[octaveCount]->lines_.yCors.data();
		unsigned int *pLineSID   = edLineVec_[octaveCount]->lines_.sId.data();
		ostringstream imagename;
		imagename<<"cvImage"<<octaveCount<<".png";
		IplImage* cvFinalImg = cvLoadImage(imagename.str().c_str());
		CvScalar s;
		tempValue = scale[octaveCount];
		while(lineID<numOfFinalLine&&octaveLines[lineID].octaveCount==octaveCount){
			tempID = octaveLines[lineID].lineIDInScaleLineVec;
			s.val[0] = b[tempID]; s.val[1] = g[tempID];  s.val[2] = r[tempID];
			lineIDInOctave = octaveLines[lineID].lineIDInOctave;
			offsetInLineArray = pLineSID[lineIDInOctave];
			startPoint = cvPoint(pLineXCors[offsetInLineArray],pLineYCors[offsetInLineArray]);
			for(; offsetInLineArray < pLineSID[lineIDInOctave+1]; offsetInLineArray++){
				cvSet2D(cvFinalImg,pLineYCors[offsetInLineArray],pLineXCors[offsetInLineArray],s);
				cvSet2D(cvColorImage,pLineYCors[offsetInLineArray]*tempValue,pLineXCors[offsetInLineArray]*tempValue,s);
			}
			char buf[10];
			sprintf( buf,   "%d ",  octaveLines[lineID].lineIDInOctave);
			cvPutText(cvFinalImg,buf,startPoint,&font,CV_RGB(r[tempID],g[tempID],b[tempID]));
//			startPoint = cvPoint(int(keyLines[lineID].sPointInOctaveX),int(keyLines[lineID].sPointInOctaveY));
//			endPoint   = cvPoint(int(keyLines[lineID].ePointInOctaveX),int(keyLines[lineID].ePointInOctaveY));
//			cvLine( cvFinalImg,startPoint,endPoint,CV_RGB(250,0,0));
//			startPoint = cvPoint(int(keyLines[lineID].startPointX),int(keyLines[lineID].startPointY));
//			endPoint   = cvPoint(int(keyLines[lineID].endPointX),int(keyLines[lineID].endPointY));
//			cvLine( cvColorImage,startPoint,endPoint,CV_RGB(250,0,0));
			lineID++;
		}
		cvSaveImage(imagename.str().c_str(),cvFinalImg);
		cvReleaseImage(&cvFinalImg);
	}
  cvNamedWindow("LineColorImage", CV_WINDOW_AUTOSIZE);
  cvShowImage("LineColorImage", cvColorImage);
  cvSaveImage("cvImageAll.png",cvColorImage);
  cvWaitKey(0);
  delete [] r;
  delete [] g;
  delete [] b;
  cvReleaseMat(&cvImage);
  cvReleaseMat(&cvColorImage);
#endif

  delete [] scale;
  return 1;
}

/*The definitions of line descriptor,mean values of {g_dL>0},{g_dL<0},{g_dO>0},{g_dO<0} of each row in band
 *and std values of sum{g_dL>0},sum{g_dL<0},sum{g_dO>0},sum{g_dO<0} of each row in band.
 * With overlap region. */
int LineDescriptor::ComputeLBD_(ScaleLines &keyLines)
{
	//the default length of the band is the line length.
	short numOfFinalLine = keyLines.size();
	float *dL = new float[2];//line direction cos(dir), sin(dir)
	float *dO = new float[2];//the clockwise orthogonal vector of line direction.
	short heightOfLSP = widthOfBand_*numOfBand_;//the height of line support region;
	short descriptorSize = numOfBand_ * 8;//each band, we compute the m( pgdL, ngdL,  pgdO, ngdO) and std( pgdL, ngdL,  pgdO, ngdO);
	float pgdLRowSum;//the summation of {g_dL |g_dL>0 } for each row of the region;
	float ngdLRowSum;//the summation of {g_dL |g_dL<0 } for each row of the region;
	float pgdL2RowSum;//the summation of {g_dL^2 |g_dL>0 } for each row of the region;
	float ngdL2RowSum;//the summation of {g_dL^2 |g_dL<0 } for each row of the region;
	float pgdORowSum;//the summation of {g_dO |g_dO>0 } for each row of the region;
	float ngdORowSum;//the summation of {g_dO |g_dO<0 } for each row of the region;
	float pgdO2RowSum;//the summation of {g_dO^2 |g_dO>0 } for each row of the region;
	float ngdO2RowSum;//the summation of {g_dO^2 |g_dO<0 } for each row of the region;

	float *pgdLBandSum  = new float[numOfBand_];//the summation of {g_dL |g_dL>0 } for each band of the region;
	float *ngdLBandSum  = new float[numOfBand_];//the summation of {g_dL |g_dL<0 } for each band of the region;
	float *pgdL2BandSum = new float[numOfBand_];//the summation of {g_dL^2 |g_dL>0 } for each band of the region;
	float *ngdL2BandSum = new float[numOfBand_];//the summation of {g_dL^2 |g_dL<0 } for each band of the region;
	float *pgdOBandSum  = new float[numOfBand_];//the summation of {g_dO |g_dO>0 } for each band of the region;
	float *ngdOBandSum  = new float[numOfBand_];//the summation of {g_dO |g_dO<0 } for each band of the region;
	float *pgdO2BandSum = new float[numOfBand_];//the summation of {g_dO^2 |g_dO>0 } for each band of the region;
	float *ngdO2BandSum = new float[numOfBand_];//the summation of {g_dO^2 |g_dO<0 } for each band of the region;

	short numOfBitsBand = numOfBand_*sizeof(float);
	short lengthOfLSP; //the length of line support region, varies with lines
	short halfHeight = (heightOfLSP-1)/2;
	short halfWidth;
	short bandID;
	float coefInGaussion;
	float lineMiddlePointX, lineMiddlePointY;
	float sCorX, sCorY,sCorX0, sCorY0;
	short tempCor, xCor, yCor;//pixel coordinates in image plane
	short dx, dy;
	float gDL;//store the gradient projection of pixels in support region along dL vector
	float gDO;//store the gradient projection of pixels in support region along dO vector
	short imageWidth, imageHeight, realWidth;
	short *pdxImg, *pdyImg;
	float *desVec;

	short sameLineSize;
	short octaveCount;
	OctaveSingleLine *pSingleLine;
	for(short lineIDInScaleVec = 0; lineIDInScaleVec<numOfFinalLine; lineIDInScaleVec++){
		sameLineSize = keyLines[lineIDInScaleVec].size();
		for(short lineIDInSameLine = 0; lineIDInSameLine<sameLineSize; lineIDInSameLine++){
			pSingleLine = &(keyLines[lineIDInScaleVec][lineIDInSameLine]);
			octaveCount = pSingleLine->octaveCount;
			pdxImg = edLineVec_[octaveCount]->dxImg_->data.s;
			pdyImg = edLineVec_[octaveCount]->dyImg_->data.s;
			realWidth = edLineVec_[octaveCount]->imageWidth;
			imageWidth  = realWidth -1;
			imageHeight = edLineVec_[octaveCount]->imageHeight-1;
			//initialization
			memset(pgdLBandSum,  0, numOfBitsBand);
			memset(ngdLBandSum, 0, numOfBitsBand);
			memset(pgdL2BandSum,  0, numOfBitsBand);
			memset(ngdL2BandSum, 0, numOfBitsBand);
			memset(pgdOBandSum,  0, numOfBitsBand);
			memset(ngdOBandSum, 0, numOfBitsBand);
			memset(pgdO2BandSum,  0, numOfBitsBand);
			memset(ngdO2BandSum, 0, numOfBitsBand);
			lengthOfLSP = keyLines[lineIDInScaleVec][lineIDInSameLine].numOfPixels;
			halfWidth   = (lengthOfLSP-1)/2;
			lineMiddlePointX = 0.5 * (pSingleLine->sPointInOctaveX +  pSingleLine->ePointInOctaveX);
			lineMiddlePointY = 0.5 * (pSingleLine->sPointInOctaveY +  pSingleLine->ePointInOctaveY);
			/*1.rotate the local coordinate system to the line direction
			 *2.compute the gradient projection of pixels in line support region*/
			dL[0] = cos(pSingleLine->direction);
			dL[1] = sin(pSingleLine->direction);
			dO[0] = -dL[1];
			dO[1] = dL[0];
			sCorX0= -dL[0]*halfWidth + dL[1]*halfHeight + lineMiddlePointX;//hID =0; wID = 0;
			sCorY0= -dL[1]*halfWidth - dL[0]*halfHeight + lineMiddlePointY;
			//      BIAS::Matrix<float> gDLMat(heightOfLSP,lengthOfLSP);
			for(short hID = 0; hID <heightOfLSP; hID++){
				//initialization
				sCorX = sCorX0;
				sCorY = sCorY0;

				pgdLRowSum = 0;
				ngdLRowSum = 0;
				pgdORowSum = 0;
				ngdORowSum = 0;

				for(short wID = 0; wID <lengthOfLSP; wID++){
					tempCor = round(sCorX);
					xCor = (tempCor<0)?0:(tempCor>imageWidth)?imageWidth:tempCor;
					tempCor = round(sCorY);
					yCor = (tempCor<0)?0:(tempCor>imageHeight)?imageHeight:tempCor;
					/* To achieve rotation invariance, each simple gradient is rotated aligned with
					 * the line direction and clockwise orthogonal direction.*/
					dx = pdxImg[yCor*realWidth+xCor];
					dy = pdyImg[yCor*realWidth+xCor];
					gDL = dx * dL[0] + dy * dL[1];
					gDO = dx * dO[0] + dy * dO[1];
					if(gDL>0){
						pgdLRowSum  += gDL;
					}else{
						ngdLRowSum  -= gDL;
					}
					if(gDO>0){
						pgdORowSum  += gDO;
					}else{
						ngdORowSum  -= gDO;
					}
					sCorX +=dL[0];
					sCorY +=dL[1];
					//					gDLMat[hID][wID] = gDL;
				}
				sCorX0 -=dL[1];
				sCorY0 +=dL[0];
				coefInGaussion = gaussCoefG_[hID];
				pgdLRowSum = coefInGaussion * pgdLRowSum;
				ngdLRowSum = coefInGaussion * ngdLRowSum;
				pgdL2RowSum = pgdLRowSum * pgdLRowSum;
				ngdL2RowSum = ngdLRowSum * ngdLRowSum;
				pgdORowSum = coefInGaussion * pgdORowSum;
				ngdORowSum = coefInGaussion * ngdORowSum;
				pgdO2RowSum = pgdORowSum * pgdORowSum;
				ngdO2RowSum = ngdORowSum * ngdORowSum;
				//compute {g_dL |g_dL>0 }, {g_dL |g_dL<0 },
				//{g_dO |g_dO>0 }, {g_dO |g_dO<0 } of each band in the line support region
				//first, current row belong to current band;
				bandID = hID/widthOfBand_;
				coefInGaussion = gaussCoefL_[hID%widthOfBand_+widthOfBand_];
				pgdLBandSum[bandID] +=  coefInGaussion * pgdLRowSum;
				ngdLBandSum[bandID] +=  coefInGaussion * ngdLRowSum;
				pgdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdL2RowSum;
				ngdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdL2RowSum;
				pgdOBandSum[bandID] +=  coefInGaussion * pgdORowSum;
				ngdOBandSum[bandID] +=  coefInGaussion * ngdORowSum;
				pgdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdO2RowSum;
				ngdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdO2RowSum;
				/* In order to reduce boundary effect along the line gradient direction,
				 * a row's gradient will contribute not only to its current band, but also
				 * to its nearest upper and down band with gaussCoefL_.*/
				bandID--;
				if(bandID>=0){//the band above the current band
					coefInGaussion = gaussCoefL_[hID%widthOfBand_ + 2*widthOfBand_];
					pgdLBandSum[bandID] +=  coefInGaussion * pgdLRowSum;
					ngdLBandSum[bandID] +=  coefInGaussion * ngdLRowSum;
					pgdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdL2RowSum;
					ngdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdL2RowSum;
					pgdOBandSum[bandID] +=  coefInGaussion * pgdORowSum;
					ngdOBandSum[bandID] +=  coefInGaussion * ngdORowSum;
					pgdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdO2RowSum;
					ngdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdO2RowSum;
				}
				bandID = bandID+2;
				if(bandID<numOfBand_){//the band below the current band
					coefInGaussion = gaussCoefL_[hID%widthOfBand_];
					pgdLBandSum[bandID] +=  coefInGaussion * pgdLRowSum;
					ngdLBandSum[bandID] +=  coefInGaussion * ngdLRowSum;
					pgdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdL2RowSum;
					ngdL2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdL2RowSum;
					pgdOBandSum[bandID] +=  coefInGaussion * pgdORowSum;
					ngdOBandSum[bandID] +=  coefInGaussion * ngdORowSum;
					pgdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * pgdO2RowSum;
					ngdO2BandSum[bandID] +=  coefInGaussion * coefInGaussion * ngdO2RowSum;
				}
			}
			//			gDLMat.Save("gDLMat.txt");
			//			return 0;
			//construct line descriptor
			pSingleLine->descriptor.resize(descriptorSize);
			desVec = pSingleLine->descriptor.data();
			short desID;
			/*Note that the first and last bands only have (lengthOfLSP * widthOfBand_ * 2.0) pixels
			 * which are counted. */
			float invN2 = 1.0/(widthOfBand_ * 2.0);
			float invN3 = 1.0/(widthOfBand_ * 3.0);
			float invN, temp;
			for(bandID = 0; bandID<numOfBand_; bandID++){
				if(bandID==0||bandID==numOfBand_-1){	invN = invN2;
				}else{ invN = invN3;}
				desID = bandID * 8;
				temp = pgdLBandSum[bandID] * invN;
				desVec[desID]   = temp;//mean value of pgdL;
				desVec[desID+4] = sqrt(pgdL2BandSum[bandID] * invN - temp*temp);//std value of pgdL;
				temp = ngdLBandSum[bandID] * invN;
				desVec[desID+1] = temp;//mean value of ngdL;
				desVec[desID+5] = sqrt(ngdL2BandSum[bandID] * invN - temp*temp);//std value of ngdL;

				temp = pgdOBandSum[bandID] * invN;
				desVec[desID+2] = temp;//mean value of pgdO;
				desVec[desID+6] = sqrt(pgdO2BandSum[bandID] * invN - temp*temp);//std value of pgdO;
				temp = ngdOBandSum[bandID] * invN;
				desVec[desID+3] = temp;//mean value of ngdO;
				desVec[desID+7] = sqrt(ngdO2BandSum[bandID] * invN - temp*temp);//std value of ngdO;
			}
			//normalize;
			float tempM, tempS;
			tempM = 0;
			tempS = 0;
			desVec = pSingleLine->descriptor.data();
			for(short i=0; i<numOfBand_; i++){
				tempM += (*desVec) * *(desVec++);//desVec[8*i+0] * desVec[8*i+0];
				tempM += (*desVec) * *(desVec++);//desVec[8*i+1] * desVec[8*i+1];
				tempM += (*desVec) * *(desVec++);//desVec[8*i+2] * desVec[8*i+2];
				tempM += (*desVec) * *(desVec++);//desVec[8*i+3] * desVec[8*i+3];
				tempS += (*desVec) * *(desVec++);//desVec[8*i+4] * desVec[8*i+4];
				tempS += (*desVec) * *(desVec++);//desVec[8*i+5] * desVec[8*i+5];
				tempS += (*desVec) * *(desVec++);//desVec[8*i+6] * desVec[8*i+6];
				tempS += (*desVec) * *(desVec++);//desVec[8*i+7] * desVec[8*i+7];
			}
			tempM = 1/sqrt(tempM);
			tempS = 1/sqrt(tempS);
			desVec = pSingleLine->descriptor.data();
			for(short i=0; i<numOfBand_; i++){
				(*desVec) = *(desVec++) * tempM;//desVec[8*i] =  desVec[8*i] * tempM;
				(*desVec) = *(desVec++) * tempM;//desVec[8*i+1] =  desVec[8*i+1] * tempM;
				(*desVec) = *(desVec++) * tempM;//desVec[8*i+2] =  desVec[8*i+2] * tempM;
				(*desVec) = *(desVec++) * tempM;//desVec[8*i+3] =  desVec[8*i+3] * tempM;
				(*desVec) = *(desVec++) * tempS;//desVec[8*i+4] =  desVec[8*i+4] * tempS;
				(*desVec) = *(desVec++) * tempS;//desVec[8*i+5] =  desVec[8*i+5] * tempS;
				(*desVec) = *(desVec++) * tempS;//desVec[8*i+6] =  desVec[8*i+6] * tempS;
				(*desVec) = *(desVec++) * tempS;//desVec[8*i+7] =  desVec[8*i+7] * tempS;
			}
			/*In order to reduce the influence of non-linear illumination,
			 *a threshold is used to limit the value of element in the unit feature
			 *vector no larger than this threshold. In Z.Wang's work, a value of 0.4 is found
			 *empirically to be a proper threshold.*/
			desVec = pSingleLine->descriptor.data();
			for(short i=0; i<descriptorSize; i++ ){
				if(desVec[i]>0.4){
					desVec[i]=0.4;
				}
			}
			//re-normalize desVec;
			temp = 0;
			for(short i=0; i<descriptorSize; i++){
				temp += desVec[i] * desVec[i];
			}
			temp = 1/sqrt(temp);
			for(short i=0; i<descriptorSize; i++){
				desVec[i] =  desVec[i] * temp;
			}
		}//end for(short lineIDInSameLine = 0; lineIDInSameLine<sameLineSize; lineIDInSameLine++)
	}//end for(short lineIDInScaleVec = 0; lineIDInScaleVec<numOfFinalLine; lineIDInScaleVec++)

	delete [] dL;
	delete [] dO;
	delete [] pgdLBandSum;
	delete [] ngdLBandSum;
	delete [] pgdL2BandSum;
	delete [] ngdL2BandSum;
	delete [] pgdOBandSum;
	delete [] ngdOBandSum;
	delete [] pgdO2BandSum;
	delete [] ngdO2BandSum;
}



int LineDescriptor::GetLineDescriptor(BIAS::Image<unsigned char> & image, ScaleLines &keyLines)
{
	if(!OctaveKeyLines(image,keyLines)){
		cout<<"OctaveKeyLines failed"<<endl;
		return -1;
	}

//	BIAS::TimeMeasure timer;
//	timer.Start();
	ComputeLBD_(keyLines);
//	timer.Stop();
//	cout<<"Time to compute line descriptors is "<<endl;
//	timer.Print();
//	timer.Reset();
	return 1;
}



/*Match line by their descriptors.
 *The function will use opencv FlannBasedMatcher to mathc lines. */
int LineDescriptor::MatchLineByDescriptor(ScaleLines &keyLinesLeft, 	ScaleLines &keyLinesRight,
		std::vector<short> &matchLeft, std::vector<short> &matchRight,
		int criteria)
{
	int leftSize = keyLinesLeft.size();
	int rightSize = keyLinesRight.size();
	if(leftSize<1||rightSize<1){
		return -1;
	}

	matchLeft.clear();
	matchRight.clear();

	int desDim = keyLinesLeft[0][0].descriptor.size();
	float *desL, *desR, *desMax, *desOld;
	if(criteria==NearestNeighbor){
		float minDis,dis,temp;
		int corresId;
		for(int idL=0; idL<leftSize; idL++){
			short sameLineSize = keyLinesLeft[idL].size();
			minDis = 100;
			for(short lineIDInSameLines = 0; lineIDInSameLines<sameLineSize; lineIDInSameLines++){
				desOld = keyLinesLeft[idL][lineIDInSameLines].descriptor.data();
				for(int idR=0; idR<rightSize; idR++){
					short sameLineSizeR = keyLinesRight[idR].size();
					for(short lineIDInSameLinesR = 0; lineIDInSameLinesR<sameLineSizeR; lineIDInSameLinesR++){
						desL = desOld;
						desR = keyLinesRight[idR][lineIDInSameLinesR].descriptor.data();
						desMax = desR+desDim;
						dis = 0;
						while(desR<desMax){
							temp = *(desL++) - *(desR++);
							dis += temp*temp;
						}
						dis = sqrt(dis);
						if(dis<minDis){
							minDis = dis;
							corresId = idR;
						}
					}
				}//end for(int idR=0; idR<rightSize; idR++)
			}//end for(short lineIDInSameLines = 0; lineIDInSameLines<sameLineSize; lineIDInSameLines++)
			if(minDis<LowestThreshold){
				matchLeft.push_back(idL);
				matchRight.push_back(corresId);
			}
		}// end for(int idL=0; idL<leftSize; idL++)
	}
}






















