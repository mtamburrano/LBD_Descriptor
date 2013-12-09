/*
 * LineStructure.hh
 *
 *  Created on: Dec 5, 2011
 *      Author: lz
 */

#ifndef LINESTRUCTURE_HH_
#define LINESTRUCTURE_HH_

#include <vector>
// A 2D line (normal equation parameters).
struct SingleLine
{
	//note: rho and theta are based on coordinate origin, i.e. the top-left corner of image
	double rho;//unit: pixel length
	double theta;//unit: rad
	double linePointX;// = rho * cos(theta);
	double linePointY;// = rho * sin(theta);
	//for EndPoints, the coordinate origin is the top-left corner of image.
	double startPointX;
	double startPointY;
	double endPointX;
	double endPointY;
	//direction of a line, the angle between positive line direction (dark side is in the left) and positive X axis.
	double direction;
	//mean gradient magnitude
	double gradientMagnitude;
	//mean gray value of pixels in dark side of line
	double darkSideGrayValue;
	//mean gray value of pixels in light side of line
	double lightSideGrayValue;
	//the length of line
	double lineLength;
	//the width of line;
	double width;
	//number of pixels
	int numOfPixels;
	//the decriptor of line
	std::vector<double> descriptor;
};

// Specifies a vector of lines.
typedef std::vector<SingleLine> Lines_list;

struct OctaveSingleLine
{
	/*endPoints, the coordinate origin is the top-left corner of the original image.
	 *startPointX = sPointInOctaveX * (factor)^octaveCount;	*/
	float startPointX;
	float startPointY;
	float endPointX;
	float endPointY;
	//endPoints, the coordinate origin is the top-left corner of the octave image.
	float sPointInOctaveX;
	float sPointInOctaveY;
	float ePointInOctaveX;
	float ePointInOctaveY;
	//direction of a line, the angle between positive line direction (dark side is in the left) and positive X axis.
	float direction;
	//the summation of gradient magnitudes of pixels on lines
	float salience;
	//the length of line
	float lineLength;
	//number of pixels
	unsigned int numOfPixels;
	//the octave which this line is detected
	unsigned int octaveCount;
	//the decriptor of line
	std::vector<float> descriptor;
};

// Specifies a vector of lines.
typedef std::vector<OctaveSingleLine> LinesVec;

typedef std::vector<LinesVec> ScaleLines;//each element in ScaleLines is a vector of lines which corresponds the same line detected in different octave images.

#endif /* LINESTRUCTURE_HH_ */
