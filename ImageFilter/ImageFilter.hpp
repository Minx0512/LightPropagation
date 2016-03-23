/*
 * ImageFilter.hpp
 *
 *  Created on: Feb 4, 2016
 *      Author: matthias
 */

#ifndef IMAGEFILTER_HPP_
#define IMAGEFILTER_HPP_

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkImageShiftScale.h>
#include <vtkImageMapToColors.h>
#include <vtkImageLogarithmicScale.h>
#include <vtkColorTransferFunction.h>


#define IF_LOG10 0
#define IF_LIN 1
#define IF_COORDS_CARTESIAN 0
#define IF_COORDS_RADIAL 1
#define COLORMAP_C2W 0
#define COLORMAP_C2WL 1
#define COLORMAP_RAINBOW 2
#define COLORMAP_GRAYSCALE 3
#define COLORMAP_BLACKBODY 4


class ImageFilter {

	vtkSmartPointer<vtkImageData> imageIn;
	vtkSmartPointer<vtkImageData> imageOut;
	vtkSmartPointer<vtkImageMapToColors> colorimageS;
	vtkSmartPointer<vtkLookupTable> lookupTable;
	vtkSmartPointer<vtkImageShiftScale> shiftScaleFilter;
	vtkSmartPointer<vtkImageLogarithmicScale> imLogS;
	vtkSmartPointer<vtkColorTransferFunction> colors;



	float constant;
	float rangeMin, rangeMax;
	int wToBar, barWidth;
	int imExtend[6];
	int dims[3];
	double shift,max,newRange, scale;
	int colormap;
	int range, oldRange;
	int scaleMode;
	int coordSystem;
	double colorToBar[3];


public:

	void SetInput(vtkSmartPointer<vtkImageData> imageIn);
	void SetConstant(float c);
	void SetNewRange(float newRange);
	void SetToLog10Scale();
	void SetToLinearScale();

	void RadialExtent(vtkSmartPointer<vtkImageData> imageInput);
	void RadialExtent(vtkSmartPointer<vtkImageData> imageInput, vtkSmartPointer<vtkImageData> imageOutput);

	void CartesianExtent(vtkSmartPointer<vtkImageData> imageInput);
	void CartesianExtent(vtkSmartPointer<vtkImageData> imageInput, vtkSmartPointer<vtkImageData> imageOutput);

	void SetColorMap(int mapID);
	void SetCoordinatesystem(int mode);
	void SetCoordinatesystemToRadial();
	void SetCoordinatesystemToCartesian();

	void SetColorToBar(double  r, double g, double b);


	void SetColorMapToCool2Warm();
	void SetColorMapToCool2WarmLab();
	void SetColorMapToRainbow();
	void SetColorMapToGrayScale();
	void SetColorMapToBlackbody();


	void SetScaleString(std::string scaleString);
	void SetCoordsystemString(std::string coordsysString);
	void SetColormapString(std::string colormapString);


	void SetColorBarWidth(int width);
	void SetBlackWidth(int width);

	vtkImageData* GetOutput();
	int GetCoordinatesystem();
	int GetScaleMode();
	int GetColormapMode();



	void Update();
	void ColorMap();


	void PrintLogScale();
	void PrintShiftScale();




	ImageFilter();
	virtual ~ImageFilter();




private:

	void SetRGBColor(int x,int y,double color[3],vtkSmartPointer<vtkImageData> image);

	void InitLogImage(vtkSmartPointer<vtkImageData> imageInput,vtkSmartPointer<vtkImageData> imageOutput);
	void InitLogImage(vtkSmartPointer<vtkImageData> imageInput);
	void InitLinearImage(vtkSmartPointer<vtkImageData> imageInput,vtkSmartPointer<vtkImageData> imageOutput);
	void InitLinearImage(vtkSmartPointer<vtkImageData> imageInput);

	void EvalScaleString(std::string scaleString);
	void EvalCoordsystemString(std::string coordsysString);
	void EvalColormapString(std::string colormapString);


	void LUTHot(int range);
	void Cool2Warm(vtkSmartPointer<vtkColorTransferFunction> cool2warm);
	void Cool2WarmLab(vtkSmartPointer<vtkColorTransferFunction> cool2warmLab);
	void Rainbow(vtkSmartPointer<vtkColorTransferFunction> rainbow);
	void GrayScale(vtkSmartPointer<vtkColorTransferFunction> grayscale);
	void Blackbody(vtkSmartPointer<vtkColorTransferFunction> blackbody);


	void Log10Scale(vtkSmartPointer<vtkImageData> imageInput);

	void Log10ScaleRadial(vtkSmartPointer<vtkImageData> imageInput);
	void Log10ScaleCartesian(vtkSmartPointer<vtkImageData> imageInput);

	void LinearScale(vtkSmartPointer<vtkImageData> imageInput);

	void LinearScaleRadial(vtkSmartPointer<vtkImageData> imageInput);
	void LinearScaleCartesian(vtkSmartPointer<vtkImageData> imageInput);


	void ShiftScale(vtkSmartPointer<vtkImageData> imageInput);
	void ColorImage(vtkSmartPointer<vtkImageData> imageInput);

	void DrawGrid(vtkSmartPointer<vtkImageData> imageInput);




};

#endif /* IMAGEFILTER_HPP_ */
