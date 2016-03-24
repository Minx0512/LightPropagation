/*
 * ImageFilter.cpp
 *
 *  Created on: Feb 4, 2016
 *      Author: matthias
 */

#include "ImageFilter.hpp"
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkLookupTable.h>
#include <vtkImageShiftScale.h>
#include <vtkImageMapToColors.h>
#include <vtkImageLogarithmicScale.h>



ImageFilter::ImageFilter() {
	// TODO Auto-generated constructor stub

	imageOut = vtkSmartPointer<vtkImageData>::New();
	imageIn = vtkSmartPointer<vtkImageData>::New();
	imLogS = vtkSmartPointer<vtkImageLogarithmicScale>::New();
	shiftScaleFilter = vtkSmartPointer<vtkImageShiftScale>::New();
	lookupTable = vtkSmartPointer<vtkLookupTable>::New();
	colorimageS = vtkSmartPointer<vtkImageMapToColors>::New();
	colors = vtkSmartPointer<vtkColorTransferFunction>::New();
	constant = 20;
	wToBar = 10;
	barWidth = 10;
	oldRange = 0;
	range = 255;
	colormap = 0;
	newRange = 255;
	scale = 1;
	max = 255;
	shift = 0;
	rangeMin = 0;
	rangeMax = 255;
	scaleMode = IF_LOG10;

	coordSystem = IF_COORDS_RADIAL;

	colorToBar[0] = 0;
	colorToBar[1] = 0;
	colorToBar[2] = 0;



}

ImageFilter::~ImageFilter() {
	// TODO Auto-generated destructor stub
}

void ImageFilter::SetInput(vtkSmartPointer<vtkImageData> imageInput){


	imageIn->SetExtent(imageInput->GetExtent());
	imageIn->SetDimensions(imageInput->GetDimensions());
	imageIn = imageInput;


}


void ImageFilter::SetColorBarWidth(int width){
	barWidth = width;
}

void ImageFilter::SetBlackWidth(int width){
	wToBar = width;

}

void ImageFilter::SetColorToBar(double  r, double g, double b){

	colorToBar[0] = r;
	colorToBar[1] = g;
	colorToBar[2] = b;

}


void ImageFilter::SetColorMap(int mapID){
	colormap = mapID;

}

void ImageFilter::SetCoordinatesystem(int mode){
	coordSystem = mode;
}

void ImageFilter::SetCoordinatesystemToRadial(){
	coordSystem = IF_COORDS_RADIAL;
}
void ImageFilter::SetCoordinatesystemToCartesian(){
	coordSystem = IF_COORDS_CARTESIAN;
}

void ImageFilter::SetColorMapToCool2Warm(){
	colormap = COLORMAP_C2W;
}

void ImageFilter::SetColorMapToCool2WarmLab(){
	colormap = COLORMAP_C2WL;
}
void ImageFilter::SetColorMapToRainbow(){
	colormap = COLORMAP_RAINBOW;

}
void ImageFilter::SetColorMapToGrayScale(){
	colormap = COLORMAP_GRAYSCALE;

}
void ImageFilter::SetColorMapToBlackbody(){
	colormap = COLORMAP_BLACKBODY;
}

void ImageFilter::SetScaleString(std::string scaleString){

	EvalScaleString(scaleString);
}
void ImageFilter::SetCoordsystemString(std::string coordsysString){
	EvalCoordsystemString(coordsysString);
}
void ImageFilter::SetColormapString(std::string colormapString){
	EvalColormapString(colormapString);
}


void ImageFilter::SetToLog10Scale(){

	scaleMode = IF_LOG10;

}

void ImageFilter::SetToLinearScale(){
	scaleMode = IF_LIN;
}

void ImageFilter::Log10Scale(vtkSmartPointer<vtkImageData> imageInput){

	switch(coordSystem){
		case IF_COORDS_CARTESIAN:
			Log10ScaleCartesian(imageInput);
			break;
		case IF_COORDS_RADIAL:
			Log10ScaleRadial(imageInput);
			break;

		default: Log10ScaleCartesian(imageInput);
			break;
	}



}


void ImageFilter::RadialExtent(vtkSmartPointer<vtkImageData> imageInput){

	RadialExtent(imageInput,imageOut);

}

void ImageFilter::RadialExtent(vtkSmartPointer<vtkImageData> imageIn, vtkSmartPointer<vtkImageData> imageOutput){
	imageIn->GetExtent(imExtend);
	imageIn->GetDimensions(dims);

	imageOutput->SetExtent(int(-dims[0]),int(dims[0])+wToBar+barWidth,int(-dims[0]),int(dims[0]),0,0);
	//	imageOut->SetOrigin(0,0,0);
	imageOutput->AllocateScalars(VTK_UNSIGNED_CHAR,3);

	imageOutput->GetExtent(imExtend);

}
void ImageFilter::CartesianExtent(vtkSmartPointer<vtkImageData> imageInput){

	CartesianExtent(imageInput,imageOut);

}
void ImageFilter::CartesianExtent(vtkSmartPointer<vtkImageData> imageIn, vtkSmartPointer<vtkImageData> imageOutput){
		imageIn->GetExtent(imExtend);
		//imageInput->GetDimensions(dims);

		//imageOutput->SetDimensions(dims[0]+wToBar+barWidth,dims[1],1);

	//	imageOutput->SetExtent(int(-dims[0]/2),int(dims[0]/2)+wToBar+barWidth,int(-dims[0]/2),int(dims[0]/2),0,1);
		imageOutput->SetExtent(imExtend[0],imExtend[1]+wToBar+barWidth,imExtend[2],imExtend[3],0,0);

	//	imageOutput->SetOrigin(0,0,0);

		imageOutput->AllocateScalars(VTK_UNSIGNED_CHAR,3);

		imageOutput->GetExtent(imExtend);

}

void ImageFilter::EvalScaleString(std::string scaleString){


	if(scaleString=="Log10"){
		this->SetToLog10Scale();
	}else if(scaleString=="Linear"){
		this->SetToLinearScale();
	}else{
		this->SetToLog10Scale();
	}

}
void ImageFilter::EvalCoordsystemString(std::string coordsysString){

	if(coordsysString=="cartesian"){
		this->SetCoordinatesystemToCartesian();
	}else if(coordsysString=="radial"){
		this->SetCoordinatesystemToRadial();
	}else{
		this->SetCoordinatesystemToRadial();
	}

}

void ImageFilter::EvalColormapString(std::string colormapString){

	if(colormapString=="c2w"){
		this->SetColorMapToCool2Warm();
	}else if(colormapString=="rainbow"){
		this->SetColorMapToRainbow();
	}else if(colormapString=="c2wl"){
		this->SetColorMapToCool2WarmLab();
	}else if(colormapString=="grayscale"){
		this->SetColorMapToGrayScale();
	}else if(colormapString=="blackbody"){
		this->SetColorMapToBlackbody();
	}else{
		this->SetColorMapToCool2Warm();
	}

}


void ImageFilter::SetRGBColor(int x,int y,double color[3],vtkSmartPointer<vtkImageData> image){

	image->SetScalarComponentFromDouble(x,y,0,0,color[0]*255);
	image->SetScalarComponentFromDouble(x,y,0,1,color[1]*255);
	image->SetScalarComponentFromDouble(x,y,0,2,color[2]*255);

}
void ImageFilter::InitLogImage(vtkSmartPointer<vtkImageData> imageInput){
	InitLogImage(imageInput,imageOut);


}
void ImageFilter::InitLogImage(vtkSmartPointer<vtkImageData> imageInput,vtkSmartPointer<vtkImageData> imageOutput){

		double dimY = imExtend[3] - imExtend[2];
		int a = imExtend[1]-wToBar-barWidth;
		int b = imExtend[1]-barWidth;

		shift = -1.0f * imageInput->GetScalarRange()[0];

		shift = (shift==0) ? 1 : shift;

		max = imageInput->GetScalarRange()[1]+shift;

		for(int y=imExtend[2];y<=imExtend[3];y++){

			for(int x=imExtend[0];x<=imExtend[1];x++){

				 if(x > a && x <= b){
					 // space to colorbar
					 SetRGBColor(x,y,colorToBar,imageOutput);


				 }else if(x > b){
					// colorbar
					double v = log10(y+abs(imExtend[2])+shift);
					double *colorV = colors->GetColor((double)v/(log10(dimY)));

					SetRGBColor(x,y,colorV,imageOutput);

				}else if(x <= a){
					double *colorZ = colors->GetColor(0.0);
					SetRGBColor(x,y,colorZ,imageOutput);


				}

			}
		}



}

void ImageFilter::InitLinearImage(vtkSmartPointer<vtkImageData> imageInput){
	InitLinearImage(imageInput,imageOut);


}
void ImageFilter::InitLinearImage(vtkSmartPointer<vtkImageData> imageInput,vtkSmartPointer<vtkImageData> imageOutput){

		double dimY = imExtend[3] - imExtend[2];
		int a = imExtend[1]-wToBar-barWidth;
		int b = imExtend[1]-barWidth;

		for(int y=imExtend[2];y<=imExtend[3];y++){

			for(int x=imExtend[0];x<=imExtend[1];x++){

				 if(x > a && x <= b){
					 // space to colorbar
					 SetRGBColor(x,y,colorToBar,imageOutput);

				 }else if(x > b){
					// colorbar
					double v = y-imExtend[2];
					double *colorV = colors->GetColor((double)v/dimY);

					SetRGBColor(x,y,colorV,imageOutput);

				}else if(x <= a){
					double *colorZ = colors->GetColor(0.0);
					SetRGBColor(x,y,colorZ,imageOutput);

				}

			}
		}



}

void ImageFilter::Log10ScaleRadial(vtkSmartPointer<vtkImageData> imageInput){

	RadialExtent(imageInput);

	InitLogImage(imageInput);

	max = imageInput->GetScalarRange()[1];

	for(int phiIdx=0;phiIdx<dims[1];phiIdx++){

		for(int r=0;r<dims[0]/2;r++){

			double v = log10(imageInput->GetScalarComponentAsDouble(r,phiIdx,0,0) + shift);
			double *color2 = colors->GetColor(v/log10(max));
			int cosPhi = int(round(r*cos(phiIdx*360/dims[1])));
			int sinPhi = int(round(r*sin(phiIdx*360/dims[1])));

			SetRGBColor(cosPhi,sinPhi,color2,imageOut);

		}
	}

	range = int(ceil(imageOut->GetScalarRange()[1] - imageOut->GetScalarRange()[0]));


}

void ImageFilter::Log10ScaleCartesian(vtkSmartPointer<vtkImageData> imageInput){

	CartesianExtent(imageInput);

	InitLogImage(imageInput);

	for(int y=imExtend[2];y<=imExtend[3];y++){
		for(int x=imExtend[0];x<=imExtend[1]-wToBar-barWidth;x++){

			double v = log10(imageInput->GetScalarComponentAsDouble(x,y,0,0) + shift);
			double *color2 = colors->GetColor(v/log10(max));
			///if(v>0){

				//cout << "V(" << x << ", " << y << "): " << v  << " | "  << v/log(max)<< endl;
			//}
			SetRGBColor(x,y,color2,imageOut);

		}
	}

	range = int(ceil(imageOut->GetScalarRange()[1] - imageOut->GetScalarRange()[0]));


}


void ImageFilter::LinearScale(vtkSmartPointer<vtkImageData> imageInput){

	switch(coordSystem){
			case IF_COORDS_CARTESIAN:
				LinearScaleCartesian(imageInput);
				break;
			case IF_COORDS_RADIAL:
				LinearScaleRadial(imageInput);
				break;

			default: LinearScaleCartesian(imageInput);
				break;
		}


}

void ImageFilter::LinearScaleRadial(vtkSmartPointer<vtkImageData> imageInput){

	RadialExtent(imageInput);
	InitLinearImage(imageInput);

	max = imageInput->GetScalarRange()[1];

	for(int phiIdx=0;phiIdx<dims[1];phiIdx++){
		for(int r=0;r<dims[0]/2;r++){

			double v = imageInput->GetScalarComponentAsDouble(r,phiIdx,0,0);
			double *color2 = colors->GetColor(v/max);
			int cosPhi = int(round(r*cos(phiIdx*360/dims[1])));
			int sinPhi = int(round(r*sin(phiIdx*360/dims[1])));

			SetRGBColor(cosPhi,sinPhi,color2,imageOut);

		}
	}


}
void ImageFilter::LinearScaleCartesian(vtkSmartPointer<vtkImageData> imageInput){
	CartesianExtent(imageInput);
	InitLinearImage(imageInput);

	max = imageInput->GetScalarRange()[1];

	for(int y=imExtend[2];y<=imExtend[3];y++){
		for(int x=imExtend[0];x<=imExtend[1]-wToBar-barWidth;x++){

			double v = imageInput->GetScalarComponentAsDouble(x,y,0,0);
			double *color2 = colors->GetColor(v/max);

			SetRGBColor(x,y,color2,imageOut);
		}
	}



}




/*
void ImageFilter::LUTHot(int range){

	colors->SetColorSpaceToDiverging();
	colors->AddRGBPoint(0.0, 0.230, 0.299, 0.754);
	colors->AddRGBPoint(1.0, 0.706, 0.016, 0.150);



	lookupTable->SetNumberOfTableValues(256);

	//lookupTable->SetTableRange(0,range);
	//lookupTable->SetHueRange(0.667,0.0);
	//lookupTable->SetNumberOfColors(256);


	for(int i= 0;i<256;i++){

		double *color = colors->GetColor(i/255);
		lookupTable->SetTableValue(i,color[0],color[1],color[2],1);
	}

//
	lookupTable->SetTableValue(0, 0, 0,1,1);  //Blue
	lookupTable->SetTableValue(127, 0, 1,0,1);  //green
	lookupTable->SetTableValue(255,1,0,0,1); // red
//
	lookupTable->Build();



}
*/

void ImageFilter::ColorMap(){


	switch(colormap){

		case COLORMAP_C2W:
			Cool2Warm(colors);
			break;
		case COLORMAP_C2WL:
			Cool2WarmLab(colors);
			break;

		case COLORMAP_RAINBOW:
			Rainbow(colors);
			break;

		case COLORMAP_GRAYSCALE:
			GrayScale(colors);
			break;

		case COLORMAP_BLACKBODY:
			Blackbody(colors);
			break;
		default:
			Cool2Warm(colors);
			break;
	}




}


void ImageFilter::Cool2Warm(vtkSmartPointer<vtkColorTransferFunction> cool2warm){

	cool2warm->SetColorSpaceToDiverging();
	cool2warm->SetScaleToLog10();
	cool2warm->AddRGBPoint(0.0, 0.230, 0.299, 0.754);
	cool2warm->AddRGBPoint(1.0, 0.706, 0.016, 0.150);


}

void ImageFilter::Cool2WarmLab(vtkSmartPointer<vtkColorTransferFunction> cool2warmLab){
	cool2warmLab->SetColorSpaceToLab();
	cool2warmLab->AddRGBPoint(0.0, 0.230, 0.299, 0.754);
	cool2warmLab->AddRGBPoint(0.5, 0.865, 0.865, 0.865);
	cool2warmLab->AddRGBPoint(1.0, 0.706, 0.016, 0.150);
}

void ImageFilter::Rainbow(vtkSmartPointer<vtkColorTransferFunction> rainbow){
	rainbow->SetColorSpaceToHSV();
	rainbow->HSVWrapOff();
	rainbow->AddHSVPoint(0.0, 0.66667, 1.0, 1.0);
	rainbow->AddHSVPoint(0.75, 0.120, 1.0, 1.0);
	rainbow->AddHSVPoint(1.0, 0.0, 1.0, 1.0);
}

void ImageFilter::GrayScale(vtkSmartPointer<vtkColorTransferFunction> grayscale){
	grayscale->SetColorSpaceToRGB();
	grayscale->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	grayscale->AddRGBPoint(1.0, 1.0, 1.0, 1.0);

}

void ImageFilter::Blackbody(vtkSmartPointer<vtkColorTransferFunction> blackbody){
	blackbody->SetColorSpaceToRGB();
	blackbody->AddRGBPoint(0.0, 0.0, 0.0, 0.0);
	blackbody->AddRGBPoint(0.4, 0.9, 0.0, 0.0);
	blackbody->AddRGBPoint(0.8, 0.9, 0.9, 0.0);
	blackbody->AddRGBPoint(1.0, 1.0, 1.0, 1.0);
}

void ImageFilter::ShiftScale(vtkSmartPointer<vtkImageData> imageInput){


	//shiftScaleFilter->SetOutputScalarTypeToUnsignedChar();
	shiftScaleFilter->SetInputData(imageInput);

//	shiftScaleFilter->SetShift(-1.0f * imageInput->GetScalarRange()[0]); // brings the lower bound to 0

	oldRange = imageInput->GetScalarRange()[1] - imageInput->GetScalarRange()[0];
	newRange = 255; // We want the output [1,256]


	shiftScaleFilter->SetScale(newRange/oldRange);
	shiftScaleFilter->Update();


}

void ImageFilter::ColorImage(vtkSmartPointer<vtkImageData> imageInput){

	colorimageS->SetLookupTable(lookupTable);
	colorimageS->PassAlphaToOutputOff();
	colorimageS->SetInputData(imageInput);


	colorimageS->Update();

	imageOut = colorimageS->GetOutput();

}

void ImageFilter::DrawGrid(vtkSmartPointer<vtkImageData> imageOutput){

	int ext[6];
	imageOutput->GetExtent(ext);
	double colorGray[3];
	colorGray[0] = 0;
	colorGray[1] = 0;
	colorGray[2] = 0;


	double dimY = ext[3] - ext[2];
	int a = ext[1]-wToBar-barWidth;
	int b = ext[1]-barWidth;

	for(int x=ext[0];x<=a;x++){
		for(int y=ext[2];y<=ext[3];y++){

			if(abs(x)<10){

				SetRGBColor(x,y,colorGray,imageOutput);

			}

		}
	}

}

void ImageFilter::Update(){

	ColorMap();


	switch(scaleMode){

	case IF_LOG10:
		Log10Scale(imageIn);
		break;

	case IF_LIN:
		LinearScale(imageIn);
		break;
	default:
		Log10Scale(imageIn);
		break;


	}

}



vtkImageData* ImageFilter::GetOutput(){


	return imageOut;



}

int ImageFilter::GetCoordinatesystem(){

	return coordSystem;

}
int ImageFilter::GetScaleMode(){

	return scaleMode;

}
int ImageFilter::GetColormapMode(){

	return colormap;

}


void ImageFilter::PrintLogScale(){

	cout << "################### LogScale ###################################" << endl;
	cout << "Image Extent (): " << imExtend[0] << ", " << imExtend[1];
	cout << ", " << imExtend[2] << ", " << imExtend[3] << ", ";
	cout << imExtend[4] << ", " << imExtend[5] << endl;
	cout << "LogScale Scalar Range old (min, max): (" << imageIn->GetScalarRange()[0] << ", ";
	cout <<  imageIn->GetScalarRange()[1] << ")" << endl;
	cout << "Shift: " << shift << endl;
	cout << "Scale: " << scale << endl;
	cout << "Max: " << max << endl;
	cout << "LogScale Scalar Range (min, max): (" << imageOut->GetScalarRange()[0] << ", ";
	cout <<  imageOut->GetScalarRange()[1] << ")" << endl;
	//cout << "Range: " << range << endl;

}

void ImageFilter::PrintShiftScale(){

	cout << "################### ShiftScale ###################################" << endl;
	cout << "Sclar Range (min, max): (" << shiftScaleFilter->GetOutput()->GetScalarRange()[0] << ", ";
	cout <<  shiftScaleFilter->GetOutput()->GetScalarRange()[1] << ")" << endl;
	cout << "Old Range: " << oldRange << endl;
	cout << "Scale: " << shiftScaleFilter->GetScale() << endl;
	cout << "######################################################" << endl;
	imageOut = shiftScaleFilter->GetOutput();

	cout << "Sclar Range (min, max): (" << imageOut->GetScalarRange()[0] << ", ";
	cout <<  imageOut->GetScalarRange()[1] << ")" << endl;

}

