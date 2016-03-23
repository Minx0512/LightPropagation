/*
 * CalculatePropagation.cpp
 *
 *  Created on: Mar 15, 2015
 *      Author: matthias
 */
#define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include "CalculatePropagation.hpp"
#include <CL/cl.hpp>
#include "../PSNR/PSNR.hpp"

//#include "SDKUtil/CLUtil.hpp"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkDoubleArray.h"
#include "vtkCell.h"
#include <vector>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#define SUCCESS 0
#define FAILURE 1
#define EXPECTED_FAILURE 2

CalculatePropagation::CalculatePropagation() {
	// TODO Auto-generated constructor stub

	verbose = false;
	showPercent = true;
	//PolyImage = vtkPolyData::New();

	useCPU = false;
	useGPU = true;
	isContextSet = false;
	numberGenerator = 0;
	seed = 1;
	MBMZ.assign(2,0);
	inextInextp.assign(3,0);

	numPhotonInteractions = 100;
	absNumberPhotons = 1;
	numberPhotons = 1;
	chance = 0.5;
	w_threshold = 0.001;
	pixelPerUnitFactor = 100; //pixel zu cm
	isReadOnly = true;
	status = 0;
	filename = "MCLightPropagation.cl";
	kernelname = "MCLightPropagation";
	sourceStr = "";

	buildOptions = "-x clc++ -cl-no-signed-zeros";
	//buildOptions = "-x clc++";

	//platform = NULL;
	//devices = NULL;

	poolSize = 100;

	iter = 1;
	fetalTissueID = 1;
	DetGridDims.resize(6,1);
	eP0XYZ.resize(12,0);
	DimsXYZ.resize(3,0);

	lockedRecord.assign(2,0);
	SetPOI(0,0,0);
	SetExVector(1,0,0);
	SetEyVector(0,1,0);
	SetEzVector(0,0,1);

	SetDetectionGridDimensions(1,1,1,1,1,1);

	photonsPerLoop = 6400;
	numTest = 10000000;


}

CalculatePropagation::~CalculatePropagation() {
	// TODO Auto-generated destructor stub
	cout << "Destructing CalculatePropagation" << endl;


	inextInextp.~vector();
	MBMZ.~vector();
	g.~vector();
	mus.~vector();
	mua.~vector();
	n.~vector();
	triN.~vector();
	normals.~vector();
	triangles.~vector();
	dimensions.~vector();
	RecordAxz.~vector();
	RecordAyz.~vector();
	RandNmbersLog.~vector();
//	SurfRxy.~vector();
//	FetalRxy.~vector();
//	TransTxy.~vector();
	pStart.~vector();
	vStart.~vector();
	vStartOcc.~vector();

//	devices->~vector();
//	cout << "Delete PolyImage" << endl;
	//PolyImage->Delete();


	cout << "after deleting" << endl;





}

void CalculatePropagation::InitLogVectors(){
	if(verbose){cout << "Init LogVectors..."  << endl;}
	// PositionLog (aka Photontrace)

//	cout << "Dimensions: "<< dimensions[0] << ", " << dimensions[1] << ", " << dimensions[2]<< endl;


	int dimx = int(dimensions[0]/DetGridDims[0]);
	int dimy = int(dimensions[1]/DetGridDims[1]);
	int dimz = int(dimensions[2]/DetGridDims[2]);
	int dimAlpha = int(90/DetGridDims[5]);
	int dimR = int(ceil(sqrt(pow(dimx/2,2)+pow(dimy/2,2))));




	SurfFetalTransRxyImage = vtkSmartPointer<vtkImageData>::New();
	SurfFetalTransRxyImage->SetExtent(-dimx/2,dimx/2,-dimy/2,dimy/2,0,0);
	SurfFetalTransRxyImage->SetSpacing(1, 1, 1);
	SurfFetalTransRxyImage->SetOrigin(0, 0, 0);
	SurfFetalTransRxyImage->AllocateScalars(VTK_FLOAT,3);

	int ext[6];
	SurfFetalTransRxyImage->GetExtent(ext);
	if(verbose){
		cout << "Extent SFT: [" << ext[0] << ", " << ext[1]<<"],";
		cout << "[" << ext[2] << ", " << ext[3]<<"],";
		cout << "[" << ext[4] << ", " << ext[4]<<"]" << endl;
	}
	for(int i=ext[0];i<=ext[1];i++){
		for(int j=ext[2];j<=ext[3];j++){

			SurfFetalTransRxyImage->SetScalarComponentFromFloat(i,j,0,0,0);
			SurfFetalTransRxyImage->SetScalarComponentFromFloat(i,j,0,1,0);
			SurfFetalTransRxyImage->SetScalarComponentFromFloat(i,j,0,2,0);

		}
	}

	RTdxyImage = vtkSmartPointer<vtkImageData>::New();
	RTdxyImage->SetExtent(0,dimR-1,0,dimAlpha-1,0,0);
	RTdxyImage->SetSpacing(1, 1, 1);
	RTdxyImage->SetOrigin(0, 0, 0);
	RTdxyImage->AllocateScalars(VTK_FLOAT,3);

	RTdxyImage->GetExtent(ext);
	if(verbose){
		cout << "Extent RTdxyImage: [" << ext[0] << ", " << ext[1]<<"],";
		cout << "[" << ext[2] << ", " << ext[3]<<"],";
		cout << "[" << ext[4] << ", " << ext[4]<<"]" << endl;
	}
	for(int i=ext[0];i<=ext[1];i++){
		for(int j=ext[2];j<=ext[3];j++){

			RTdxyImage->SetScalarComponentFromFloat(i,j,0,0,0);//Rd
			RTdxyImage->SetScalarComponentFromFloat(i,j,0,1,0); //Fd
			RTdxyImage->SetScalarComponentFromFloat(i,j,0,2,0); //Td

		}
	}

	ArzImage = vtkSmartPointer<vtkImageData>::New();
	ArzImage->SetExtent(-dimR,dimR,-dimz,0,0,0);
	ArzImage->SetSpacing(1, 1, 1);
	ArzImage->SetOrigin(0, 0, 0);
	ArzImage->AllocateScalars(VTK_FLOAT,2);

	ArzImage->GetExtent(ext);

	if(verbose){
		cout << "Extent ArzImage: [" << ext[0] << ", " << ext[1]<<"],";
		cout << "[" << ext[2] << ", " << ext[3]<<"],";
		cout << "[" << ext[4] << ", " << ext[4]<<"]" << endl;
	}
		for(int i=ext[0];i<=ext[1];i++){
			for(int j=ext[2];j<=ext[3];j++){

				ArzImage->SetScalarComponentFromFloat(i,j,0,0,0);//Arz
				ArzImage->SetScalarComponentFromFloat(i,j,0,1,0); //Phirz


			}
		}

	SetDimensionsXYZ(dimx,dimy,dimz);

	// (phiBins*R)*3
	//int phiBins = int(360/DetGridDims[4]);
	//int R = floor(sqrt(pow(dimx/2,2)+pow(dimy/2,2))/DetGridDims[3]);


	//alphaRa.assign(dimAlpha,0);

	//SurfFetalTransRxy.assign(3*dimx*dimy,0);


	int dims[3];
	ArzImage->GetDimensions(dims);

	//int axyzDim = (dimx<dimy) ? (dimy+1)*dimz*2 : (dimx+1)*dimz*2;

	RecordAxzAyz.assign(3*dims[0]*dims[1],0);

//	cout << "RecordAxzAyz.size: " << RecordAxzAyz.size() << endl;
//	cout << "ArzImage Dims : " <<dims[0]<<", " << dims[1]<< endl;
//	RecordAxz.assign((dimx+1)*dimz,0);
//	RecordAyz.assign((dimy+1)*dimz,0);

//	RAxzCoords.assign(2*(dimx+1)*dimz,0);
//	RAyzCoords.assign(2*(dimy+1)*dimz,0);

//	RandNmbersLog.assign(100,0);

//	if(verbose){
	//cout << "Size SurfRxy: "<< SurfRxy.size() << endl;
//	cout << "Size FetalRxy: "<< FetalRxy.size() << endl;
//	cout << "Size RecordAxzAyz: "<< RecordAxzAyz.size() << endl;
//	cout << "Size RecordAxz: "<< RecordAxz.size() << endl;
//	cout << "Size RecordAyz: "<< RecordAyz.size() << endl;
///	cout << "Size RandNmbersLog: "<< RandNmbersLog.size() << endl;
//	}


}


/*
 * Setter functions
 */
void CalculatePropagation::SetVerbosity(bool verb){
	verbose = verb;
}
void CalculatePropagation::SetShowPercentDone(bool s){
	this->showPercent = s;
}

void CalculatePropagation::SetShowPercentDoneTrue(){
	this->showPercent = true;
}
void CalculatePropagation::SetShowPercentDoneFalse(){
	this->showPercent= false;
}

void CalculatePropagation::SetVerbosityTrue(){
	verbose = true;
}
void CalculatePropagation::SetVerbosityFalse(){
	verbose = false;
}
int CalculatePropagation::SetSourceString(const char* fname){

	status = convertToString(fname,sourceStr);
	if (status != 0) {
		std::cout << "Failed to open " << filename << std::endl;
		status = 1;
	}

	return status;

}

int CalculatePropagation::SetTestSourceString(const char* fname){

	status = convertToString(fname,testSourceStr);
	if (status != 0) {
		std::cout << "Failed to open " << fname << std::endl;
		status = 1;
	}

	return status;

}

void CalculatePropagation::SetBuildOptions(const char* options){

	this->buildOptions = options;

}

int CalculatePropagation::SetKernelName(const char* kname){

	kernelname = kname;

	return status;
}

void CalculatePropagation::SetDeviceType(int device){
// Type either "cpu" or "gpu"
// 0: CPU
// 1: GPU

	switch(device){
	case 0:

		devType = CL_DEVICE_TYPE_CPU;


		break;
	case 1:

	//	cout << "use GPU" << endl;
		devType = CL_DEVICE_TYPE_GPU;

		break;

	case 2:
		devType = CL_DEVICE_TYPE_DEFAULT;

		break;

	default:
		devType = CL_DEVICE_TYPE_DEFAULT;
		break;

	}

}




void CalculatePropagation::SetNumberOfPhotonInteractions(int pnum){

	this->numPhotonInteractions = pnum;
}
void CalculatePropagation::SetNumberOfPhotons(int pnum){

	this->absNumberPhotons = pnum;
	this->numberPhotons = pnum;
}

void CalculatePropagation::SetStartTrajectories(std::vector<float> &trajectories, std::vector<int> &trajOcc){
	this->vStart.assign(trajectories.begin(),trajectories.end());

	this->vStartOcc.assign(trajOcc.begin(),trajOcc.end());


}
void CalculatePropagation::SetStartTrajectory(std::vector<float> &trajectory){

	this->vStart.assign(trajectory.begin(),trajectory.end());
}

void CalculatePropagation::SetStartTrajectory(float x, float y, float z){
	this->vStart.assign(3,0);
	this->vStart[0] = x;
	this->vStart[1] = y;
	this->vStart[2] = z;

}

void CalculatePropagation::SetStartTrajectory(vtkDoubleArray* trajectory,int dims){

	std::vector<float> tr;
	float v;

	for(int i = 0;i<dims;i++){
		v = float(trajectory->GetValue(i));
		tr[i] = v;
	}

	this->vStart = tr;

}

void CalculatePropagation::SetStartPoint(std::vector<float> &startPoint){

	this->pStart.assign(startPoint.begin(),startPoint.end());

}

void CalculatePropagation::SetStartPoint(float x, float y, float z){

	this->pStart.assign(3,0);
	this->pStart[0] = x;
	this->pStart[1] = y;
	this->pStart[2] = z;

}

void CalculatePropagation::SetStartPoint(vtkDoubleArray* startPoint,int dims){

	std::vector<float> sp;
	float v;

	for(int i = 0;i<dims;i++){
		v = float(startPoint->GetValue(i));
		sp[i] = v;
	}

	this->pStart.assign(sp.begin(),sp.end());

}


void CalculatePropagation::SetChanceForRoulette(float c){
	this->chance = c;
}

void CalculatePropagation::SetWeightThreshold(float w){
	this->w_threshold = w;
}
void CalculatePropagation::SetPixelPerUnitFactor(float factor) {
	this->pixelPerUnitFactor = factor;
}


void CalculatePropagation::SetTissueProperties(std::vector<float> &refractionInd, std::vector<float> &absorbtionCoeffs, std::vector<float> &anisotropies, std::vector<float> &scatterCoeffs){
/*
	this->n.assign(refractionInd.begin(),refractionInd.end());
	this->mua.assign(absorbtionCoeffs.begin(),absorbtionCoeffs.end());
	this->mus.assign(scatterCoeffs.begin(),scatterCoeffs.end());
	this->g.assign(anisotropies.begin(),anisotropies.end());
	*/

	int s = refractionInd.size();

	this->nmuasG.assign(4*s,0);

	for(int i=0;i<s;i++){


		this->nmuasG[4*i] = refractionInd[i];
		this->nmuasG[4*i+1] = absorbtionCoeffs[i];
		this->nmuasG[4*i+2] = scatterCoeffs[i];
		this->nmuasG[4*i+3] = anisotropies[i];
	}




}

void CalculatePropagation::SetRefractionIndices(std::vector<float> &ind){

	this->n.assign(ind.begin(),ind.end());
}

void CalculatePropagation::SetRefractionIndices(vtkDoubleArray* refractions,int dims){

	std::vector<float> ind;
	float v;

	for(int i = 0;i<dims;i++){
		v = float(refractions->GetValue(i));
		ind[i] = v;
	}

	this->n = ind;
}

void CalculatePropagation::SetRefractionIndices(double* refractions,int dims){

	std::vector<float> ind;

	for(int i = 0;i<dims;i++){

		ind[i] = float(refractions[i]);
	}

	this->n = ind;
}

void CalculatePropagation::SetAbsorbtioncoefficients(std::vector<float> &absorbtions){
	this->mua.assign(absorbtions.begin(),absorbtions.end());
}

void CalculatePropagation::SetAbsorbtioncoefficients(vtkDoubleArray* absorbtions,int dims){
	std::vector<float> mua;
	float v;

	for(int i = 0;i<dims;i++){
		v = float(absorbtions->GetValue(i));
		mua[i] = v;
	}


	this->mua = mua;
}

void CalculatePropagation::SetScattercoefficients(std::vector<float> &scatter){
	this->mus.assign(scatter.begin(),scatter.end());
}

void CalculatePropagation::SetScattercoefficients(vtkDoubleArray* scatter,int dims){

	std::vector<float> mus;
		float v;

		for(int i = 0;i<dims;i++){
			v = float(scatter->GetValue(i));
			mus[i] = v;
		}

	this->mus = mus;
}

void CalculatePropagation::SetAnisotropy(std::vector<float> &anisotropies){
	this->g.assign(anisotropies.begin(),anisotropies.end());
}

void CalculatePropagation::SetAnisotropy(vtkDoubleArray* anisotropies,int dims){
	std::vector<float> g;
			float v;

			for(int i = 0;i<dims;i++){
				v = float(anisotropies->GetValue(i));
				mus[i] = v;
			}

	this->g = g;

}


void CalculatePropagation::SetDimensions(std::vector<float> &dims){

	dimensions.assign(dims.begin(),dims.end());


}
void CalculatePropagation::SetTriangles(vtkPolyData* meshimage){
	//Triangles: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*9



	int numCells = meshimage->GetNumberOfCells();
	std::vector<float> points;
	points.assign(numCells*9,0);


	for(int i = 0;i<numCells;i++){

		double p0[3];
		double p1[3];
		double p2[3];


		meshimage->GetCell(i)->GetPoints()->GetPoint(0,p0);
		meshimage->GetCell(i)->GetPoints()->GetPoint(1,p1);
		meshimage->GetCell(i)->GetPoints()->GetPoint(2,p2);


		points[9*i] = float(p0[0]); //p0_x
		points[9*i+1] = float(p0[1]); //p0_y
		points[9*i+2] = float(p0[2]); //p0_z

		points[9*i+3] = float(p1[0]); //p1_x
		points[9*i+4] = float(p1[1]); //p1_y
		points[9*i+5] = float(p1[2]); //p1_z

		points[9*i+6] = float(p2[0]); //p2_x
		points[9*i+7] = float(p2[1]); //p2_y
		points[9*i+8] = float(p2[2]); //p2_z

	}

	this->triangles.assign(points.begin(),points.end());

	//printVector("Triangles: ",this->triangles);

}
void CalculatePropagation::SetNormals(vtkDataArray* normals){
	 // Normals per Triangle: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3

	int num = normals->GetSize();
		std::vector<float> normal;
		normal.assign(num*3,0);

		for(int i = 0;i<num;i++){

			double n0[3];

			normals->GetTuple(i,n0);

			normal[3*i] = float(n0[0]); //n0_x
			normal[3*i+1] = float(n0[1]); //n0_y
			normal[3*i+2] = float(n0[2]); //n0_z

		}


		this->normals.assign(normal.begin(),normal.end());
}

void CalculatePropagation::SetInnerAndOuterTissueIndices(vtkDoubleArray* mediums, double* contours, int numContours){
	// index of n and n_outer for every triangle {(n0_in,n0_out),(n1_in,n1_out),...}, Size: NumTriangles*2;

	int num = mediums->GetNumberOfTuples();
	std::vector<int> medium;
	if(verbose){
	cout << "in SetInnerOuter: NumTuples: " << num << endl;
	}
	medium.assign(2*num,0);


	for(int i = 0;i<num;i++){

		double n0[2];
		int indInner, indOuter;

		mediums->GetTuple(i,n0);


		FindIndex(n0[0],contours,numContours,indInner);
		FindIndex(n0[1],contours,numContours,indOuter);


		medium[2*i] = indInner; //n0_inner
		medium[2*i+1] = indOuter; //n0_outer


	}

	this->triN.assign(medium.begin(),medium.end());

}


void CalculatePropagation::SetInnerAndOuterTissueIndices(vtkDoubleArray* mediums){
	// index of n and n_outer for every triangle {(n0_in,n0_out),(n1_in,n1_out),...}, Size: NumTriangles*2;

	int num = mediums->GetNumberOfTuples();
	std::vector<int> medium;
	if(verbose){
	cout << "in SetInnerOuter: NumTuples: " << num << endl;
	}
	medium.assign(2*num,0);


	for(int i = 0;i<num;i++){

		double n0[2];
	//	int indInner, indOuter;

		mediums->GetTuple(i,n0);



		medium[2*i] = int(n0[0]); //n0_inner
		medium[2*i+1] = int(n0[1]); //n0_outer


	}

	this->triN.assign(medium.begin(),medium.end());

}

void CalculatePropagation::FindIndex(double value,double* valuesArray,int arrSize, int &index){
	// Looks for value in ValuesArray and sets index to index where values was found

	for(int i = 0; i<arrSize;i++){

			if(value == valuesArray[i]){
				index = i;
				break;
			}else{
			    index = 0;
			}
	}
}

void CalculatePropagation::SetRand3Array(std::vector<float> &rands){

	this->rand3Arr.assign(rands.begin(),rands.end());

}

void CalculatePropagation::SetRand3Array(vtkDoubleArray* rands,int dims){
	std::vector<float> r;
			float v;

			for(int i = 0;i<dims;i++){
				v = float(rands->GetValue(i));
				r[i] = v;
			}

	this->rand3Arr = r;

}


void CalculatePropagation::SetSeed(int s){

	this->seed = s;

}

void CalculatePropagation::SetNumberGenerator(int psnr){
	this->numberGenerator = psnr;
}

void CalculatePropagation::SetDimensionsXYZ(float xdim, float ydim, float zdim){
	DimsXYZ[0] = xdim;
	DimsXYZ[1] = ydim;
	DimsXYZ[2] = zdim;

}
void CalculatePropagation::SetDimensionsXYZ(std::vector<float> &dims){
	DimsXYZ.assign(dims.begin(),dims.end());
}

void CalculatePropagation::SetDetectionGridDimensions(float xdim, float ydim, float zdim, float rdim, float phiDim, float rhoDim){

	DetGridDims[0] = xdim;
	DetGridDims[1] = ydim;
	DetGridDims[2] = zdim;
	DetGridDims[3] = rdim;
	DetGridDims[4] = phiDim;
	DetGridDims[5] = rhoDim;

}

void CalculatePropagation::SetDetectionGridDimensions(std::vector<float> &detGrid){
	DetGridDims.assign(detGrid.begin(),detGrid.end());
}

void CalculatePropagation::SetPOI(float x, float y, float z){
	eP0XYZ[0] = x;
	eP0XYZ[1] = y;
	eP0XYZ[2] = z;
}

void CalculatePropagation::SetPOI(std::vector<float> &point){
	eP0XYZ[0] = point[0];
	eP0XYZ[1] = point[1];
	eP0XYZ[2] = point[2];
}

void CalculatePropagation::SetBasisVectors(std::vector<float> &veX,std::vector<float> &veY,std::vector<float> &veZ){

	SetExVector(veX[0],veX[1],veX[2]);
	SetEyVector(veY[0],veY[1],veY[2]);
	SetEzVector(veZ[0],veZ[1],veZ[2]);


}

void CalculatePropagation::SetExVector(float x,float y, float z){
	eP0XYZ[3] = x;
	eP0XYZ[4] = y;
	eP0XYZ[5] = z;
}
void CalculatePropagation::SetEyVector(float x,float y, float z){
	eP0XYZ[6] = x;
	eP0XYZ[7] = y;
	eP0XYZ[8] = z;
}
void CalculatePropagation::SetEzVector(float x,float y, float z){
	eP0XYZ[9] = x;
	eP0XYZ[10] = y;
	eP0XYZ[11] = z;
}

void CalculatePropagation::SetFetalTissueID(int id){
	fetalTissueID = id;
}

/*
 * Getter functions
 */

int CalculatePropagation::GetNumberOfPhotons(){
	return this->numberPhotons;
}

std::vector<float> CalculatePropagation::GetStartTrajectory(){

	return this->vStart;
}

void CalculatePropagation::GetStartTrajectory(std::vector<float> &trajectory){

	trajectory.assign(vStart.begin(),vStart.end());

}

std::vector<float> CalculatePropagation::GetStartPoint(){

	return this->pStart;
}

void CalculatePropagation::GetStartPoint(std::vector<float> &startPoint){

	startPoint.assign(pStart.begin(),pStart.end());

}

float CalculatePropagation::GetChanceForRoulette(){

	return this->chance;
}

float CalculatePropagation::GetWeightThreshold(){

	return this->w_threshold;
}

std::vector<float> CalculatePropagation::GetRefractionindices(){

	return this->n;
}

void CalculatePropagation::GetRefractionindices(std::vector<float> &refractions){

	refractions.assign(n.begin(),n.end());

}

std::vector<float> CalculatePropagation::GetAbsorbtioncoefficients(){

	return this->mua;
}

void CalculatePropagation::GetAbsorbtioncoefficients(std::vector<float> &absorbtions){

	absorbtions.assign(mua.begin(),mua.end());

}

std::vector<float> CalculatePropagation::GetScattercoefficients(){

	return this->mus;
}

void CalculatePropagation::GetScattercoefficients(std::vector<float> &scatter){

	scatter.assign(mus.begin(),mus.end());

}

std::vector<float> CalculatePropagation::GetAnisotropy(){

	return this->g;
}

void CalculatePropagation::GetAnisotropy(std::vector<float> &Anisotropy){

	Anisotropy.assign(g.begin(),g.end());

}



// Get Logs
/*
void CalculatePropagation::GetLogs(std::vector<float> &randoms,
		std::vector<float> &absorbtionsXZ,std::vector<float> &absorbtionsYZ,
		std::vector<float> &surfaceXY,std::vector<float> &fetalXY,
		std::vector<float> &transmittenceXY){

	randoms.assign(RandNmbersLog.begin(),RandNmbersLog.end());
	absorbtionsXZ.assign(RecordAxz.begin(),RecordAxz.end());
	absorbtionsYZ.assign(RecordAyz.begin(),RecordAyz.end());
	surfaceXY.assign(SurfRxy.begin(),SurfRxy.end());
	transmittenceXY.assign(TransTxy.begin(),TransTxy.end());
	fetalXY.assign(FetalRxy.begin(),FetalRxy.end());

}

void CalculatePropagation::GetLogs(std::vector<float> &randoms,
		std::vector<float> &absorbtionsXZ,std::vector<float> &absorbtionsYZ,
		std::vector<float> &surfFetalTransXY){


	randoms.assign(RandNmbersLog.begin(),RandNmbersLog.end());
	absorbtionsXZ.assign(RecordAxz.begin(),RecordAxz.end());
	absorbtionsYZ.assign(RecordAyz.begin(),RecordAyz.end());
	//surfFetalTransXY.assign(SurfFetalTransRxy.begin(),SurfFetalTransRxy.end());


}


*/

vtkSmartPointer<vtkImageData> CalculatePropagation::GetSFTImage(){
	return SurfFetalTransRxyImage;
}

vtkSmartPointer<vtkImageData> CalculatePropagation::GetSFTImage(int component){
	int ext[6];


	SurfFetalTransRxyImage->GetExtent(ext);
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

	image->SetExtent(ext);
	image->SetOrigin(SurfFetalTransRxyImage->GetOrigin());
	image->SetSpacing(SurfFetalTransRxyImage->GetSpacing());
	image->AllocateScalars(VTK_FLOAT,1);


	for(int i=ext[0];i<=ext[1];i++){
			for(int j=ext[2];j<=ext[3];j++){

				image->SetScalarComponentFromFloat(i,j,0,0,SurfFetalTransRxyImage->GetScalarComponentAsFloat(i,j,0,component));


			}
		}


	return image;
}

vtkSmartPointer<vtkImageData> CalculatePropagation::GetArzImage(){
	int ext[6];
// not scaled yet

	ArzImage->GetExtent(ext);
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

	image->SetExtent(ext);
	image->SetOrigin(ArzImage->GetOrigin());
	image->SetSpacing(ArzImage->GetSpacing());
	image->AllocateScalars(VTK_FLOAT,1);


	for(int r=ext[0];r<=ext[1];r++){
		for(int z=ext[2];z<=ext[3];z++){
			float v = ArzImage->GetScalarComponentAsFloat(r,z,0,0);
			image->SetScalarComponentFromFloat(r,z,0,0,v);

		}
	}


	return image;
}







vtkSmartPointer<vtkImageData> CalculatePropagation::GetPhirzImage(){
	int ext[6];
	// not scaled yet

	ArzImage->GetExtent(ext);
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

	image->SetExtent(ext);
	image->SetOrigin(ArzImage->GetOrigin());
	image->SetSpacing(ArzImage->GetSpacing());
	image->AllocateScalars(VTK_FLOAT,1);


	for(int i=ext[0];i<=ext[1];i++){
			for(int j=ext[2];j<=ext[3];j++){
				float v = ArzImage->GetScalarComponentAsFloat(i,j,0,1);
				image->SetScalarComponentFromFloat(i,j,0,0,v);


			}
		}


	return image;
}

float CalculatePropagation::ScalePerAngle1D(int ia){


	float scale = 1;
	float da = DetGridDims[5]; // in [°]

	scale  = 2.0*M_PI*da*numberPhotons;

	return 1/(sin((ia+0.5)*da)*scale);

}

float CalculatePropagation::Scale2D(int ir, int ia){

	 double dr = DetGridDims[3];
	 double da = DetGridDims[5]*M_PI/180; // in [rad]

	float scale1 = 4.0*M_PI*M_PI*dr*sin(da/2)*dr*numberPhotons;
		/* The factor (ir+0.5)*sin(2a) to be added. */
	float scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);

return scale2;

}

float CalculatePropagation::ScaleArz2D(int ir){
// Usag: /=scale;


	  double dz = DetGridDims[2];
	  double dr = DetGridDims[3];

	  double scale1;

	  /* Scale A_rz. */
	  scale1 = 2.0*M_PI*dr*dr*dz*numberPhotons;
		/* volume is 2*pi*(ir+0.5)*dr*dr*dz.*/
	  return (ir+0.5)*scale1;

}


vtkSmartPointer<vtkImageData> CalculatePropagation::GetRTdImage(){
	return RTdxyImage;
}

vtkSmartPointer<vtkImageData> CalculatePropagation::GetRTdImage(int component){
	int ext[6];


	RTdxyImage->GetExtent(ext);
	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

	image->SetExtent(ext);
	image->SetOrigin(RTdxyImage->GetOrigin());
	image->SetSpacing(RTdxyImage->GetSpacing());
	image->AllocateScalars(VTK_FLOAT,1);


	for(int i=ext[0];i<=ext[1];i++){
			for(int j=ext[2];j<=ext[3];j++){

				float value = RTdxyImage->GetScalarComponentAsFloat(i,j,0,component);

				image->SetScalarComponentFromFloat(i,j,0,0,value);


			}
		}


	return image;
}



void CalculatePropagation::GetAbsorbtionLog(std::vector<float> &AXZ, std::vector<float> &AYZ){

	AXZ.assign(RecordAxz.begin(),RecordAxz.end());
	AYZ.assign(RecordAyz.begin(),RecordAyz.end());

}
void CalculatePropagation::GetAbsorbtionLogCoords(std::vector<int> &XZcoords, std::vector<int> &YZcoords){

	XZcoords.assign(RAxzCoords.begin(),RAxzCoords.end());
	YZcoords.assign(RAyzCoords.begin(),RAyzCoords.end());

}

std::vector<float> CalculatePropagation::GetRandomNumbersLog(){

	return this->RandNmbersLog;
}

std::vector<float> CalculatePropagation::GetRand3Array(){

	return this->rand3Arr;
}

void CalculatePropagation::GetRand3Array(std::vector<float> &psnrs){

	psnrs.assign(rand3Arr.begin(),rand3Arr.end());

}

int CalculatePropagation::GetSeed(){

	return this->seed;
}
int CalculatePropagation::GetNumberGenerator(){

	return this->numberGenerator;
}


void FillArray(int num, std::vector<int> &arr){

	for(int i=0;i<num;i++){

		arr.insert(arr.end(),i);


	}

}

void CalculatePropagation::CreateTestBuffer(){


	FillArray(numTest,A);
	FillArray(numTest,B);
	C.assign(numTest,0);


	d_A = cl::Buffer(context,A.begin(),A.end(),isReadOnly);
	d_B = cl::Buffer(context,B.begin(),B.end(),isReadOnly);
	d_C = cl::Buffer(context,C.begin(),C.end(),!isReadOnly);


}

void CalculatePropagation::CreateBufferObjects(){


	//Photon properties
	 d_vStart = cl::Buffer(vStart.begin(), vStart.end(),isReadOnly);
	// d_pStart = cl::Buffer(pStart.begin(), pStart.end(),isReadOnly);
	 d_vStartOcc = cl::Buffer(vStartOcc.begin(), vStartOcc.end(),isReadOnly);
 	 d_DetGridDims = cl::Buffer(DetGridDims.begin(), DetGridDims.end(),isReadOnly);
	 d_eP0XYZ = cl::Buffer(eP0XYZ.begin(), eP0XYZ.end(),isReadOnly);
	 d_pStart = {pStart[0],pStart[1], pStart[2]};
	//Logs

	 d_SurfFetalTransRxy = cl::Buffer(CL_MEM_READ_WRITE, photonsPerLoop*sizeof(cl_float8),NULL,NULL);
//	 d_TransTxy = cl::Buffer(TransTxy.begin(), TransTxy.end(),!isReadOnly);
//	 d_SurfRxy = cl::Buffer(SurfRxy.begin(), SurfRxy.end(),!isReadOnly);
//	 d_Fetalxy = cl::Buffer(FetalRxy.begin(), FetalRxy.end(),!isReadOnly);

	 //	 d_RecordAxz = cl::Buffer(RecordAxz.begin(), RecordAxz.end(),!isReadOnly);
	 //d_RecordAyz = cl::Buffer(RecordAyz.begin(), RecordAyz.end(),!isReadOnly);
	 d_RandNmbersLog = cl::Buffer(RandNmbersLog.begin(), RandNmbersLog.end(),!isReadOnly);

	// Tissue properties and geometry

	 d_triangles = cl::Buffer(triangles.begin(), triangles.end(),isReadOnly);
	 d_normals = cl::Buffer(normals.begin(), normals.end(),isReadOnly);
	 d_triN = cl::Buffer(triN.begin(), triN.end(),isReadOnly);
//	 d_n = cl::Buffer(n.begin(), n.end(),isReadOnly);
//	 d_mua = cl::Buffer(mua.begin(), mua.end(),isReadOnly);
//	 d_mus = cl::Buffer(mus.begin(), mus.end(),isReadOnly);
//	 d_g = cl::Buffer(g.begin(), g.end(),isReadOnly);

	 d_nmuasG = cl::Buffer(nmuasG.begin(),nmuasG.end(),isReadOnly);
	// PRNG
	 d_rand3Array = cl::Buffer(rand3Arr.begin(), rand3Arr.end(),!isReadOnly);
	 d_randPool = cl::Buffer(CL_MEM_READ_WRITE,poolSize*sizeof(cl_float),NULL,NULL);

	 d_MBMZ = cl::Buffer(MBMZ.begin(), MBMZ.end(),isReadOnly);
	 d_inextInextp = cl::Buffer(inextInextp.begin(), inextInextp.end(),!isReadOnly);


}

void CalculatePropagation::CreateBufferObjects(cl::Context con){

	if(verbose){cout << "Create BufferObjects..."  << endl;}
	int dimx = int(dimensions[0]/DetGridDims[0]);
	int dimy = int(dimensions[1]/DetGridDims[1]);
	int dimz = int(dimensions[2]/DetGridDims[2]);
	int dimR = int(ceil(sqrt(pow(dimx/2,2)+pow(dimy/2,2))));


	//Photon properties
		 d_vStart = cl::Buffer(con,vStart.begin(), vStart.end(),isReadOnly,NULL,NULL);
	//	 d_pStart = cl::Buffer(con,pStart.begin(), pStart.end(),isReadOnly,NULL,NULL);
		 d_vStartOcc = cl::Buffer(con,vStartOcc.begin(), vStartOcc.end(),isReadOnly,NULL,NULL);
	 	 d_DetGridDims = cl::Buffer(con,DetGridDims.begin(), DetGridDims.end(),isReadOnly,NULL,NULL);
		 d_eP0XYZ = cl::Buffer(con,eP0XYZ.begin(), eP0XYZ.end(),isReadOnly,NULL,NULL);
		 d_DimsXYZ = cl::Buffer(con,DimsXYZ.begin(), DimsXYZ.end(),isReadOnly,NULL,NULL);

		 d_pStart = {pStart[0],pStart[1], pStart[2]};

//cout << "befor Logs..."<< endl;

		//Logs
		 int s = RecordAxzAyz.size();

		 d_SurfFetalTransRxy = cl::Buffer(con,CL_MEM_READ_WRITE, photonsPerLoop*sizeof(cl_float8),NULL,NULL);
		 //d_SurfFetalTransRxy = cl::Buffer(con,SurfFetalTransRxy.begin(), SurfFetalTransRxy.end(),!isReadOnly,true,NULL);

		 d_RecordAxzAyz = cl::Buffer(con,RecordAxzAyz.begin(),RecordAxzAyz.end(),!isReadOnly,NULL,NULL);
		// d_RecordAxzAyz = cl::Buffer(con,CL_MEM_READ_WRITE, s*sizeof(cl_float2),NULL,NULL);
//		 d_TransTxy = cl::Buffer(con,TransTxy.begin(), TransTxy.end(),!isReadOnly,true,NULL);
//		 d_SurfRxy = cl::Buffer(con,SurfRxy.begin(), SurfRxy.end(),!isReadOnly,true,NULL);
//		 d_Fetalxy = cl::Buffer(con,FetalRxy.begin(), FetalRxy.end(),!isReadOnly,true,NULL);
	//	 d_RecordAxz = cl::Buffer(con,RecordAxz.begin(), RecordAxz.end(),!isReadOnly,true,NULL);
	//	 d_RecordAyz = cl::Buffer(con,RecordAyz.begin(), RecordAyz.end(),!isReadOnly,true,NULL);
//	 d_RandNmbersLog = cl::Buffer(con,RandNmbersLog.begin(), RandNmbersLog.end(),!isReadOnly,true,NULL);

	//	 d_lockedRecord = cl::Buffer(con,lockedRecord.begin(), lockedRecord.end(),!isReadOnly,true,NULL);
		// Tissue properties and geometry

		 d_triangles = cl::Buffer(con,triangles.begin(), triangles.end(),isReadOnly,false,NULL);
		 d_normals = cl::Buffer(con,normals.begin(), normals.end(),isReadOnly,false,NULL);
		 d_triN = cl::Buffer(con,triN.begin(), triN.end(),isReadOnly,false,NULL);
//		 d_n = cl::Buffer(con,n.begin(), n.end(),isReadOnly,false,NULL);
//		 d_mua = cl::Buffer(con,mua.begin(), mua.end(),isReadOnly,false,NULL);
//		 d_mus = cl::Buffer(con,mus.begin(), mus.end(),isReadOnly,false,NULL);
//		 d_g = cl::Buffer(con,g.begin(), g.end(),isReadOnly,false,NULL);
		 d_nmuasG = cl::Buffer(con,nmuasG.begin(), nmuasG.end(),isReadOnly,false,NULL);
		// PRNG
		 d_rand3Array = cl::Buffer(con,rand3Arr.begin(), rand3Arr.end(),!isReadOnly,true,NULL);

		 // d_randPool = cl::Buffer(con,CL_MEM_READ_WRITE,poolSize*sizeof(cl_float),NULL,NULL);

		 d_MBMZ = cl::Buffer(con,MBMZ.begin(), MBMZ.end(),isReadOnly,false,NULL);
		 d_inextInextp = cl::Buffer(con,inextInextp.begin(), inextInextp.end(),!isReadOnly,NULL,NULL);

}


int CalculatePropagation::CreateCommandQueue(){

	//cl_command_queue_properties cqpss = CL_QUEUE_PROFILING_ENABLE;

	cmd_queue = cl::CommandQueue(context,selectedDevice,CL_QUEUE_PROFILING_ENABLE,NULL);

	testcmd_queue = cl::CommandQueue(context,selectedDevice,CL_QUEUE_PROFILING_ENABLE | CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE,NULL);


	//cmd_queue.enqueueReadBuffer(d_SurfRxy,CL_TRUE,0,sizeof(cl_float)*SurfRxy.size(),SurfRxy);



	return SUCCESS;
}

cl::Platform CalculatePropagation::getPlatform(cl_device_type type) {
    // Get available platforms

	status = cl::Platform::get(&platforms);

	//cout << "Number of Platforms: " << platforms.size() << endl;

    if(platforms.size() == 0)
        throw cl::Error(1, "No OpenCL platforms were found");

    int platformID = -1;
    int deviceID = -1;

        for(unsigned int i = 0; i < platforms.size(); i++) {
        	try {

                platforms[i].getDevices(type, &devices);
                platformID = i;

                int numd = devices.size();
              //  cout << "Number of Devices: " << numd << endl;
                for(int j=0;j<numd;j++){
                //	cout << "DEvice Type: " << devices[j].getInfo<CL_DEVICE_TYPE>() << endl;

                	if(devices[j].getInfo<CL_DEVICE_TYPE>()==type){

                		selectedDevice = devices[j];
                		deviceID = j;
                		//cout << "Selected Device Type: " << selectedDevice.getInfo<CL_DEVICE_TYPE>() << endl;
                	}
                }


                break;
            } catch(cl::Error e) {
               continue;
            }
        }


    if(platformID == -1)
        throw cl::Error(1, "No compatible OpenCL platform found");

    cl::Platform platform = platforms[platformID];
//    std::cout << "Using platform vendor: " << platform.getInfo<CL_PLATFORM_VENDOR>() << std::endl;
//    std::cout << "Using platform name: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
  //  cout << "Device Type(" << deviceID << "): " << devices[deviceID].getInfo<CL_DEVICE_TYPE>() << endl;


//    cout << "CL_DEVICE_TYPE_ALL:" << CL_DEVICE_TYPE_ALL << endl;
//    cout << "CL_DEVICE_TYPE_CPU:" << CL_DEVICE_TYPE_CPU << endl;
//    cout << "CL_DEVICE_TYPE_GPU:" << CL_DEVICE_TYPE_GPU << endl;
//    cout << "CL_DEVICE_TYPE_CUSTOM:" << CL_DEVICE_TYPE_CUSTOM << endl;
//    cout << "CL_DEVICE_TYPE_DEFAULT:" << CL_DEVICE_TYPE_DEFAULT << endl;
//    cout << "CL_DEVICE_TYPE_ACCELERATOR:" << CL_DEVICE_TYPE_ACCELERATOR << endl;



    return platform;
}

int CalculatePropagation::createContext(){
	if(verbose){cout << "create Context..."  << endl;}
	 cl_int status = 0;

	 cl::Platform platform = getPlatform(devType);

	 //platform.getDevices(dType,&devices);

//cout << "dType: " << devType << endl;

	 // Use the preferred platform and create a context
	    cl_context_properties cps[] = {
	        CL_CONTEXT_PLATFORM,
	        (cl_context_properties)(platform)(),
	        0
	    };

	    try {
	        context = cl::Context(devType,cps);

	        this->isContextSet = true;

	        //cout << "COntext Devices: " << context.getInfo<CL_CONTEXT_DEVICES>() << endl;

	        //cout << "(GPU: " << CL_DEVICE_TYPE_GPU << ", CPU: " << .getInfo<CL_DEVICE_N>() << ")" << endl;

	    } catch(cl::Error error) {
	        throw cl::Error(1, "Failed to create an OpenCL context!");
	    }


	 //   wgSize = selectedDevice.getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();




	    return status;
}

int CalculatePropagation::BuildKernel(){

	if(verbose){cout << "Build Kernel..."  << endl;}
	/*
	 * create and build program
	 */
	//cl::Error e;
		SetSourceString(filename);
		SetTestSourceString("vecAdd.cl");

		if(isContextSet){
			program = cl::Program(context,sourceStr);
			Testprogram = cl::Program(context,testSourceStr);
			devices = context.getInfo<CL_CONTEXT_DEVICES>();
		}else{
			program = cl::Program(sourceStr);
		}

//		for(auto &arr : devices){
//
//			cout << "Device name: " << arr.getInfo<CL_DEVICE_NAME>() << endl;
//
//		}


		try{

			program.build(devices,buildOptions); //for build options
			Testprogram.build(devices,buildOptions);

		}catch (cl::Error e) {
			std::cout << e.what() << std::endl;
			std::cout << "Build Status: " << program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(cl::Device::getDefault()) << std::endl;
			std::cout << "Build Options:\t" << program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(cl::Device::getDefault()) << std::endl;
			std::cout << "Build Log:\t " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(cl::Device::getDefault()) << std::endl;

			return FAILURE;

		}



	//	photonsPerLoop=64*8;
	//	cout << "Build done..." << endl;

		return SUCCESS;
}

int CalculatePropagation::ExecKernel(){

	if(verbose){cout << "Exec Kernel..."  << endl;}
/*
	int dimx = int(dimensions[0]/DetGridDims[0]);
	int dimy = int(dimensions[1]/DetGridDims[1]);
	int dimz = int(dimensions[2]/DetGridDims[2]);
*/




	cl::Event e;
	e = 0;
	cl::Event eTest;

	  // set arguments for kernel, and execute it.
		//	    cl::NDRange ndrg(numberPhotons);
			    //cl::NDRange ndrl(64);



			    cl::EnqueueArgs argTest(testcmd_queue,cl::NDRange(numTest));


	// create kernel as a functor
	//kernel = cl::Kernel(program,kernelname);

	KernelType MCKernel(program,kernelname);
	TestKernel Adder(Testprogram,"addition");



//	int info = devices[0].getInfo<CL_DEVICE_MAX_WORK_GROUP_SIZE>();
//	cout << "WGS: " << info << endl;


		    cl_int numtris = int (triangles.size()/9);
		    if(verbose){
		    cout << "Number od Triangles: " << numtris << endl;
		    cout << "OpenCL-File: " << this->filename << endl;
		    cout << "OCL Kernelname: " << this->kernelname << endl;
		    printVector("vStart: ", vStart);
		    printVector("pStart: ", pStart);
		    cout << "Size Triangles: " << triangles.size() << endl;
		    cout << "Size N: " << nmuasG.size()/4 << endl;
		    cout << "Size triN: " << triN.size() << endl;
		    cout << "Size rand3Arr: " << rand3Arr.size() << endl;
		    cout << "Size RecordAxzyz: "<< RecordAxzAyz.size() << endl;
		    }



//kernel.getWorkGroupInfo(selectedDevice,CL_KERNEL_WORK_GROUP_SIZE,&lWGSize);
//kernel.getWorkGroupInfo(selectedDevice,CL_KERNEL_GLOBAL_WORK_SIZE,&gWSize);
		    neededTime = 0;

		    int numLoops = int(numberPhotons/photonsPerLoop);

		    if(numberPhotons%photonsPerLoop !=0){
		    	numLoops++;
		    }

		    cl::NDRange ndrg(numberPhotons);
		  //  cl::NDRange ndrl(128);



		    std::vector<float> tempSFT;
		    std::vector<float> tempAxyz;
		//    std::vector<float> tempRands;
		 //   tempRands.assign(56,0);

		  //  C.assign(numTest,0);
		  // tempSFT.assign(photonsPerLoop*4,0);



		    for(int i=0;i<numLoops;i++){

		    	if(photonsPerLoop >= numberPhotons){

		    		ndrg = numberPhotons;
		    		tempSFT.assign(numberPhotons*8,0);


		    	}else if(numLoops>1 && i==(numLoops-1)){

		    		int photons = numberPhotons-i*photonsPerLoop;
		    		ndrg = photons;
		    		tempSFT.assign(photons*8,0);


		    	}else {

		    		ndrg = photonsPerLoop;
		    		tempSFT.assign(photonsPerLoop*8,0);

		    	}


		    	int dimx = int(dimensions[0]/DetGridDims[0]);
		    	int dimy = int(dimensions[1]/DetGridDims[1]);
		    	int dimz = int(dimensions[2]/DetGridDims[2]);
		    	int dimR = int(ceil(sqrt(pow(dimx/2,2)+pow(dimy/2,2))));
		    	int dims[3];
		    	ArzImage->GetDimensions(dims);

		    		//int axyzDim = (dimx<dimy) ? (dimy+1)*dimz*2 : (dimx+1)*dimz*2;


		    	tempAxyz.assign(2*dims[0]*dims[1],0);


		    	 //cl::NDRange ndrl(64);
		    	//cl::NDRange offset = i*photonsPerLoop;
		    	cl::EnqueueArgs arg(cmd_queue,ndrg);
		    	// execute the kernel by calling the kernel functor

		    	 if(showPercent){

		    		 PrintPercentDone(double(i*100)/numLoops);
		    	 }




		    	e = MCKernel(arg,
		    			numtris,
		    			pixelPerUnitFactor,
		    			d_vStart,
		    			d_vStartOcc,
		    			d_pStart,
		    			chance,
		    			w_threshold,
		    			fetalTissueID,
		    			d_SurfFetalTransRxy,
		    			d_RecordAxzAyz,
		    			d_DimsXYZ,
		    			d_eP0XYZ,
		    			d_DetGridDims,
		    			d_triangles,
		    			d_normals,
		    			d_triN,
		    			d_nmuasG,
		    			d_rand3Array,
		    			d_inextInextp);

		    	e.wait();



		    cl::copy(cmd_queue,d_SurfFetalTransRxy,tempSFT.begin(), tempSFT.end());




		    	e.getProfilingInfo(CL_PROFILING_COMMAND_START,&startTime);
		    	e.getProfilingInfo(CL_PROFILING_COMMAND_END,&endTime);
		    	neededTime +=(endTime-startTime);
		    	int s = ndrg.operator const unsigned long int *()[0];
		    	int ex[6];
		    	SurfFetalTransRxyImage->GetExtent(ex);

		    	for(int j = 0;j<s;j++){

		    		int x,y,z;


		    		x = int(tempSFT[8*j]/DetGridDims[0]);
		    		y = int(tempSFT[8*j+1]/DetGridDims[1]);
		    		z = int(tempSFT[8*j+2]/DetGridDims[2]);



		    		y = (y<ex[2]) ? ex[2]:y;
		    		y = (y>ex[3]) ? ex[3]:y;


		    		if(x<ex[0]){
		    			 x = ex[0];
		    		}else if(x>ex[1]){
		    			x = ex[1];
		    		}

		    //		cout << "SFT x: " << x << ", y: " << y << endl;

		    		float surfRxyVal = SurfFetalTransRxyImage->GetScalarComponentAsFloat(x,y,0,0);
		    		float fetalRxyVal = SurfFetalTransRxyImage->GetScalarComponentAsFloat(x,y,0,1);
		    		float transRxyVal = SurfFetalTransRxyImage->GetScalarComponentAsFloat(x,y,0,2);

		    		if(x < ex[0] || x > ex[1] || y < ex[2] || y > ex[3] ){
		    		cout << "j of " << s-1 << ": " << j << " x: " << x << ", y: " << y << ", z: " << z << endl;
		    		cout << "surf: " << surfRxyVal << " | " << "fetal: " << fetalRxyVal << " | " << transRxyVal << endl;
		    		}


		    		SurfFetalTransRxyImage->SetScalarComponentFromFloat(x,y,0,0,surfRxyVal+tempSFT[8*j+3]);
		    		SurfFetalTransRxyImage->SetScalarComponentFromFloat(x,y,0,1,fetalRxyVal+tempSFT[8*j+4]);
		    		SurfFetalTransRxyImage->SetScalarComponentFromFloat(x,y,0,2,transRxyVal+tempSFT[8*j+5]);


		    		int alphaRaIdx = int((tempSFT[8*j+6]*180/M_PI)/DetGridDims[5]); // angle from tempSFT in abs(rad)
		    		int radiusIdx = int(tempSFT[8*j+7]/DetGridDims[3]);
		    	//	radiusIdx = int(sqrt(pow(x,2)+pow(y,2)));
		    		if(verbose){
		    			cout << "RadiusIdx: idx(temp):  " <<  tempSFT[8*j+7] << " |  idx: " << radiusIdx << endl;

		    		}

		    	//	cout << "RTdxy" << endl;
		    		int ext[6];
		    		RTdxyImage->GetExtent(ext);


		    		if(radiusIdx>ext[1]){
		    			radiusIdx = ext[1];
		    		}else if(radiusIdx<ext[0]){
		    			radiusIdx = 0;
		    		}

		    		if(alphaRaIdx>ext[3]){
		    			alphaRaIdx = ext[3];
		    		}else if(alphaRaIdx<0){
		    			alphaRaIdx = 0;
		    		}
		    		if(radiusIdx < ext[0] || radiusIdx > ext[1] || alphaRaIdx < ext[2] || alphaRaIdx > ext[3] ){
		    		cout << "RadiusIdx: " << radiusIdx << ", alphaIdx: " << alphaRaIdx << endl;
		    		}
		    		float RdxyVal = RTdxyImage->GetScalarComponentAsFloat(radiusIdx,alphaRaIdx,0,0);
		    		float FdxyVal = RTdxyImage->GetScalarComponentAsFloat(radiusIdx,alphaRaIdx,0,1);
		    		float TdxyVal = RTdxyImage->GetScalarComponentAsFloat(radiusIdx,alphaRaIdx,0,2);

		    		RTdxyImage->SetScalarComponentFromFloat(radiusIdx,alphaRaIdx,0,0,tempSFT[8*j+3]+RdxyVal);
		    		RTdxyImage->SetScalarComponentFromFloat(radiusIdx,alphaRaIdx,0,1,tempSFT[8*j+4]+FdxyVal);
		    		RTdxyImage->SetScalarComponentFromFloat(radiusIdx,alphaRaIdx,0,2,tempSFT[8*j+5]+TdxyVal);





		    	}
		    	if(verbose){
		    		cout << "Size RecordAxzAyz: " << tempAxyz.size() << endl;
		    		cout << "Size RecordAxz: " << tempAxyz.size() << endl;
		    	 //		cout << "Size RAxzCoords: " << tempAxyz.size() << endl;
		    	}

		    	 cl::copy(cmd_queue,d_RecordAxzAyz,tempAxyz.begin(),tempAxyz.end());

		    	 int tas = tempAxyz.size();
		    	 int ext[6];

		    	 ArzImage->GetExtent(ext);

//cout << "TempAyz-Siz:" << tas << endl;

		    	 for(int i=0;i<tas/3;i++){

		    		 //int idx = (int) ((fabs(coords.z)+0.5f)*dimsXYZ.y+coords.y+fabs(coords.z));


		    		if(tempAxyz[3*i]!=0){
		    			//cout << "!=0" << endl;
		    			int cz = int(i/(2*dimR+1));
		    			int Rmax = (dimR-1);
		    			int radius = i-cz*(2*dimR+1)-dimR;



		    			float a = tempAxyz[3*i+1];///ScaleArz2D(abs(radius)); //absorbedWeight


		    			float mua = tempAxyz[3*i+2];
		    			mua = (mua==0) ? 1:mua;



		    			float ArzVal = ArzImage->GetScalarComponentAsFloat(radius,-cz,0,0);
		    			float PhirzVal = ArzImage->GetScalarComponentAsFloat(radius,-cz,0,1);

		    		// cz
		    		// |__ radius

		    	/*	if(radius<ext[0] || radius>ext[1] || cz < ext[2] || cz > ext[3] || a!=0){
		    			cout << "i: " << i << ", cz:" << cz << ", Rmax: " << Rmax << ", dimR: " << dimR << ", radius: " << radius;
		    			cout << " | a: " << a << " , mua: " << mua << endl;


		    		}
		    		*/

		    			ArzImage->SetScalarComponentFromFloat(radius,-cz,0,0,a + ArzVal);
		    			ArzImage->SetScalarComponentFromFloat(radius,-cz,0,1,PhirzVal+a/mua);

		    		}


		    	 }



		    	 //	if(radius>ext[1]){
		    	 //		radius = ext[1];
		    	 //	}else if(radius<ext[0]){
		    	 //		radius = 0;
		    	 //	}

		    	 //	if(cz>ext[3]){
		    	 //		cz = ext[3];
		    	 //	}else if(cz<=0){
		    	 //		cz = 0;
		    	 //	}








		    	if(showPercent){

		    		PrintPercentDone(double((i+1)*100)/numLoops);


		    	}



		    }







		    eTest = Adder(argTest,numTest,d_A,d_B,d_C);

		    eTest.wait();

		    CopyTestDataToHost(d_C,C);

		    eTest.getProfilingInfo(CL_PROFILING_COMMAND_START,&startTimeTest);
		    eTest.getProfilingInfo(CL_PROFILING_COMMAND_END,&endTimeTest);



		   // int wg[3];
	//	    int kcwgs = <CL_KERNEL_COMPILE_WORK_GROUP_SIZE>();


   	//	 cout << "CL_KERNEL_COMPILE_WORK_GROUP_SIZE: " <<  kcwgs << endl;


		   // CopyDataToHost(d_Fetalxy,FetalRxy);
		  ///  CopyDataToHost(d_TransTxy,TransTxy);

		//    CopyTestDataToHost(d_C,C);
		//    CopyDataToHost(d_RecordAxz, RecordAxz);
		//    CopyDataToHost(d_RecordAyz,RecordAyz);
		    // CopyDataToHost(d_WeightLog,WeightLog);


		   // afterCopyToHost = clock();

		  //  printVector("Randoms:",rand3Arr);
	return SUCCESS;

}
void CalculatePropagation::CopyDataToHost(cl::Buffer &buf, std::vector<float> &vec){


	cl::copy(cmd_queue,buf,vec.begin(), vec.end());

}

void CalculatePropagation::CopyTestDataToHost(cl::Buffer &buf, std::vector<int> &vec){


	cl::copy(testcmd_queue,buf,vec.begin(), vec.end());

}

int CalculatePropagation::convertToString(const char *fname, std::string& s)
{
    size_t size;
    char*  str;

    // create a file stream object by filename
    std::fstream f(fname, (std::fstream::in | std::fstream::binary));


    if(!f.is_open())
    {
     	return 1;
    }
    else
    {
        size_t fileSize;
        f.seekg(0, std::fstream::end);
        size = fileSize = (size_t)f.tellg();
        f.seekg(0, std::fstream::beg);

        str = new char[size+1];
        if(!str)
        {
            f.close();
            return 1;
        }

        f.read(str, fileSize);
        f.close();
        str[size] = '\0';

        s = str;
        delete[] str;
        return 0;
    }
}

/*
int CalculatePropagation::fillRandom(std::vector<float> &vec,const int width,const int height, const float rangeMin, const float rangeMax,unsigned int seed){

	if(vec.empty())
    {
        std::cout << "Cannot fill vector." << std::endl;
        return FAILURE;
    }

    // set seed
    if(!seed)
        seed = (unsigned int)time(NULL);

    srand(seed);

    // set the range
    double range = double(rangeMax - rangeMin) + 1.0;

  //  random initialisation of input

    for(int i = 0; i < height; i++)
        for(int j = 0; j < width; j++)
        {
            int index = i*width + j;
            vec[index] = rangeMin + float(range*rand()/(RAND_MAX + 1.0));
        }

    return SUCCESS;
}

*/

void CalculatePropagation::printVector(std::string header, const std::vector<float> vec){
    std::cout<< header<<"\n";

    // print all the elements of the data
    for(std::vector<float>::size_type ix = 0; ix != vec.size(); ++ix)
    {
        std::cout<<vec[ix]<<" ";
    }
    std::cout<<"\n";
}
void CalculatePropagation::printVector(std::string header, const std::vector<int> vec){
    std::cout<<header<<"\n";

    // print all the elements of the data
    for(std::vector<int>::size_type ix = 0; ix != vec.size(); ++ix)
    {
        std::cout<<vec[ix]<<" ";
    }
    std::cout<<"\n";
}
void CalculatePropagation::printVector(std::string header, const std::vector<long> vec){
    std::cout<<header<<"\n";

    // print all the elements of the data
    for(std::vector<long>::size_type ix = 0; ix != vec.size(); ++ix)
    {
        std::cout<<vec[ix]<<" ";
    }
    std::cout<<"\n";
}
int CalculatePropagation::setup(){
	/* setup: (fft: L.92 - 199)
		 * 	device: cpu | gpu
		 * 	SetContextPlatform
		 * 	clCommandQueue
		 * 	CreateBuffer
		 * 	create CL program from source
		 * 	build program
		 * 	create kernel from program
	*/

	// Init
	InitLogVectors();

	// Get all necessary data for computation


		createContext();

		//PSNR
		PSNR *rand = new PSNR;

		rand->InitRand3();
		rand->GetRan3Array(rand3Arr);
		//rand->GetRan3(rand3Arr,wgSize);
		rand->GetInextInextp(inextInextp);//inextInxtpPoolID
		rand->GetMBMZ(MBMZ);


		//printVector("Randoms:",rand3Arr);
		//printVector("inextInextp:",inextInextp);
		//printVector("MBMZ:",MBMZ);
		// cout << "PoolSize: " << poolSize << endl;


	//	CreateBufferObjects();
		CreateBufferObjects(context);
		CreateTestBuffer();
		 //create queue to which we will push commands for the device.
		CreateCommandQueue();


		BuildKernel();

		return SUCCESS;

}

int CalculatePropagation::run(){
// in loop

	/* run: (204 - 332)
	 * 	copy variables to device
	 * 	flush commandqueue
	 * 	wait for data transfer to device
	 * 	set workgroup sizes
	 * 	kernel args
	 * 	global and local threads
	 * 	wait for events
	 * 	read output  data:
	 * 		wait for data transfer to host
	 */



		ExecKernel();


	return SUCCESS;
}

int CalculatePropagation::cleanup(){
	/*
	 * cleanup: (609 - 650)
	 * 	release memory objects
	 * 	release Kernel
	 * 	release program
	 * 	release commandQueue
	 * 	release context
	 */


	d_vStart.~Memory();
	d_vStartOcc.~Memory();
	//d_pStart.~Memory();

	d_RandNmbersLog.~Memory();
	d_DimsXYZ.~Memory();
	d_SurfRxy.~Memory();
	//d_RecordAxz.~Memory();
	//d_RecordAyz.~Memory();
	d_RecordAxzAyz.~Memory();
	d_TransTxy.~Memory();
	d_Fetalxy.~Memory();
	d_SurfFetalTransRxy.~Memory();
	d_eP0XYZ.~Memory();
	d_DetGridDims.~Memory();
	d_lockedRecord.~Memory();
	d_triangles.~Memory();
	d_normals.~Memory();
	//d_n.~Memory();
	d_triN.~Memory();
//	d_mua.~Memory();
//	d_mus.~Memory();
//	d_g.~Memory();
	d_nmuasG.~Memory();
	d_randPool.~Memory();
	d_rand3Array.~Memory();
	d_MBMZ.~Memory();
	d_inextInextp.~Memory();


	//Test

	d_A.~Memory();
	d_B.~Memory();
	d_C.~Memory();


	kernel.~Kernel();
	program.~Wrapper();
//	cmd_queue.~Wrapper();

//	context.~Context();


	return SUCCESS;
}


int CalculatePropagation::Init(){
	// Setup
	if(setup() != SUCCESS){
		return FAILURE;
	}

}

int CalculatePropagation::Propagation(){


	        // Run
	        if(run() != SUCCESS)
	        {
	            return FAILURE;
	        }
	        // Cleanup
	        if(cleanup() != SUCCESS)
	        {
	            return FAILURE;
	        }



	    return SUCCESS;
}

void CalculatePropagation::PrintSimulationParams(){

	cout << "########### Simulation Parameter ###########" << endl;
	cout << "Device Name: " << selectedDevice.getInfo<CL_DEVICE_NAME>() << endl;
	//cout << "Max. work group size: " << wgSize << endl;
	//cout << "Local work group size: " << lWGSize << endl;
	//cout << "Global work size: " << gWSize[0] << endl;
//	cout << "Number of Workgroups: " << numberPhotons/lWGSize << endl;
	cout << "Number of Photons: " << numberPhotons << endl;
	cout << "Number of Photons/Loop: " << photonsPerLoop << endl;








}

void CalculatePropagation::PrintPercentDone(double p){

//	cout << "Calculate Propagation...." << endl;
	if(p<100){
		cout << "\r\r| " << p << "% " << std::flush;
	}else{
		cout << "\r\r| 100% \t" << std::endl;
	}



}

void CalculatePropagation::PrintTimeSpent(){


    double timeTestInNS = double(endTimeTest - startTimeTest); // calc time
    cout << "############# Calculation Time #############" << endl;
    cout << "Time in ms: " << double(neededTime)/1000000 << " \t| Time in s: " <<double(neededTime)/1000000000 << endl;
    cout << "Time in ms: " << timeTestInNS/1000000 << " \t| Time in s: " <<timeTestInNS/1000000000 << " (numTest: " << numTest << ")" <<  endl;
}

void CalculatePropagation::PrintMediumIDs(){
	printVector("MediumIDs: ", this->triN);

}
void CalculatePropagation::PrintRefractionIndices(){
	printVector("N, mua, mus, g : ", this->nmuasG);
}
