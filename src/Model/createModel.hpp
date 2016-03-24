/*
 * createModel.hpp
 *
 *  Created on: Mar 3, 2015
 *      Author: matthias
 */

#ifndef CREATEMODEL_HPP_
#define CREATEMODEL_HPP_
#include  "../IO/ReadXML.hpp"

#include "vtkCommonCoreModule.h"
#include "vtkSmartPointer.h"
#include "vtkCellArray.h"
#include "vtkArray.h"
#include "vtkTypedArray.h"
#include "vtkIntArray.h"
#include "vtkDoubleArray.h"
#include "vtkDenseArray.h"
#include "vtkObject.h"
#include "vtkMapper.h"
#include <vector>

#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"

#include "vtkImageImport.h"
#include "vtkImageData.h"
// #include "vtkImageResample.h"
#include "vtkSmartPointer.h"
#include "vtkStripper.h"
#include "vtkActor.h"


#include "vtkCellCenters.h"
#include "vtkGlyph3D.h"

// #include "vtkAbstractVolumeMapper.h"
// #include "vtkVolume.h"
// #include "vtkVolumeProperty.h"
// #include "vtkFixedPointVolumeRayCastMapper.h"




class createModel {


	vtkIntArray* dimensions_array;
	vtkDoubleArray* densities_array;
	vtkIntArray* layerdepths_array;
	vtkDoubleArray* rotation_array;
	int modelID_int;// 0: Layer, 1: Spheres, 2: 3D (DICOM?)
	int modeID_int; // 0: voxel, 1: mesh
	double arrowScale;
	vtkSmartPointer<vtkImageData> image;
//	char* filename;
	vtkSmartPointer<vtkPolyData> PolyImage;
	vtkDoubleArray* mediumIDs;
	vtkPolyDataNormals *normals;
	vtkCellCenters *Cellcenter;
	vtkSmartPointer<vtkActor> actorArrow;
	vtkStripper *stripper;
	vtkPolyDataMapper *mapper;
	vtkActor *actor;



public:

	// vtkTypeMacro(createModel, vtkObject);
 	static createModel *New();
 	void GenMediumIDs();
	// Getter functions
	vtkDataArray* GetNormals();
	vtkPolyData* GetNormalsOutput();
	vtkPolyData* GetMeshImage();
	vtkActor* GetActor();
	vtkSmartPointer<vtkActor> GetNormalArrowsActor();
	vtkCellCenters* GetCellcenters();
	void SetCellcenters();
	void SetArrowScale(double);
	//ReadXML *reader;

 	//void setFilename(char*);
	void setModel(int);
	//void setInputVarsLayer();
	void setInputVarsLayer(vtkDoubleArray*,vtkIntArray*, vtkDoubleArray*);
	void setInputVarsLayer(std::vector<float> &densities,std::vector<int> &layerdepths, std::vector<float> &rotations);
	vtkDoubleArray* GetMediumIDs();
	void setSize(int,int,int);
	void setSize(vtkIntArray*);
	void setMode(int);
	void Update();
//	vtkDoubleArray* GetDensities();
//	double GetDensities(int);
	vtkSmartPointer<vtkImageData> GetOutput();
	vtkSmartPointer<vtkPolyData> GetOutputPolyImage();
	 void PrintVariable();
	 void PrintNormalsCenters();

	createModel();
	~createModel(){};


// protected:


private:

	int getLayer(int, int, int);
	unsigned int getLayerDensity(int, int, int);

	void setDensities(vtkDoubleArray*, unsigned);
	void setLayerDepth(vtkIntArray*, unsigned);

	void setRotation(vtkDoubleArray*);
	void setRotation(double, double, double);

	// void rotateImage(); // evtl. nicht nötig

	// void GetVtkImage(int*);
	void CreateLayerModel(int,int,int);
	void CreateLayerModel();
	void CreateSphereModel(int,int,int); // not yet implemented
	void Create3DModel(int,int,int);  // not yet implemented
	void Rotate(std::vector<double> pIn, std::vector<double> &pOut, double* rotPoint);



	// createModel(const createModel&);
	// void operator=(const createModel&);

};


#endif /* CREATEMODEL_HPP_ */
