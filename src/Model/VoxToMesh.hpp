/*
 * VoxToMesh.hpp
 *
 *  Created on: Mar 9, 2015
 *      Author: matthias
 */

#ifndef VOXTOMESH_HPP_
#define VOXTOMESH_HPP_

#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkStripper.h"
#include "vtkImageToPolyDataFilter.h"

#include "vtkImageData.h"
#include "vtkVolumeProperty.h"
#include "vtkPiecewiseFunction.h"

#include "vtkSmartPointer.h"

#include "vtkDecimatePro.h"
#include "vtkQuadricDecimation.h"
#include "vtkDelaunay3D.h"
#include "vtkDataObject.h"

#include "vtkActor.h"
#include "vtkOutlineFilter.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"

#include "vtkCellCenters.h"
#include "vtkGlyph3D.h"
#include "vtkArrowSource.h"

#include "vtkDoubleArray.h"

#include "vtkSmartVolumeMapper.h"

#include "vtkVolume.h"

class VoxToMesh {


	vtkContourFilter *contour;
	vtkPolyDataNormals *normals;
	vtkSmartPointer<vtkImageData> VoxImage;
	vtkSmartPointer<vtkImageData> IDVoxelImage;
	vtkPolyData *PolyImage;
	vtkStripper *stripper;

	vtkPolyDataMapper *mapper;
	vtkActor *actor;
	double density_double;
	vtkDoubleArray *density_doubleArray;
	double reductionfactor;
	vtkDoubleArray *reductionfactor_doubleArray;
	double featureAngle;
	vtkOutlineFilter *outlineData;
	vtkPolyDataMapper *mapOutline;
	vtkActor *outline;
	int iter;
	double converge;
	double relax;
	int densID;
	double opacity;
	double arrowScale;

	// vtkDoubleArray* normals_array;
	vtkCellCenters *Cellcenter;
	//vtkDoubleArray* normalcenters_array;

	vtkVolume* VoxVolume;
	vtkSmartPointer<vtkActor> actorArrow;

	vtkDoubleArray* mediumIDs;




public:



	void convertToMesh();

	// Setter funcitons
	void SetDensityValue(double);
	void SetDensityValue(vtkDoubleArray*);
	void SetDensityValue(std::vector<float> &densities);
	void SetVoxImage(vtkImageData*);
	void SetTargetReduction(double);
	void SetTargetReduction(vtkDoubleArray*);
	void SetFeatureAngle(double);
	void SetIterations(int);
	void SetConvergenceFactor(double);
	void SetRelaxationFactor(double);
	void SetVolumeOpacity(double);
	void SetArrowScale(double);
	// Getter functions
	vtkDataArray* GetNormals();
	vtkPolyData* GetNormalsOutput();
	vtkPolyData* GetMeshImage();
	vtkActor* GetActor();
	// vtkPolyDataMapper GetMapper();
	vtkSmartPointer<vtkImageData> GetVoxImage();
	vtkSmartPointer<vtkImageData> GetIDVoxelImage();

	int GetNumberOfContours();
	void GetContourValues(double* values);
	double* GetContourValues();


	void SetDensID(int id);

	vtkSmartPointer<vtkActor> GetNormalArrowsActor();
	vtkVolume* GetVoxelVolume();
	vtkDoubleArray* GetNormalCenters();
	void PrintNormalsCenters();

	void GenMediumIDs(int mode);
	vtkCellCenters* GetCellcenters();
	void SetCellcenters();

	void GetMediumIDs(double**);
	vtkDoubleArray* GetMediumIDs();




	VoxToMesh();
	virtual ~VoxToMesh();




private:

	void GenMediumIDs0();
	void GenMediumIDs1();
	void AddPointsToPoints(vtkPoints* inputToAdd, vtkPoints* output);
	void BarycentricToCartesian(vtkPoints*, vtkDoubleArray*, vtkPoints*);
	void MediumPoints(vtkPoints* cartesian, double* normalVector, vtkPoints* mediumPoints);
	vtkPoints* MediumPoints(vtkPoints* cartesian, double* normalVector);
	void Cross(double* a, double* b, double* c);
	double Dot(double* a, double* b);
	void Vec(double* a, double* b, double* vec);
	void Normalize(double* v);
	void PointInPlaneCoords(double* p0, double* point, double* ex, double* ey,double* planeCoords);
	void PointInPlaneCoords(double* p0, double* point, double* ex, double* ey, int i, vtkDoubleArray* planeCoords);
	void PlaneCoordsTo3D(double x, double y, double* p0, double* ex, double* ey, double* coords);
	vtkDoubleArray* GetDensityDoubleArray();

	int GetDensID();




};

#endif /* VOXTOMESH_HPP_ */
