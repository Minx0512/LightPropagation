/*
 * ImportExport.hpp
 *
 *  Created on: Mar 22, 2015
 *      Author: matthias
 */

#ifndef IMPORTEXPORT_HPP_
#define IMPORTEXPORT_HPP_


#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include "vtkPolyData.h"
#include "vtkXMLPolyDataReader.h"
#include "vtkCleanPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataNormals.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkStripper.h"
#include "vtkActor.h"
#include "vtkCellCenters.h"
#include "vtkGlyph3D.h"

class ImportExport {


	vtkPolyData* importmesh;
	//vtkSmartPointer<vtkXMLPolyDataReader> reader;
	//vtkSmartPointer<vtkXMLPolyDataWriter> writer;
	//std::vector<float> absorbtion;
	//std::vector<float> scatter;
	//std::vector<float> refraction;
	//std::vector<float> anisotropy;
	std::vector<float> normal_vector;


	std::vector<float> refractionInd;
	std::vector<float> absorbtionCoeffs;
	std::vector<float> anisotropies;
	std::vector<float> scatterCoeffs;
	std::vector<float> densities;
	std::vector<int> layerdepths;
	std::vector<float> rotations;
	std::vector<float> view_array;
	std::vector<std::string> tissueNames;
	// double dims[3];
	std::vector<float> dims;
	std::vector<float> detGridDims;
	std::vector<int> surfaceIDs;
	std::vector<std::string> surfaceIDnames;
	std::vector<std::string> surfaceIDnamesConfig;


	std::vector<int> picDims;

	bool verbose;

	std::vector<float> angles;
	std::vector<float> relRadiantIntensity;

	std::vector<float> startTrajectoryX;
	std::vector<float> startTrajectoryY;
	std::vector<float> startTrajectoryZ;

	std::vector<float> startPointX;
	std::vector<float> startPointY;
	std::vector<float> startPointZ;

	std::vector<float> emitterHeading;


	float ppu;
	std::string configfile;
	std::string modelfile;
	std::string ModelName;
	std::string emitterFileName;
	std::string saveDir;
	std::string unit;

	std::string surfaceFilename;
	std::string transmissionFilename;
	std::string fetalFilename;
	std::string recordAxzFilename;
	std::string recordAyzFilename;

	std::string scaleString;
	std::string coordsysString;
	std::string colormapString;

	double arrowScale;
	int numPhotons;
	int numPhotonsinteractions;
	double chance;
	double w_thresh;
	bool saveFile;
	std::vector<float> emitterPoint;
	float emitterDistance;
	int fetalID;
	std::vector<int> tissueOuterIDs;



	vtkSmartPointer<vtkPoints> points;
	vtkSmartPointer<vtkCellArray> Vertices;
	vtkDoubleArray* mediumIDs;
	vtkPolyDataNormals *normals;

	vtkCellCenters *Cellcenter;
	vtkSmartPointer<vtkActor> actorArrow;
	vtkStripper *stripper;
	vtkPolyDataMapper *mapper;
	vtkActor *actor;


public:

	void SetFilename(std::string);
	void SetConfig(std::string filename);
	void ParseFile();
	void SetVerbosity(bool verb);
	void SetVerbosityTrue();
	void SetVerbosityFalse();


	std::string GetSurfaceFilename();
	std::string GetFetalFilename();
	std::string GetTransmissionFilename();
	std::string GetRecordAxzFilename();
	std::string GetRecordAyzFilename();

	std::string GetScaleString();
	std::string GetCoordSystemString();
	std::string GetColormapString();

	//void GenMediumIDs();
	//vtkDoubleArray* GetMediumIDs();
		// Getter functions
	vtkDataArray* GetNormals();
	vtkPolyData* GetNormalsOutput();
	vtkPolyData* GetMeshImage();

	vtkActor* GetActor();
	vtkSmartPointer<vtkActor> GetNormalArrowsActor();
	vtkCellCenters* GetCellcenters();

	void SetCellcenters();
	void SetArrowScale(double);

	void ExportToVTP(vtkPolyData*, std::string);

	void ImportModel(std::string fname);
	void ImportXMLModel(std::string fname);
	void ImportObjModel(std::string fname);
	void ImportEmitterModel(std::string file_name);
	void ExportModel(std::string fname);

	// Get
	void GetAbsorbtionCoefficients(std::vector<float> &absorbtions);
	void GetScatterCoefficients(std::vector<float> &scatter);
	void GetAnisotropies(std::vector<float> &anisotropies);
	void GetRefractionIndices(std::vector<float> &refractions);
	float GetPixelPerUnitFactor();
	vtkDoubleArray* GetMediumIDs();
	void GetMediumIDs(vtkSmartPointer<vtkDoubleArray> &medIDsOut);
	vtkSmartPointer<vtkPoints> GetPoints();
	void GetPoints(vtkSmartPointer<vtkPoints> &points);
	vtkSmartPointer<vtkCellArray> GetVertices();
	void GetVertices(vtkSmartPointer<vtkCellArray> &verts);
	void GetNormalVectors(std::vector<float> &normalVector);
	std::string GetUnit();
	std::string GetModelName();
	vtkIntArray* GetDimensions();
	void GetDimensions(std::vector<float> &dims);
	void GetDetectionGrid(std::vector<float> &detGrid);
	double GetViewAzimuth();
	double GetViewElevation();
	double GetViewRoll();
	int GetNumberOfPhotons();
	int GetNumberOfPhotonInteractions();
	bool GetSaveFile();
	std::string GetSaveDir();
	double GetArrowScaleFactor();
	double GetChance();
	double GetWeightThreshold();
	void GetSurfaceNamesVector(std::vector<std::string> &surfNamesVectorOut);
	void GetSurfaceConfigNames(std::vector<std::string> &surfNamesVectorOut);
	void GetSurfaceIDVector(std::vector<int> &surfIDVectorOut);

	void GetPictureDimensions(std::vector<int> &dims);

	float GetEmitterdistanceToSurface();
	void GetEmitterAngles(std::vector<float> &angles);
	void GetEmitterRelativeRadiantIntensities(std::vector<float> &relRadInt);
	void GetEmitterPointInPlane(std::vector<float> &emitterPiP);
	void GetEmitterHeading(std::vector<float> &heading);

	int GetFetalID();
	//Set
	void SetAbsorbtionCoefficients(std::vector<float> &absorbtions);
	void SetScatterCoefficients(std::vector<float> &scatter);
	void SetAnisotropies(std::vector<float> &anisotropies);
	void SetRefractionIndices(std::vector<float> &refractions);
	void SetPixelPerUnitFactor(float ppu);
	void SetMediumIDs(vtkDoubleArray*);
	void SetPoints(vtkSmartPointer<vtkPoints> points);
	void SetVertices(vtkSmartPointer<vtkCellArray> verts);
	void SetNormalVectors(std::vector<float> &normalVector);
	void SetModelName(char* modelname);
	void SetUnit(char* unit);

	void GenMediumIDsFromObjModel(std::string file_name);
	void PrintNormalsCenters();

	void PrintCells();
	void PrintNormals();
	void PrintMediumIDs();
	void PrintImport();

	ImportExport();
	virtual ~ImportExport();

private:

	void ExtractMediumsList(vtkXMLDataElement *mediumsList);
	void ExtractName(vtkXMLDataElement *name);
	void ExtractUnit(vtkXMLDataElement *unit);
	void ExtractMesh(vtkXMLDataElement *mesh);
	void ExtractPointsList(vtkXMLDataElement *pointsList);
	void ExtractCellsList(vtkXMLDataElement *cellsList);

	void ExtractPixelPerUnit(vtkXMLDataElement *ppu);
	void Clean();

	//int sumVector(std::vector<int> &vec, int idx);



};

#endif /* IMPORTEXPORT_HPP_ */
