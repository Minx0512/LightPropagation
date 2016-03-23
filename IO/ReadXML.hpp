/*
 * ReadXML.hpp
 *
 *  Created on: Mar 7, 2015
 *      Author: matthias
 */

#ifndef READXML_HPP_
#define READXML_HPP_


#include "vtkXMLReader.h"
#include "vtkXMLTreeReader.h"
#include "vtkXMLDataParser.h"
#include "vtkTree.h"
#include "vtkXMLDataElement.h"
#include "vtkXMLDataParser.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"
#include "vtkDoubleArray.h"
#include "vtkObject.h"
#include <vector>


class ReadXML {


	vtkXMLDataParser *configParser;
	vtkXMLDataParser *ModelXMLParser;


	// Model properties
	vtkDoubleArray* dimensions_array;
	vtkDoubleArray* rotation_array;
	vtkDoubleArray* view_array;
	vtkDoubleArray* reductionfactors_array;


	// Tissue properties
	vtkStringArray* tissueNames_array;
	vtkIntArray* layerdepths_array; // Layerdepth in px
	vtkDoubleArray* thickness_array; // Layerthickness in cm
	vtkDoubleArray* densities_array; // densities in grey scale [0,255]
	vtkDoubleArray* anisotropies_array;
	vtkDoubleArray* absorbtion_array;
	vtkDoubleArray* scatter_array;
	vtkDoubleArray* refractions_array;


	std::vector<int> surfaceIDs;
	std::vector<std::string> surfaceIDnames;
	std::vector<std::string> surfaceIDnamesConfig;
	int modelID_int;
	int modeID_int;
	float ppunit;
	std::string filename;
	std::string modelFileName;
	std::string saveDir;
	std::string emitterFileName;

	std::string surfaceFilename;
	std::string transmissionFilename;
	std::string fetalFilename;
	std::string recordAxzFilename;
	std::string recordAyzFilename;


	std::string scaleString;
	std::string coordsysString;
	std::string colormapString;




	int numTissues;
	int currTissueID;
	int fetalID;
	int saveFile;
	int numPhotons;
	int numPhotonsinteractions;
	double arrScaleFactor;
	double chance;
	double w_thresh;

	float emitterDistance;
	std::vector<float> emitterPoint;
	std::vector<float> detGridDims;
	std::vector<float> picDims;
	std::vector<float> emitterHeading;

	bool verbose;



public:


	//static ReadXML *New();

		void ParseModelXMLFile(std::string fname);
		void ParseConfig(std::string fname);
		void SetFilename(std::string);
		void SetVerbosityTrue();
		void SetVerbosityFalse();
		// Emitter
		std::string GetEmitterFile();

		// Get Model properties
		std::string GetModelFilename();
		int GetModelID();
		int GetModeID();
		float GetPixelPerUnitFactor();
		vtkDoubleArray* GetReductionFactors();
		double GetReductionFactors(int);
		vtkDoubleArray* GetDimensions();
		void GetDimensions(std::vector<float> &dims);
		void GetDetectionGrid(std::vector<float> &detGrid);
		vtkDoubleArray* GetRotations();
		void GetRotations(std::vector<float> &rot);

		// Get Tissue properties
		std::vector<std::string> GetTissueNames();
		void GetTissueNames(std::vector<std::string> &lNames);
		//vtkIntArray* GetLayerdepth();
		std::vector<int> GetLayerdepth();
		void GetLayerdepth(std::vector<int> &ldpth);
		std::vector<float> GetThicknesses();
		void GetThicknesses(std::vector<float> &thicks);
		//vtkDoubleArray* GetDensities();
		//std::vector<float> GetDensities();
		void GetDensities(std::vector<float> &denses);
		double GetDensities(int);
		std::vector<float> GetDensities();

		std::vector<float> GetAnisotropies();
		void GetAnisotropies(std::vector<float> &anisotropies);
		std::vector<float> GetAbsorbtionCoefficients();
		void GetAbsorbtionCoefficients(std::vector<float> &absorbtions);
		std::vector<float> GetScatterCoefficients();
		void GetScatterCoefficients(std::vector<float> &scatter);
		std::vector<float> GetRefractionIndices();
		void GetRefractionIndices(std::vector<float> &refractions);

		void GetEmitterHeading(std::vector<float> &heading);

		double GetViewAzimuth();
		double GetViewElevation();
		double GetViewRoll();

		void GetViewRotations(std::vector<float> &rot);

		int GetNumberOfPhotons();
		int GetNumberOfPhotonInteractions();
		bool GetSaveFile();

		std::string GetSaveDir();

		std::string GetSurfaceFilename();
		std::string GetFetalFilename();
		std::string GetTransmissionFilename();
		std::string GetRecordAxzFilename();
		std::string GetRecordAyzFilename();

		std::string GetScaleString();
		std::string GetCoordSystemString();
		std::string GetColormapString();


		void GetPictureDimensions(std::vector<int> &dims);
		int GetFetalID();

		double GetArrowScaleFactor();
		double GetChance();
		double GetWeightThreshold();

		void GetSurfaceNamesVector(std::vector<std::string> &surfNamesVectorOut);
		void GetSurfaceConfigNames(std::vector<std::string> &surfNamesVectorOut);
		void GetSurfaceIDVector(std::vector<int> &surfIDVectorOut);

		float GetEmitterdistanceToSurface();
		void GetEmitterPointInPlane(std::vector<float> &emitterPiP);

		void PrintArray(std::string header, vtkDoubleArray* arr);
		void PrintArray(std::string header, vtkIntArray* arr);
		void PrintArray(std::string header, vtkStringArray* arr);

		void PrintVars();

		ReadXML();
		virtual ~ReadXML();


private:

		void SetDimensions(double*);
		void SetDensities(double*);
		void SetLayerdepth(int*);
		void SetRotations(double*);
		void SetModelID(int);
		void SetModeID(int);


		void InitTissueArrays(int);
		void ExtractVars(vtkXMLDataElement*);
		void ExtractConfigVars(vtkXMLDataElement*);
		void ExtractReductionfactors(vtkXMLDataElement*);
		void ExtractDimensions(vtkXMLDataElement*);
		void ExtractDensities(vtkXMLDataElement*);
		void ExtractLayerdepths(vtkXMLDataElement*);
		void ExtractRotations(vtkXMLDataElement*);
		void ExtractModelID(vtkXMLDataElement*);
		void ExtractModeID(vtkXMLDataElement*);
		void ExtractPixelPerUnitFactor(vtkXMLDataElement*);
		void ExtractTissueID(vtkXMLDataElement*);
		void ExtractTissueName(vtkXMLDataElement*);
		void ExtractThickness(vtkXMLDataElement*);
		void ExtractRefractionIndices(vtkXMLDataElement*);
		void ExtractAbsorbtionCoefficients(vtkXMLDataElement*);
		void ExtractScatterCoefficients(vtkXMLDataElement*);
		void ExtractAnisotropies(vtkXMLDataElement*);


		void ExtractModelFile(vtkXMLDataElement*);
		void ExtractSaveLogFiles(vtkXMLDataElement*);
		void ExtractSaveDir(vtkXMLDataElement*);
		void ExtractNumberOfPhotonInteractions(vtkXMLDataElement*);
		void ExtractNumberOfPhotons(vtkXMLDataElement*);
		void ExtractInitViewRotation(vtkXMLDataElement*);
		void ExtractArrowScaleFactor(vtkXMLDataElement*);
		void ExtractTopBottom(vtkXMLDataElement*);
		void ExtractSimulationParams(vtkXMLDataElement*);
		void ExtractTissueList(vtkXMLDataElement*);
		void ExtractTissue(vtkXMLDataElement*);
		void ExtractModel(vtkXMLDataElement*);
		void ExtractAnalyse(vtkXMLDataElement*);
		void ExtractAnalyseBoundaries(vtkXMLDataElement*);
		void ExtractAnalyseSurfaceNames(vtkXMLDataElement*);
		void ExtractEmittermodel(vtkXMLDataElement *emittermodel);
		void ExtractEmitterfile(vtkXMLDataElement *emitterfile);

		void ExtractLogs(vtkXMLDataElement *logsroot);
		void ExtractLogsDetectionGrid(vtkXMLDataElement *logsDetGrid);
		void ExtractFilenames(vtkXMLDataElement *filenamesroot);
		void ExtractLogsFileDimensions(vtkXMLDataElement *fileDimensions);



};

#endif /* READXML_HPP_ */
