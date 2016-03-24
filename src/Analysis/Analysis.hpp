/*
 * Analysis.hpp
 *
 *  Created on: Jul 5, 2015
 *      Author: matthias
 */

#ifndef ANALYSIS_HPP_
#define ANALYSIS_HPP_

#include <vector>
#include <string>

#include "vtkPolyData.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include <vtkImageData.h>
#include <vtkSmartPointer.h>


class Analysis {

	bool verbose;

	vtkPolyData* mesh;
	vtkDoubleArray* normals;

	std::vector<float> histogram_arr;
	std::vector<float> random_arr;
	std::vector<double> RdTd;

	//std::vector<float> absorbtionLog;
	// std::vector<float> positions;

	std::vector<float> absorbYZLog;
	std::vector<float> absorbXZLog;
	std::vector<float> surfRxyLog;
	std::vector<float> transTxyLog;
	std::vector<float> fetalRxyLog;

	std::vector<float> SurfFetalTransRxy;

	float fetalRxySum, surfRxySum, transTxySum;




	std::vector<int> surfaceIDs; // IDs der Cells | size: numTris
	std::vector<std::string> surfaceIDnames; // ID aus object file | size: numTris
	std::vector<std::string> surfaceIDnamesConfig; // ID names aus config | size: z.b. 2
	std::vector<double> surfIDSums; // z.B. {nid,Rd, Td};
	std::vector<double> surfIDSumsNum; // z.B. {nid,Rd, Td};


	float chiSquared;
	int numSamples;
	int numBins;
	int numberPhotons;









public:

	void CalcHistogram(std::vector<float> &randomArrIn, int bins, std::vector<float> &histArrOut);
	void CalcHistogram(std::vector<float> &randomArrIn, int bins);
	void CalcHistogram(int bins);

	void CalcChiSquared(std::vector<float> &histArrIn, int numSamples);
	void CalcChiSquared();


	void SetNumberOfPhotons(int num);
	void SetVerbosity(bool verb);
	void SetVerbosityTrue();
	void SetVerbosityFalse();

	void SetRxyLogs(std::vector<float> &surfRxy, std::vector<float> &fetalRxy, std::vector<float> &transRxy);
	void SetRxyLogs(std::vector<float> &surffetaltransRxy);
	void SetRTdxyLogs(vtkSmartPointer<vtkImageData> image);
	void SetRandomArray(std::vector<float> &randomArr);
	void SetHistogramArray(std::vector<float> &histArrIn);
	void SetAbsorbtionLogs(std::vector<float> &absorbtionXZLogIn,std::vector<float> &absorbtionYZLogIn);
	void SetTransmittanceLog(std::vector<float> &transmittanceLogIn);
	void SetMesh(vtkPolyData* meshIn);
	void SetNormals(vtkDataArray* normalsIn);

	void SetSurfaceIDs(std::vector<int> &surfIDsIn);
	void SetSurfaceIDNames(std::vector<std::string> &surfIDNamesIn);
	void SetSurfaceIDNamesConfig(std::vector<std::string> &surfIDNamesIn);

	void GetHistogram(std::vector<float> &histArrOut);
	std::vector<float> GetHistogram();
	float GetChiSquared(std::vector<float> &histArrIn, int numSamples);
	float GetChiSquared();


	void CalculateSurfaceSums();
	void GenerateSurfaceImage();




	double CalculateRd(std::vector<float> &vec, int numPhotons, float top);
	double CalculateTd(std::vector<float> &vec, int numPhotons, float bottom);
	void CalculateRTd(std::vector<float> &vec,std::vector<float> &pz, float top, float bottom, std::vector<double> &out);
	void CalculateRTd(std::vector<float> &vec,std::vector<float> &pz, float top, float bottom);

	void PrintChiSquared();
	void PrintRTd();
	void PrintSurfaceSums();
	void PrintAnalysis();



	Analysis();
	~Analysis();




private:

	bool IsInCell(double position[3], double cellpoint[3],double normal[3], double maxDistance);
	int GetSurfaceID(std::string surfaceNameIn);

};

#endif /* ANALYSIS_HPP_ */
