/*
 * Logs.hpp
 *
 *  Created on: Jul 26, 2015
 *      Author: matthias
 */

#ifndef LOGS_HPP_
#define LOGS_HPP_


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>

#include "../ImageFilter/ImageFilter.hpp"

#define Rxy_Log_SURF 0
#define Rxy_Log_FETAL 1
#define Rxy_Log_TRANS 2

#define Axz_Log 0
#define Ayz_Log 1


class Logs {

	std::vector<float> surfRxyLog;
	std::vector<float> transTxyLog;
	std::vector<float> fetalRxyLog;

	std::vector<float> AxzLog;
	std::vector<float> AyzLog;

	std::vector<int> AxzCoords;
	std::vector<int> AyzCoords;

	std::vector<int> AxyzDims;
	vtkSmartPointer<vtkImageData> imageSFT;
	vtkSmartPointer<vtkImageData> imageRFTd;


	std::string scaleString;
	std::string coordsysString;
	std::string colormapString;

	std::string scaleStringFname;
	std::string coordsysStringFname;
	std::string colormapStringFname;
	float pixelPerUnitFactor;
	int numberPhotons;


public:

	std::string saveDir;
	std::vector<float> dimensions;
	std::vector<float> gridDims; // x,y,z,r,phi,rho



	void SetSaveDir(std::string saveDir);
	void SetDimensions(std::vector<float> &dims);
	void SetGridDimensions(std::vector<float> &gdims);
	void SetPixelPerUnitFactor(float factor);
	void SetNumberOfPhotons(int number);
	void SetRTdxyLogs(vtkSmartPointer<vtkImageData> image);
	void SetSFTLogs(vtkSmartPointer<vtkImageData> image);
	void SetRxyLogs(std::vector<float> &surffetaltransRxy);
	void SetAxzLog(std::vector<float> &Axz);
	void SetAyzLog(std::vector<float> &Ayz);

	void SetAxzCoords(std::vector<int> &Axzcoords);
	void SetAyzCoords(std::vector<int> &Ayzcoords);

	void SetScaleString(std::string);
	void SetCoordSystemString(std::string);
	void SetColormapString(std::string);

	void EvalScaleForFname(int);
	void EvalCoordsysForFname(int);
	void EvalColormapForFname(int);




	template<typename T> void WriteArrInFile(std::string fname, std::vector<T> &vec){
		std::ofstream myfile;
		std::string filename = saveDir;
		filename.append(fname.begin(),fname.end());
		std::cout << "WriteArrInFile Filename: " << fname << "| Number ofElements: " << vec.size() <<  std::endl;
		myfile.open (filename.c_str());

		for(T &arr : vec){
			if(arr>=0){
			myfile << arr << std::endl;
			}
		}
		myfile.close();

	}

	void SaveSurfFetTransArr(std::string fnameSurf,std::string fnameFetal,std::string fnameTrans);

	void SaveSurfaceArr(std::string fname);

	void SaveSurfaceArr(std::string fname, int LOG);

	void SaveSurfaceArr(std::string fname, vtkSmartPointer<vtkImageData> image);

	void SaveAbsorbtionLog(std::string fname,std::vector<float> &absobArr, std::vector<int> &coords);
	void SaveRecordRd(std::string fname,vtkSmartPointer<vtkImageData> image);
	float ScalePerRadius(int ir);
	float ScalePerAngle1D(int ia);
	float ScalePerRadius1D(int ir);
	float Scale2D(int ir, int ia);

	void ScaleRdTt(vtkSmartPointer<vtkImageData> imageRdRA, std::vector<float> &RdR, std::vector<float> &RdA);

	void CalcReflectancePerRadius(std::string fname, int LOG);
	void CalcReflectancePerRadius(std::string fname,vtkSmartPointer<vtkImageData> image);

	void CalcReflectancePerAngle(std::string fname, int LOG);
	void CalcReflectancePerAngle(std::string fname, vtkSmartPointer<vtkImageData> image);

	void CreateSFT_XYPicture(std::string fnameS,vtkSmartPointer<vtkImageData> imageData);
	void CreateRdTt_XYPictureScaled(std::string fnameS,vtkSmartPointer<vtkImageData> imageData);

	void CreateSFT_XYPicture(std::string fnameS,int LOG, int dimR, int dimPhi);

	void CreateAbsorbtionPicture(std::string fname,std::vector<float> &absobArr, std::vector<int> &coords, int dim0, int dim1);
	void CreateAbsorbtionPicture(std::string fname,int Aplane);
	void CreateAbsorbtionPicture(std::string fnameS,vtkSmartPointer<vtkImageData> imageData);


	void FluenceZ(std::string fname,vtkSmartPointer<vtkImageData> imageData);




	void WritePositionLog(std::string saveDir,std::string fname, std::vector<float> &vec);
	template<typename T> void WritePositionLog(std::string fname, std::vector<T> &vec){

		std::ofstream myfile;
		 std::string filename = saveDir;
		 filename.append(fname.begin(),fname.end());
		 std::cout << "WritePositionLog Filename: " << fname << std::endl;
		 myfile.open (filename.c_str());

		     int sizeVec = vec.size();

		     for(int i = 0;i<sizeVec;i+=3){
		    	 myfile << vec[i] << "\t" << vec[i+1] << "\t" << vec[i+2]<< std::endl;
		     }
		   myfile.close();

	}

	void WriteAbsorbtionLog(std::string saveDir,std::string fname, std::vector<float> &vec);
	template<typename T> void WriteAbsorbtionLog(std::string fname, std::vector<T> &vec){

		std::ofstream myfile;
			 std::string filename = saveDir;
			 filename.append(fname.begin(),fname.end());
			 std::cout << "WriteAbsorbtionLog Filename: " << fname << std::endl;
			 myfile.open (filename.c_str());

			     int sizeVec = vec.size();

			     for(int i = 0;i<sizeVec;i+=4){
			    	 myfile << vec[i] << "\t" << vec[i+1] << "\t" << vec[i+2] << "\t" << vec[i+3]<< std::endl;
			     }
			   myfile.close();
	}

	void WriteConstInFile(std::string saveDir,std::string fname, std::string constName, float constant);
	template<typename T> void WriteConstInFile(std::string fname, std::string constName, T constant){

		std::ofstream myfile;
		 std::string filename = saveDir;
		 filename.append(fname.begin(),fname.end());
		 myfile.open (filename.c_str());
		 std::cout << "WriteConstInFile Filename: " << fname << std::endl;
		 myfile << constName << "\t" << constant << std::endl;

		 myfile.close();

	}

	void WritePNG(std::string fname,vtkSmartPointer<vtkImageData> imageData);


	template<typename T> void Printarray(std::string header, std::vector<T> arrIn){

		std::cout << header << ": ";
		for(T &arr : arrIn){

			std::cout << arr << " ";

		}

		std::cout << "" << std::endl;

	}

	Logs();
	virtual ~Logs();


private:





};

#endif /* LOGS_HPP_ */
