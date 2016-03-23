/*
 * CalculatePropagation.hpp
 *
 *  Created on: Mar 15, 2015
 *      Author: matthias
 */

#ifndef CALCULATEPROPAGATION_HPP_
#define CALCULATEPROPAGATION_HPP_

#include <CL/cl.hpp>
#include "../PSNR/PSNR.hpp"
//#include "SDKUtil/CLUtil.hpp"
#include "vtkPolyData.h"
#include "vtkFloatArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include <vector>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

// #include <array>



class CalculatePropagation {

	bool verbose;
	//vtkPolyData *PolyImage;
	bool showPercent;
	long startTime;
	long endTime;
	long neededTime;


	long startTimeTest;
	long endTimeTest;
	cl_device_type devType;
	int numTest;
	bool isReadOnly;
	bool useCPU;
	bool useGPU;
	bool isContextSet;
	cl_int status;
	const char *filename;
	const char* kernelname;
	const char* buildOptions;
	std::string sourceStr;

	 int photonsPerLoop;
	std::string testSourceStr;



	cl::Context context;
	std::vector<cl::Device> devices;
	cl::Device selectedDevice;
	cl::CommandQueue cmd_queue;
	cl::CommandQueue testcmd_queue;
	std::vector<cl::Platform> platforms;
	cl::Program  program;

	cl::Program  Testprogram;

/*	size_t wgSize;
	size_t lWGSize;
	size_t gWSize[3];
*/
	cl::Kernel kernel;
	typedef cl::make_kernel<cl_int, //numTriangles
			cl_float, // factor
			cl::Buffer&, // vStart
			cl::Buffer&, // vStartOcc
			cl_float3, // pStart
			cl_float, // chance
			cl_float, // w_thresh
			cl_int, // fetalTissueID
			cl::Buffer&, //d_SurfFetalTransRxy
			cl::Buffer&, //d_RecordAxzAyz
			cl::Buffer&, //d_DimsXYZ
			cl::Buffer&, //e_P0XYZRPhi(Rho?)
			cl::Buffer&, //gridDims
			cl::Buffer&, // triangles
			cl::Buffer&, // normals
			cl::Buffer&, // mediumIDs
			cl::Buffer&, // n mua, mus,g
			cl::Buffer&, // rand3Arr
			cl::Buffer&> KernelType;// inextInxtpPoolIDLock


	typedef cl::make_kernel<cl_int,
			cl::Buffer&,
			cl::Buffer&,
			cl::Buffer&> TestKernel;

	// Test
	std::vector<int> A;
	cl::Buffer d_A;
	std::vector<int> B;
	cl::Buffer d_B;
	std::vector<int> C;
	cl::Buffer d_C;

	cl_int iter;

	//Photon properties
	cl_int numPhotonInteractions;
	cl_int absNumberPhotons;
	cl_int numberPhotons;
	std::vector<float> vStart; // {(v_x,v_y,v_z)_0, (v_x,v_y,v_z)_1,...}
	cl::Buffer d_vStart;
	std::vector<int> vStartOcc; // {(v_x,v_y,v_z)_0, (v_x,v_y,v_z)_1,...}
	cl::Buffer d_vStartOcc;
	std::vector<float> pStart; // {(p_x,p_y,p_z)_0, (p_x,p_y,p_z)_1,...}
	cl_float3 d_pStart;
	cl_float chance;
	cl_float w_threshold;
	cl_float pixelPerUnitFactor;
	cl_int fetalTissueID;

	//Logs

	std::vector<float> DimsXYZ; // (x,y,z)
	cl::Buffer d_DimsXYZ;
	// DetectionGrid
	std::vector<float> DetGridDims; // (x,y,z,r,phi, rho)
	cl::Buffer d_DetGridDims;
	// Coordinate axis
	std::vector<float> eP0XYZ; // (p0x,p0y,p0z,exx,exy,exz,...)
	cl::Buffer d_eP0XYZ;


	//std::vector<float> AbsorbtionLog; // {(x0,y0,z0, w0),(x1,y1,z1,w1),...}
	//cl::Buffer d_AbsorbtionLog;


	std::vector<float> RandNmbersLog; // {(x0,y0,z0),(x1,y1,z1),...} per photon
	cl::Buffer d_RandNmbersLog;

	// Record Transmittance
	//std::vector<float> TransTxy;
	cl::Buffer d_TransTxy;

	// Record Reflectance
	//std::vector<float> SurfRxy;
	cl::Buffer d_SurfRxy;

	//std::vector<float> FetalRxy;
	cl::Buffer d_Fetalxy;




	//std::vector<float> SurfFetalTransRxy;
	vtkSmartPointer<vtkImageData> SurfFetalTransRxyImage;
	cl::Buffer d_SurfFetalTransRxy;

	vtkSmartPointer<vtkImageData> RTdxyImage;
	vtkSmartPointer<vtkImageData> ArzImage;


	// Record Absorbtion in xz- and yz-Planes
	std::vector<float> RecordAxz;
	std::vector<int> RAxzCoords;
	//cl::Buffer d_RecordAxz;
	std::vector<float> RecordAyz;
	std::vector<int> RAyzCoords;
	//cl::Buffer d_RecordAyz;
	std::vector<float> RecordAxzAyz;
	cl::Buffer d_RecordAxzAyz;

	std::vector<int> lockedRecord;
	cl::Buffer d_lockedRecord;

	// Tissue properties and geometry
	std::vector<float> dimensions; //dimensions
	std::vector<float> triangles; //Triangles: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3
	cl::Buffer d_triangles;
	std::vector<float> normals; // Normals per Triangle: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3
	cl::Buffer d_normals;
	std::vector<int> triN; // index of n and n_outer for every triangle {(n0_in,n0_out),(n1_in,n1_out),...}, Size: NumTriangles*2;
	cl::Buffer d_triN;

	std::vector<float> n; //Tissue refraction indices: {n0,n1,n2,...}, Size: NumTissues; specific tissue: n(triN)
	cl::Buffer d_n;
	std::vector<float> mua;
	cl::Buffer d_mua;
	std::vector<float> mus;
	cl::Buffer d_mus;
	std::vector<float> g;
	cl::Buffer d_g;

	std::vector<float> nmuasG;
	cl::Buffer d_nmuasG;

	//PRNG
	 //nach durchlauf soll Array im Device aktualisiert werden
	std::vector<float> rand3Arr;
	cl::Buffer d_rand3Array;

	std::vector<float> randPool;
	cl::Buffer d_randPool;
	cl_int poolSize;
	cl_int seed;
	cl_int numberGenerator;
	std::vector<long> MBMZ;
	cl::Buffer d_MBMZ;
	std::vector<int> inextInextp;
	cl::Buffer d_inextInextp;


public:
	void SetVerbosity(bool verb);
	void SetVerbosityTrue();
	void SetVerbosityFalse();
	void SetShowPercentDone(bool s);
	void SetShowPercentDoneTrue();
	void SetShowPercentDoneFalse();
	// Init
	int Init();

	// Execution
	int Propagation();


	// Set Log Properties
	void SetDimensionsXYZ(float xdim, float ydim, float zdim);
	void SetDimensionsXYZ(std::vector<float> &dims);
	void SetDetectionGridDimensions(float xdim, float ydim, float zdim, float rdim, float phiDim, float rhoDim);
	void SetDetectionGridDimensions(std::vector<float> &detGrid);
	void SetBasisVectors(std::vector<float> &veX,std::vector<float> &veY,std::vector<float> &veZ);



	void SetExVector(float x,float y, float z);
	void SetEyVector(float x,float y, float z);
	void SetEzVector(float x,float y, float z);
	void SetPOI(float x, float y, float z);
	void SetPOI(std::vector<float> &point);
	void SetFetalTissueID(int id);
	// Set OpenCL Properties
	int SetKernelName(const char*);
	int SetSourceString(const char*);
	int SetTestSourceString(const char* fname);
	void SetBuildOptions(const char* options);
	void SetDeviceType(int device);


	/*
	 * Set Simulation properties
	 */
	void SetNumberOfPhotonInteractions(int);
	void SetNumberOfPhotons(int);
	void SetStartTrajectories(std::vector<float> &trajectories, std::vector<int> &trajOcc);
	void SetStartTrajectory(std::vector<float> &trajectories);
	void SetStartTrajectory(float x, float y, float z);
	void SetStartTrajectory(vtkDoubleArray* trajectory, int dims);
	void SetStartPoint(std::vector<float> &point);
	void SetStartPoint(float x, float y, float z);
	void SetStartPoint(vtkDoubleArray* startPoint, int dims);
	void SetChanceForRoulette(float);
	void SetWeightThreshold(float);
	void SetPixelPerUnitFactor(float);


	/*
	 * Get Simulation properties
	 */
	int GetNumberOfPhotons();
	std::vector<float> GetStartTrajectory();
	void GetStartTrajectory(std::vector<float> &trajectory);
	std::vector<float> GetStartPoint();
	void GetStartPoint(std::vector<float> &startPoint);
	float GetChanceForRoulette();
	float GetWeightThreshold();




	/*
	 * Get Logs
	 */
	void GetLogs(std::vector<float> &randoms,
			std::vector<float> &absorbtionsXZ,std::vector<float> &absorbtionsYZ,
			std::vector<float> &surfaceXY,std::vector<float> &fetalXY,
			std::vector<float> &transmittenceXY);

	void GetLogs(std::vector<float> &randoms,
			std::vector<float> &absorbtionsXZ,std::vector<float> &absorbtionsYZ,
			std::vector<float> &surfFetalTransXY);

	float ScalePerRadius(int ir);
	float ScalePerAngle1D(int ia);
	float ScalePerRadius1D(int ir);
	float Scale2D(int ir, int ia);
	float ScaleArz2D(int ir);



	vtkSmartPointer<vtkImageData> GetSFTImage();
	vtkSmartPointer<vtkImageData> GetSFTImage(int);
	vtkSmartPointer<vtkImageData> GetRTdImage();
	vtkSmartPointer<vtkImageData> GetRTdImage(int);

	vtkSmartPointer<vtkImageData> GetArzImage(); // Absorption
	vtkSmartPointer<vtkImageData> GetPhirzImage();//fluenz


	void GetAbsorbtionLog(std::vector<float> &AXZ, std::vector<float> &AYZ);
	void GetAbsorbtionLog(std::vector<float> &vec);
	void GetAbsorbtionLog(vtkDoubleArray* arr);
	std::vector<float> GetAbsorbtionLog();

	void GetAbsorbtionLogCoords(std::vector<int> &XZcoords, std::vector<int> &YZcoords);

	void GetRandomNumbersLog(std::vector<float> &vec);
	void GetRandomNumbersLog(vtkDoubleArray* arr);
	std::vector<float> GetRandomNumbersLog();

	void GetReflectenceLog(std::vector<float> &vec);
	void GetReflectenceLog(vtkDoubleArray* arr);
	std::vector<float> GetReflectenceLog();

	void GetTransmittenceLog(std::vector<float> &vec);
	void GetTransmittenceLog(vtkDoubleArray* arr);
	std::vector<float> GetTransmittenceLog();





	// Set Tissue Properties
	void SetTissueProperties(std::vector<float> &refractionInd, std::vector<float> &absorbtionCoeffs, std::vector<float> &anisotropies, std::vector<float> &scatterCoeffs);

	void SetRefractionIndices(std::vector<float> &refractions); //Tissue refraction indices
	void SetRefractionIndices(vtkDoubleArray* refractions,int dims);
	void SetRefractionIndices(double* refractions,int dims);

	void SetAbsorbtioncoefficients(std::vector<float> &absorbtions); //Absorbtionskoeffizienten
	void SetAbsorbtioncoefficients(vtkDoubleArray* absorbtions,int dims);

	void SetScattercoefficients(std::vector<float> &scatter); // Streuungskoeffizienten
	void SetScattercoefficients(vtkDoubleArray* scatter,int dims);

	void SetAnisotropy(std::vector<float> &anisotropies); //Anisotropie
	void SetAnisotropy(vtkDoubleArray* anisotropies,int dims);

	//Get Tissue Properties
	std::vector<float> GetRefractionindices();
	void GetRefractionindices(std::vector<float> &refractions);
	std::vector<float> GetAbsorbtioncoefficients();
	void GetAbsorbtioncoefficients(std::vector<float> &absorbtions);
	std::vector<float> GetScattercoefficients();
	void GetScattercoefficients(std::vector<float> &scatter);
	std::vector<float> GetAnisotropy();
	void GetAnisotropy(std::vector<float> &Anisotropy);


	// Set Tissue geometry
	void SetDimensions(std::vector<float> &dims);
	void SetTriangles(vtkPolyData* meshimage);
	void SetNormals(vtkDataArray* normals);
	void SetInnerAndOuterTissueIndices(vtkDoubleArray* mediums, double* contours, int numContours);
	void SetInnerAndOuterTissueIndices(vtkDoubleArray* mediums);


	// Set Pseudo Random Number Generator
	void SetNumberGenerator(int);
	std::vector<float> GetRand3Array(); //PRNs
	void GetRand3Array(std::vector<float> &psnrs);

	// Get Pseudo Random Number Generator
	int GetSeed();
	int GetNumberGenerator();



	void printVector(std::string header, const std::vector<float> vec);
	void printVector(std::string header, const std::vector<int> vec);
	void printVector(std::string header, const std::vector<long> vec);
	void PrintTimeSpent();
	void PrintMediumIDs();
	void PrintRefractionIndices();
	void PrintSimulationParams();
	void PrintPercentDone(double p);

	CalculatePropagation();
	virtual ~CalculatePropagation();






private:



	/*
	 * Pseudo Random Number Generator
	 */
	void SetRand3Array(std::vector<float> &arr); //PRNs
	void SetRand3Array(vtkDoubleArray* rands,int dims);
	void SetSeed(int);

	void InitLogVectors();


	int convertToString(const char*, std::string&);

	void FindIndex(double value,double* valuesArray,int arrSize, int &index);

	int setup();
	int run();
	int cleanup();
	cl::Platform getPlatform(cl_device_type type);
	int CreateCommandQueue();
	void CreateBufferObjects();
	void CreateBufferObjects(cl::Context con);
	int createContext();
	void CopyDataToHost(cl::Buffer &buffer, std::vector<float> &vec);
	int ExecKernel();
	int BuildKernel();

	void CopyTestDataToHost(cl::Buffer &buf, std::vector<int> &vec);
	void CreateTestBuffer();

	//int fillRandom(std::vector<float> &vec,const int width,const int height,const float rangeMin,const float rangeMax, unsigned int seed);



};

#endif /* CALCULATEPROPAGATION_HPP_ */
