/*
 * Emitter.hpp
 *
 *  Created on: Jul 11, 2015
 *      Author: matthias
 */

#ifndef EMITTER_HPP_
#define EMITTER_HPP_

#include <string.h>
#include <vector>
#include "vtkPolyData.h"
#include "vtkDataArray.h"


class Emitter {

	bool verbose;
	vtkPolyData *mesh;
	vtkDataArray* normals;
	float distance;
	float angleIn;
	bool isRadiants;

	std::vector<float> position;
	std::vector<float> angles_arr; // phi:[0:90]
	std::vector<float> theta_arr; // theta:[0:360]
	std::vector<float> relRadiantIntensity_arr;

	std::vector<float> startVector;
	std::vector<int> startVecorNumOcc; // number of occurrences per vector

	std::vector<float> emittermodel;
	std::vector<float> normal; // (n_x,n_y,n_z)
	std::vector<float> surfaceTriangle; // {(p0_x,p0_y,p0_z),(p1_x,p1_y,p1_z),(p2_x,p2_y,p2_z)}

	std::vector<float> emitterHeading;
	std::vector<float> eX,eY,eZ;

	std::vector<float> p0; // point in plane
	int numberOfPhotons;



public:
	// Set
	void SetVerbosity(bool verb);
	void SetVerbosityTrue();
	void SetVerbosityFalse();
	void SetEmitterHeading(std::vector<float> &heading);

	void SetDistanceToSurface(float distance);
	void SetAngleToSurface(float angleIn);
	void SetPointInPlane(std::vector<float> &point0);
	void AngleIsDegree();
	void AngleIsRadiants();
	void SetNumberOfPhotons(int photons);
	void SetMesh(vtkPolyData *mesh);
	void SetNormals(vtkDataArray* normalsIn);
	void SetRadiationCharacteristics(std::vector<float> &angles, std::vector<float> &radiationIntensities);

	// Get
	void GetEmittermodel(std::vector<float> &startPosition,std::vector<float> &startVectors,std::vector<int> &startVecOcc);
	void GetBasisVectors(std::vector<float> &veX,std::vector<float> &veY,std::vector<float> &veZ);


	template<typename T> void Printarray(std::string header, std::vector<T> &arrIn);
	void PrintEmitter();


	Emitter();
	virtual ~Emitter();



private:
	bool IsInCell(std::vector<float> &position, double cellpoint[3],double normal[3], double maxDistance);
	void CalculatePointOfOrigin();
	void CalculateStartVectorOccurrences();
	void SelectNormal();
	void CalculateStartVectorArray();
	void CalculateEmitterModel();
	int TriID(vtkDataArray* normalsIn,std::vector<float> &point);
	void ThetaArr();

	template<typename T> void RotateVector(std::vector<T> vIn, std::vector<T> &vOut, T);
	template<typename T> T SumVec(std::vector<T> &vecIn);

	template<typename T> T Dot(std::vector<T> &a,std::vector<T> &b);
	template<typename T> void Cross(std::vector<T> &a,std::vector<T> &b, std::vector<T> &result);
	template<typename T> void VecMultiply(std::vector<T> &a,std::vector<T> &b, std::vector<T> &result);
	template<typename T> void VecMultiply(T a,std::vector<T> &b, std::vector<T> &result);
	template<typename T> void VecAdd(std::vector<T> &a,std::vector<T> &b, std::vector<T> &result);
	template<typename T> void VecAdd(std::vector<T> &a,std::vector<T> &b,std::vector<T> &c, std::vector<T> &result);
	template<typename T> void Normalize(std::vector<T> &a);

};





#endif /* EMITTER_HPP_ */
