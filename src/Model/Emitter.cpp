/*
 * Emitter.cpp
 *
 *  Created on: Jul 11, 2015
 *      Author: matthias
 */

#include "Emitter.hpp"
#include <math.h>
#include <functional>
#include <algorithm>

#include "vtkCell.h"


Emitter::Emitter() {
	// TODO Auto-generated constructor stub

	verbose = false;
	distance = 1;
	angleIn = 0;
	numberOfPhotons = 1;

	mesh = vtkPolyData::New();
	normals = vtkDataArray::CreateDataArray(VTK_DOUBLE);
	isRadiants = true;
	p0.assign(3,0);


}

Emitter::~Emitter() {
	// TODO Auto-generated destructor stub

	cout << "Destructing Emitter" << endl;



}

void Emitter::SetVerbosity(bool verb){
	verbose = verb;
}
void Emitter::SetVerbosityTrue(){
	verbose = true;
}
void Emitter::SetVerbosityFalse(){
	verbose = false;
}
void Emitter::AngleIsDegree(){
	isRadiants = false;
}
void Emitter::AngleIsRadiants(){
	isRadiants = true;
}
void Emitter::SetNumberOfPhotons(int photons){
	this->numberOfPhotons = photons;
}


void Emitter::SetDistanceToSurface(float dist){
	this->distance = dist;
}
void Emitter::SetAngleToSurface(float angIn){
	this->angleIn = angIn;
}
void Emitter::SetPointInPlane(std::vector<float> &point0){
	this->p0.assign(point0.begin(),point0.end());
}

void Emitter::SetRadiationCharacteristics(std::vector<float> &angles, std::vector<float> &radiationIntensities){
	angles_arr.assign(angles.begin(),angles.end());
	relRadiantIntensity_arr.assign(radiationIntensities.begin(),radiationIntensities.end());

}


void Emitter::SetMesh(vtkPolyData* meshIn){

	mesh->DeepCopy(meshIn);

}
void Emitter::SelectNormal(){
	 // Normals per Triangle: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3


	std::vector<float> n;
	n.assign(3,0);

	if(verbose){
		cout << "Select Normal" << endl;

	}

	int id = TriID(normals,p0);
	double n0[3];

	normals->GetTuple(id,n0);

	n[0] = float(n0[0]); //n0_x
	n[1] = float(n0[1]); //n0_y
	n[2] = float(n0[2]); //n0_z
	if(verbose){
		cout << "Select Normal(id=" << id << ") = (" << n[0] << ", " << n[1] << ", " << n[2] << ")" << endl;
	}

	this->normal.assign(n.begin(),n.end());

}

void Emitter::SetNormals(vtkDataArray* normalsIn){
	normals->DeepCopy(normalsIn);

}

void Emitter::SetEmitterHeading(std::vector<float> &heading){
	emitterHeading.assign(heading.begin(),heading.end());
}


int Emitter::TriID(vtkDataArray* normalsIn,std::vector<float> &point){

	int id = -1;

	int numCells = mesh->GetNumberOfCells();


	for(int i = 0;i<numCells;i++){
		double n[3];
		double cellp[3];

		mesh->GetCell(i)->GetPoints()->GetPoint(0,cellp);
		normalsIn->GetTuple(i,n);

		if(IsInCell(point,cellp,n,0.01)){
			id = i;
			break;
		}

	}

	return id;

}

bool Emitter::IsInCell(std::vector<float> &position, double cellpoint[3],double normal[3], double maxDistance){

	bool r = false;

	double dist = fabs(normal[0]*(position[0]-cellpoint[0])+normal[1]*(position[1]-cellpoint[1])+normal[2]*(position[2]-cellpoint[2]));

	if(verbose){
		cout << "distance to plane: " << dist << endl;
	}


	if(dist<=maxDistance){

		r = true;
	}
	return r;
}

void Emitter::ThetaArr(){

	theta_arr.assign(360,0);

	for(int i = 0;i<360;i++){

		if(isRadiants){
			theta_arr[i] = i*M_PI/180.0;
		}else {
			theta_arr[i] = i;
		}

	}

	if(verbose){

		Printarray("Theta",theta_arr);

	}

}

void Emitter::CalculatePointOfOrigin(){

	position.assign(3,0);
	if(verbose){
		cout << "#####################################\n";
		cout << "CalculatePointOfOrigin...\n" << endl;
	}
	position[0] = -emitterHeading[0]*distance + p0[0];
	position[1] = -emitterHeading[1]*distance + p0[1];
	position[2] = -emitterHeading[2]*distance + p0[2];


	if(verbose){

		cout << "Position: (" << position[0] << ", " << position[1] << ", " << position[2] << ")\n";
		cout << "\nCalculatePointOfOrigin done\n";
		cout << "#####################################\n" << endl;
	}

}

void Emitter::CalculateStartVectorArray(){
	if(verbose){
			cout << "#####################################\n";
			cout << "CalculateStartVectorArray... \n" << endl;
		}
	std::vector<float> v0;
	v0.assign(emitterHeading.begin(),emitterHeading.end());
	//v0.resize(emitterHeading.size());
	//std::transform(normal.begin(),normal.end(),v0.begin(),std::negate<float>());

	if(verbose){
		// print variables
		cout << "Normal | inverse: (";
		for( auto &n: normal){cout << n << " ";}
		cout << ") | (";
		for(auto &inv : v0 ){cout << inv << " ";}
		cout << ")" << endl;
		cout << "anglesSize: " << angles_arr.size() << endl;
	}

		for(auto &phi : angles_arr){
			if(phi==0){
				//std::vector<float> v;
				startVector.insert(startVector.end(),v0.begin(),v0.end());

			}else{
				std::vector<float> v;
				RotateVector(v0,v,phi);
				startVector.insert(startVector.end(),v.begin(),v.end());
			}

		}


	if(verbose){
		Printarray("startVector", startVector);
		cout << "\nCalculateStartVectorArray done\n";
		cout << "#####################################\n" << endl;
	}


}

void Emitter::CalculateStartVectorOccurrences(){
	if(verbose){
		cout << "#####################################\n";
		cout << "CalculateStartVectorOccurrences... \n" << endl;
	}
float sum = SumVec(relRadiantIntensity_arr);
int s1 = 0;
	for(auto &i : relRadiantIntensity_arr){

		int occ = int(numberOfPhotons*i/sum);
		s1 += occ;
		startVecorNumOcc.insert(startVecorNumOcc.end(),occ);

	}

	if(s1< numberOfPhotons){
		startVecorNumOcc[0]+=(numberOfPhotons-s1);
	}

	if(verbose){
			// print variables
			cout << "Sum(relRadIntensities):" << sum << "" << endl;

			Printarray("relRadIntensities", relRadiantIntensity_arr);
			Printarray("startVectorOcc", startVecorNumOcc);
			cout << "\nCalculateStartVectorOccurrences done\n";
			cout << "#####################################\n" << endl;
	}
}



void Emitter::GetEmittermodel(std::vector<float> &pStart, std::vector<float> &vStart, std::vector<int> &vStartOcc){
	CalculateEmitterModel();

	pStart.assign(position.begin(),position.end());
	vStart.assign(startVector.begin(),startVector.end());
	vStartOcc.assign(startVecorNumOcc.begin(),startVecorNumOcc.end());


}
void Emitter::GetBasisVectors(std::vector<float> &veX,std::vector<float> &veY,std::vector<float> &veZ){
	veX.assign(eX.begin(),eX.end());
	veY.assign(eY.begin(),eY.end());
	veZ.assign(eZ.begin(),eZ.end());

}


void Emitter::CalculateEmitterModel(){

	// CalcP0();
	SelectNormal();

	CalculatePointOfOrigin();
	CalculateStartVectorArray();
	CalculateStartVectorOccurrences();


	// calc Basis vectors eX,eY,eZ

	eZ.assign(emitterHeading.begin(),emitterHeading.end());

	std::vector<float> a;
	a.assign(3,0);
	if(eZ[2]!=0){
		a[0] = 1;
		a[1] = 1;
		a[2] = -(eZ[0]+eZ[1])/eZ[2];
	}else if(eZ[1]!=0){
		a[0] = 1;
		a[1] = -eZ[0]/eZ[1];
		a[2] = 0;
	}else if(eZ[0]!=0){
		a[0] = 0;
		a[1] = 1;
		a[2] = 0;
	}
	Normalize(a);
	Normalize(eZ);
	//cout << "eZ: " << eZ[0] << ", " << eZ[1] << ", " << eZ[2] << endl;
	Cross(eZ,a,eY);
	//cout << "eY: " << eY[0] << ", " << eY[1] << ", " << eY[2] << endl;
	Normalize(eY);
	Cross(eY,eZ,eX);

	Normalize(eX);
	//cout << "eX: " << eX[0] << ", " << eX[1] << ", " << eX[2] << endl;

}

template<typename T> void Emitter::Printarray(std::string header, std::vector<T> &arrIn){

	cout << header << ": ";
	for(T &arr : arrIn){

		cout << arr << " ";

	}

	cout << "" << endl;

}
void Emitter::PrintEmitter(){
	cout << "########### Emitter ##################" << endl;
			Printarray("vStart: ", startVector);
			Printarray("pStart: ", position);
			Printarray("POI: ", p0);
			Printarray("eX", eX);
			Printarray("eY", eY);
			Printarray("eZ", eZ);

}
//
//void Emitter::Printarray(std::string header, std::vector<int> &arrIn){
//
//	cout << header << ": ";
//	for(auto &arr : arrIn){
//
//		cout << arr << " ";
//
//	}
//
//	cout << "" << endl;
//
//}


template<typename T> T Emitter::SumVec(std::vector<T> &vecIn){

	T sum = 0;

	for(auto &i : vecIn){
		sum+=i;
	}

	return sum;
}

template<typename T> void Emitter::RotateVector(std::vector<T> vIn, std::vector<T> &vOut,T phi){
// vIn = normal vector
	vOut.resize(vIn.size());

	if(!isRadiants){

		phi*=(M_PI/180.0); // phi(0°, 90°)
	}

	T cosPhi = cos(phi); // phi(0°, 90°)
	T sinPhi = sqrt(1-cosPhi*cosPhi);


// Rotate:
// plane: vIn[0]*x + vIn[1]*y + vIn[2]*z = 0;
	T x,y,z;
	std::vector<T> e_rot;
	e_rot.assign(3,0);
	if(vIn[2]!=0){
		z = -(vIn[0] + vIn[1])/vIn[2];
		e_rot[0] = 1;
		e_rot[1] = 1;
		e_rot[2] = z;

	}else if(vIn[1]!=0){
		y = -(vIn[0] + vIn[2])/vIn[1];
		e_rot[0] = 1;
		e_rot[1] = y;
		e_rot[2] = 1;
	}else if(vIn[0]!=0){
		x = -(vIn[1] + vIn[2])/vIn[0];
		e_rot[0] = x;
		e_rot[1] = 1;
		e_rot[2] = 1;
	}

	Normalize(e_rot);


	std::vector<T> nx;
	std::vector<T> nxERot;
	std::vector<T> vecmult;
	std::vector<T> vecMultCosT;
	std::vector<T> vecMultSinT;
	std::vector<T> vNew;

	Cross(e_rot,vIn,nx);
	Cross(nx,e_rot,nxERot);
	T nDx = Dot(e_rot,vIn);
	VecMultiply(nDx,e_rot,vecmult);
	VecMultiply(cosPhi,nxERot,vecMultCosT);
	VecMultiply(sinPhi,nx, vecMultSinT);
	VecAdd(vecmult,vecMultCosT, vecMultSinT,vNew);

	Normalize(vNew);

	if(verbose){
		cout << "\tPhi: " << phi << "" << endl;
		Printarray("\tvIn", vIn);
		Printarray("\te_rot",e_rot);
		Printarray("\tvNew",vNew);
	}

/*
	std::vector<T> v;
	v.assign(3,0);
	U R_13[3] = {cosTheta*cosPhi,-sinTheta,cosTheta*sinPhi};
	U R_23[3] = {sinTheta*cosPhi, cosTheta,sinTheta*sinPhi};
	U R_33[3] = {-sinPhi,0,cosPhi};


	v[0] = T(R_13[0]*vIn[0]+R_13[1]*vIn[1]+R_13[2]*vIn[2]);
	v[1] = T(R_23[0]*vIn[0]+R_23[1]*vIn[1]+R_23[2]*vIn[2]);
	v[2] = T(R_33[0]*vIn[0]+R_33[1]*vIn[1]+R_33[2]*vIn[2]);
*/
	vOut.assign(vNew.begin(),vNew.end());

}


template<typename T> T Emitter::Dot(std::vector<T> &a,std::vector<T> &b){

	T dotprod = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];

return dotprod;
}
template<typename T> void Emitter::Cross(std::vector<T> &a,std::vector<T> &b,std::vector<T> &result){

	//std::vector<T> vec;
		result.resize(a.size());

	result[0] = a[1]*b[2] - a[2]*b[1];
	result[1] = a[2]*b[0] - a[0]*b[2];
	result[2] = a[0]*b[1] - a[1]*b[0];

	//result.assign(vec.begin(),vec.end());
}


template<typename T> void Emitter::VecMultiply(std::vector<T> &a,std::vector<T> &b,std::vector<T> &result){

	std::vector<T> vec;
	vec.resize(a.size());
	vec[0] = a[0]*b[0];
	vec[1] = a[1]*b[1];
	vec[2] = a[2]*b[2];

	result.assign(vec.begin(),vec.end());
}
template<typename T> void Emitter::VecMultiply(T a,std::vector<T> &b,std::vector<T> &result){

	std::vector<T> vec;
	vec.resize(b.size());
	vec[0] = a*b[0];
	vec[1] = a*b[1];
	vec[2] = a*b[2];

	result.assign(vec.begin(),vec.end());
}
template<typename T> void Emitter::VecAdd(std::vector<T> &a,std::vector<T> &b,std::vector<T> &result){

	std::vector<T> vec;
	vec.resize(a.size());

	vec[0] = a[0]+b[0];
	vec[1] = a[1]+b[1];
	vec[2] = a[2]+b[2];

	result.assign(vec.begin(),vec.end());
}
template<typename T> void Emitter::VecAdd(std::vector<T> &a,std::vector<T> &b,std::vector<T> &c,std::vector<T> &result){
	std::vector<T> vec;
	vec.resize(a.size());

	vec[0] = a[0]+b[0]+c[0];
	vec[1] = a[1]+b[1]+c[0];
	vec[2] = a[2]+b[2]+c[0];

	result.assign(vec.begin(),vec.end());

}
template<typename T> void Emitter::Normalize(std::vector<T> &a){
	std::vector<T> vec;
	vec.resize(a.size());

	T f = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

	vec[0] = a[0]/f;
	vec[1] = a[1]/f;
	vec[2] = a[2]/f;

	a.assign(vec.begin(),vec.end());

}

