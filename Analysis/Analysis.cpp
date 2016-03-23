/*
 * Analysis.cpp
 *
 *  Created on: Jul 5, 2015
 *      Author: matthias
 */

#include "Analysis.hpp"
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "vtkCell.h"
#include  "vtkDataArray.h"
#include <vector>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

Analysis::Analysis() {
	// TODO Auto-generated constructor stub

	verbose = false;
	chiSquared = 0;
	numSamples = 0;
	numBins = 0;
	numberPhotons = 0;

	fetalRxySum = 0;
	surfRxySum = 0;
	transTxySum = 0;

	normals = vtkDoubleArray::New();
	mesh = vtkPolyData::New();

}

Analysis::~Analysis() {
	// TODO Auto-generated destructor stub
	cout << "Destructing Analysis" << endl;
	normals->Delete();
	mesh->Delete();

}



void Analysis::CalcHistogram(int bins){

	histogram_arr.clear();
	int s = random_arr.size();
	std::vector<float> arr;
	arr.assign(bins,0);
	if(verbose){
	std::cout << "Rand3Array Size: " << s << std::endl;
	}
	for(int i = 0;i<s;i++){
		if(random_arr[i]>=0){
		int z = int(random_arr[i]*bins);

		for(int j=0;j<bins;j++){

			if(z==j){
				arr[j]++;
			}
		}
		}
	}

	histogram_arr.assign(arr.begin(),arr.end());
	numBins = histogram_arr.size();

}


void Analysis::CalcHistogram(std::vector<float> &randomArrIn, int bins){


	SetRandomArray(randomArrIn);

	CalcHistogram(bins);


}




void Analysis::CalcHistogram(std::vector<float> &randomArrIn, int bins, std::vector<float> &histArrOut){

	CalcHistogram(randomArrIn,bins);

	histArrOut.assign(histogram_arr.begin(),histogram_arr.end());

}

void Analysis::CalcChiSquared(std::vector<float> &histArrIn, int NumSamples){

	SetHistogramArray(histArrIn);
	numSamples = NumSamples;
	CalcChiSquared();

}

void Analysis::CalcChiSquared(){


	float e = float(numSamples)/float(numBins);
	if(verbose){
	std::cout << "Bins: " << numBins << " | Samples/Bin: " << e << std::endl;
	}
	chiSquared = 0;

	for(int k = 0;k<numBins;k++){

	/*	cout << "vec[" << k << "]" << ": " << vec[k] << " | ";
		cout << "e: " << e << " | ";
		cout << "pow(): " <<  pow(vec[k]-e,2)/e << "" << endl;
*/
		chiSquared+=(pow(histogram_arr[k]-e,2)/e);

	}

}

void Analysis::GetHistogram(std::vector<float> &hist){

	hist.assign(histogram_arr.begin(),histogram_arr.end());


}

float Analysis::GetChiSquared(std::vector<float> &histArrIn, int numSamples){


	CalcChiSquared(histArrIn,numSamples);
	return chiSquared;

}

float Analysis::GetChiSquared(){
	return chiSquared;
}

void Analysis::CalculateSurfaceSums(){


	for(float &Rxy : surfRxyLog){
		surfRxySum+=Rxy;
	}
	for(float &Rxy : transTxyLog){
		transTxySum+=Rxy;
	}
	for(float &Rxy : fetalRxyLog){
		fetalRxySum+=Rxy;
	}
}

void Analysis::GenerateSurfaceImage(){

// this->surfRxyLog :
	// idx = r*360/gridDims+ phi/gridDim









}




/*
void Analysis::CalculateSurfaceSums(){



	double maxDist = 0.001;
	int numSurfIDs = surfaceIDs.size(); // indices of triangles
	int numLog = absorbtionLog.size();

	int numSurfC = surfaceIDnamesConfig.size();
	surfIDSums.assign(numSurfC,0);
	surfIDSumsNum.assign(numSurfC,0);

	if(verbose){

		cout << "numSurfIDNamesConfig: " << numSurfC << " | { ";
		for(int snc = 0;snc<numSurfC;snc++){
			cout << surfaceIDnamesConfig[snc] << " ";
		}
		cout << "}" << endl;

		cout << "numAbsorbtionLog: " << numLog << endl;
		cout << "numSurfIDs: " << numSurfIDs << " | { ";
		for(int snc = 0;snc<numSurfIDs;snc++){
			cout << surfaceIDs[snc] << " ";
		}
		cout << "}" << endl;
		cout << "##############################################" << endl;
	}

	for(int i=4;i<=numLog;i+=4){

		int LogID0 = i-1; // = (i-4)+3
		int LogID1 = i+3;


		if(verbose){
	//	cout << "LogIDs: 0: " << LogID0 << " | 1: " << LogID1 << " Absorbtions: " << absorbtionLog[LogID0] << " | " <<  absorbtionLog[LogID1] << endl;
		}

		if(absorbtionLog[LogID0]>0 && absorbtionLog[LogID1]==0){
			double position[3];
			position[0] = absorbtionLog[LogID1-3];
			position[1] = absorbtionLog[LogID1-2];
			position[2] = absorbtionLog[LogID1-1];

			for(int j = 0;j<numSurfIDs;j++){
				double n[3];
				double cellp[3];
				int id = surfaceIDs[j];

				mesh->GetCell(id)->GetPoints()->GetPoint(0,cellp);
				normals->GetTuple(id,n);

				if(verbose){
					cout << "i , j : " << i << " , " << j << "" << endl;
					cout << "Position: " << position[0] << ", " << position[1] << ", " << position[2] << endl;
					cout << "CellPoint(" << id << "): " << cellp[0] << ", " << cellp[1] << ", " << cellp[2] << endl;
					cout << "Normal: " << n[0] << ", " << n[1] << ", " << n[2] << endl;

				}


				if(IsInCell(position,cellp,n,maxDist)){

					//std::string surfName= surfaceIDnames[j];

					int sID = GetSurfaceID(surfaceIDnames[j]);
					//cout << "sID: " << sID << endl;
					if(verbose){
					cout << "surfIDSums[" << sID << "]: " << surfIDSums[sID] << endl;
					}

					surfIDSums[sID]+= absorbtionLog[LogID0];

					if(verbose){
					cout << "absorbtionLog[" << LogID0 << "]: " << absorbtionLog[LogID0] << endl;
					cout << "surfIDSums[" << sID << "]: " << surfIDSums[sID] << endl;
					cout << "surfIDSumsNum[" << sID << "]: " << surfIDSumsNum[sID] << endl;
					}
					surfIDSumsNum[sID]+=1;


					break;

				}


			}


		}


	}




}
*/

int Analysis::GetSurfaceID(std::string surfaceNameIn){


	int num = surfaceIDnamesConfig.size();
	int id = num;
	for(int i=0;i<num;i++){
		if(strcmp(surfaceNameIn.c_str(),surfaceIDnamesConfig[i].c_str())==0){
			id = i;
		}
	}
	return id;
}

bool Analysis::IsInCell(double position[3], double cellpoint[3],double normal[3], double maxDistance){

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




double Analysis::CalculateRd(std::vector<float> &vec, int numPhotons, float top){
/*
 * if(absorbtionLog[i-1][3]>0 && absorbtionLog[i][3]==0){
 *
 * 	position[3] = {absorbtionLog[i][0:2]};
 *
 * 	RTges+=absorbtionLog[i-1][3];
 *
 * 	}
 */


	int sizeVec = vec.size();
	double sumVec = 0;
	for(int i = 0;i<sizeVec;i+=4){

		if(i+4+3<=sizeVec){
			float nextWeight = vec[i+4+3];
			float nextZ = vec[i+4+2];
			float weight = vec[i+3];

			if(nextZ==top && weight>0 && nextWeight<=0){
				sumVec += double(weight);
			}
		}
	}

return sumVec/=numPhotons;
}

double Analysis::CalculateTd(std::vector<float> &vec, int numPhotons, float bottom){

	int sizeVec = vec.size();
	double sumVec = 0;
	for(int i = 0;i<sizeVec;i+=4){

		if(i+4+3<=sizeVec){
			float nextWeight = vec[i+4+3];
			float nextZ = vec[i+4+2];
			float weight = vec[i+3];

			if(nextZ==bottom && weight>0 && nextWeight<=0){
		    	sumVec += double(weight);
			}
		}
	}

return sumVec/=numPhotons;
}

void Analysis::CalculateRTd(std::vector<float> &vec,std::vector<float> &pz, float top, float bottom){

	int sizeVec = vec.size();
	double sumRd = 0;
	double sumTd = 0;
	std::vector<double> arr;
	arr.assign(2,0);
	int numR = 0;
	int numT = 0;
	for(int i = 0;i<sizeVec;i++){

		if(fabs(pz[i]-top)<=0.00001){
			numR++;
			sumRd+=vec[i];

		}else if(fabs(pz[i]-bottom)<=0.00001){
			numT++;
			sumTd+=vec[i];

		}else{
			std::cout << " else pz[" << i << "]: " << pz[i] <<std::endl;

		}
	}

	if(verbose){
		std::cout << "Top: " << top << std::endl;
		std::cout << "Bottom: " << bottom << std::endl;
		std::cout << "SizeVec: " << sizeVec << std::endl;
		std::cout << "NumR: " << numR << std::endl;
		std::cout << "NumT: " << numT << std::endl;
		std::cout << "sumRd: " << sumRd << std::endl;
		std::cout << "sumTd: " << sumTd << std::endl;
	}
	arr[0]=sumRd/sizeVec;
	arr[1]=sumTd/sizeVec;

	RdTd.assign(arr.begin(),arr.end());

}

void Analysis::CalculateRTd(std::vector<float> &vec,std::vector<float> &pz, float top, float bottom, std::vector<double> &out){

	CalculateRTd(vec,pz,top,bottom);

	out.assign(RdTd.begin(),RdTd.end());

}

void Analysis::SetNumberOfPhotons(int num){
	this->numberPhotons = num;
}
void Analysis::SetVerbosity(bool verb){
	verbose = verb;
}
void Analysis::SetVerbosityTrue(){
	verbose = true;
}
void Analysis::SetVerbosityFalse(){
	verbose = false;
}

void Analysis::SetRxyLogs(std::vector<float> &surfRxy, std::vector<float> &fetalRxy, std::vector<float> &transRxy){

	surfRxyLog.assign(surfRxy.begin(),surfRxy.end());
	fetalRxyLog.assign(fetalRxy.begin(),fetalRxy.end());
	transTxyLog.assign(transRxy.begin(),transRxy.end());


}

void Analysis::SetRxyLogs(std::vector<float> &surffetaltransRxy){

	int s = surffetaltransRxy.size();


	for(int i = 0;i<s;i+=3){

		surfRxyLog.insert(surfRxyLog.end(),surffetaltransRxy[i]);
		fetalRxyLog.insert(fetalRxyLog.end(),surffetaltransRxy[i+1]);
		transTxyLog.insert(transTxyLog.end(),surffetaltransRxy[i+2]);


	}
}
void Analysis::SetRTdxyLogs(vtkSmartPointer<vtkImageData> image){

	int ext[6];
	image->GetExtent(ext);


	for(int x=ext[0];x<=ext[1];x++){
		for(int y=ext[2];y<=ext[3];y++){

			surfRxyLog.insert(surfRxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,0));
			fetalRxyLog.insert(fetalRxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,1));
			transTxyLog.insert(transTxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,2));
		}
	}


}

void Analysis::SetAbsorbtionLogs(std::vector<float> &absorbtionXZLogIn,std::vector<float> &absorbtionYZLogIn){

	absorbXZLog.assign(absorbtionXZLogIn.begin(),absorbtionXZLogIn.end());
	absorbYZLog.assign(absorbtionYZLogIn.begin(),absorbtionYZLogIn.end());

}



void Analysis::SetRandomArray(std::vector<float> &randomArrIn){

	random_arr.assign(randomArrIn.begin(),randomArrIn.end());
	numSamples = random_arr.size();

}
void Analysis::SetHistogramArray(std::vector<float> &histArrIn){

	histogram_arr.assign(histArrIn.begin(),histArrIn.end());
	numBins = histogram_arr.size();

}
void Analysis::SetSurfaceIDs(std::vector<int> &surfIDsIn){

	surfaceIDs.assign(surfIDsIn.begin(),surfIDsIn.end());

}
void Analysis::SetSurfaceIDNames(std::vector<std::string> &surfIDNamesIn){

	surfaceIDnames.assign(surfIDNamesIn.begin(),surfIDNamesIn.end());

}
void Analysis::SetSurfaceIDNamesConfig(std::vector<std::string> &surfIDNamesIn){

	surfaceIDnamesConfig.assign(surfIDNamesIn.begin(),surfIDNamesIn.end());
	surfaceIDnamesConfig.push_back("noname");

}



void Analysis::SetMesh(vtkPolyData* meshIn){

	mesh->DeepCopy(meshIn);

}

void Analysis::SetNormals(vtkDataArray* normalsIn){
	 // Normals per Triangle: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3

	normals->DeepCopy(normalsIn);

}

void Analysis::PrintChiSquared(){

	std::cout << "Chi squared: " << GetChiSquared() << std::endl;

}
void Analysis::PrintRTd(){

	std::cout << "R_d: " << RdTd[0] << std::endl;
	std::cout << "T_d: " << RdTd[1] << std::endl;
}

void Analysis::PrintSurfaceSums(){

	cout << "Number of Photons: " << numberPhotons << endl;
	cout << "SurfRxy: " << surfRxySum  << " | in %: " << surfRxySum/numberPhotons*100<< endl;
	cout << "TransRxy: " << transTxySum << " | in %: " << transTxySum/numberPhotons*100<< endl;
	cout << "FetalRxy: " << fetalRxySum << " | in %: " << fetalRxySum/numberPhotons*100<< endl;
	cout << "Surf+TransRxy: " << surfRxySum + transTxySum << endl;
	cout << "############################################" << endl;

}

void Analysis::PrintAnalysis(){
	cout << "################ Analysis ##################" << endl;
	PrintSurfaceSums();


}


/*
void Analysis::PrintSurfaceSums(){


	int num = surfIDSums.size();

	cout << "############################################\nSurface sums:" << endl;
	cout << "NumSurfIDSums: " << num << endl;
	cout << "Number of Photons: " << numberPhotons << endl;
	cout << "NumsurfaceIDnamesConfig: " << surfaceIDnamesConfig.size() << endl;
	for(int i=0;i<num;i++){

		cout << "Num(" << surfaceIDnamesConfig[i] << "): " << surfIDSumsNum[i] << " | Sum: " << surfIDSums[i] << " | in \%: " << surfIDSums[i]/numberPhotons << endl;
	}

	cout << "############################################" << endl;

}
*/
