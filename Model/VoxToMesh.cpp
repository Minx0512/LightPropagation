/*
 * VoxToMesh.cpp
 *
 *  Created on: Mar 9, 2015
 *      Author: matthias
 */

#include "VoxToMesh.hpp"


#include <math.h>

#include "vtkSmartPointer.h"
#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkStripper.h"
#include "vtkImageToPolyDataFilter.h"

#include "vtkOutlineFilter.h"

#include "vtkUnstructuredGrid.h"
#include "vtkAppendPolyData.h"

#include "vtkVolumeProperty.h"
#include "vtkPiecewiseFunction.h"

#include "vtkDecimatePro.h"
#include "vtkQuadricDecimation.h"
#include "vtkDelaunay3D.h"
#include "vtkSmoothPolyDataFilter.h"

#include "vtkProperty.h"
#include "vtkCellData.h"
#include "vtkPoints.h"
#include "vtkCell.h"
#include "vtkGenericCell.h"
#include "vtkPointData.h"



VoxToMesh::VoxToMesh() {


	VoxImage = vtkSmartPointer<vtkImageData>::New();
	IDVoxelImage = vtkSmartPointer<vtkImageData>::New();
	normals = vtkPolyDataNormals::New();
	PolyImage = vtkPolyData::New();
	stripper = vtkStripper::New();

	mapper = vtkPolyDataMapper::New();
	actor = vtkActor::New();
	contour = vtkContourFilter::New();
	density_double = 0;
	density_doubleArray = vtkDoubleArray::New();
	featureAngle = 60;
	arrowScale = 10;
	reductionfactor = 0.99;
	reductionfactor_doubleArray = vtkDoubleArray::New();
	outlineData = vtkOutlineFilter::New();
	mapOutline = vtkPolyDataMapper::New();
	outline = vtkActor::New();
	iter = 2000;
	converge = 0.001;
	relax = 0.001;
	densID = 0;
	opacity = 0.015;
	//normalcenters_array = vtkDoubleArray::New();
	//normalcenters_array->SetNumberOfComponents(3);
	VoxVolume = vtkVolume::New();
	actorArrow = vtkSmartPointer<vtkActor>::New();

	mediumIDs = vtkDoubleArray::New(); // refraction indices (inside, outside)
	Cellcenter = vtkCellCenters::New();




}

VoxToMesh::~VoxToMesh() {




}

void VoxToMesh::convertToMesh(){


	outlineData->SetInputData(VoxImage);


	mapOutline->SetInputData(outlineData->GetOutput());



	outline->SetMapper(mapOutline);
	outline->GetProperty()->SetColor(0,0,0);


	//vtkAppendPolyData* allInOnePolyImage = vtkAppendPolyData::New();



	   contour->SetInputData(this->VoxImage);



   for(int i = 0; i<this->density_doubleArray->GetSize();i++){


		   contour->SetValue(i,this->density_doubleArray->GetValue(i));



   }



   vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();


		   vtkSmoothPolyDataFilter* sp = vtkSmoothPolyDataFilter::New();
		  	   sp->SetInputConnection(contour->GetOutputPort());
		  	   sp->SetNumberOfIterations(this->iter);
		  	   sp->SetConvergence(this->converge);
		  	   sp->SetRelaxationFactor(this->relax);
		  	   sp->SetFeatureEdgeSmoothing(0);
		  	   sp->SetBoundarySmoothing(0);
		  	   sp->Update();

		  	   cout << "Number Of Layer:" << contour->GetNumberOfContours() << endl;

		  	   cout << "Iterations: " << this->iter << endl;
		  	   cout << "COnvF: " << this->converge << endl;
		  	   cout << "RelaxF: " << this->relax << endl;

//
		  vtkDecimatePro* decimate = vtkDecimatePro::New();

		  cout << "RedFactor:" << this->reductionfactor_doubleArray->GetValue(GetDensID()) << endl;


		   decimate->SetInputConnection(sp->GetOutputPort());
		   decimate->SetTargetReduction(this->reductionfactor_doubleArray->GetValue(GetDensID()));
		   decimate->SplittingOn();
		   // decimate->BoundaryVertexDeletionOn();
		   decimate->SetFeatureAngle(this->featureAngle);
		   decimate->SetPreserveTopology(1);


		   decimate->Update();




		   cout << "CellDim: " << decimate->GetOutput()->GetCell(0)->GetCellDimension() << endl;

	   PolyImage->DeepCopy(decimate->GetOutput());
	   // PolyImage = allInOnePolyImage->GetOutput();

//	   cout << "PolyImage: NumberOfCells: " << PolyImage->GetNumberOfCells() << endl;
//	   cout << "PolyImage - NumberOfPoints:" << PolyImage->GetNumberOfPoints() << endl;

	   normals->SetInputConnection(decimate->GetOutputPort());
	   normals->ComputeCellNormalsOn();
	   normals->ComputePointNormalsOff();
	   normals->AutoOrientNormalsOn();
	   normals->ConsistencyOn();
	   normals->SplittingOff();
	   normals->SetFeatureAngle(this->featureAngle);

	   normals->Update();

//	   normals->Print(cout);




	   stripper->SetInputConnection(normals->GetOutputPort());
	   stripper->Update();

	   mapper->SetInputConnection(stripper->GetOutputPort());
	//   mapper->SetScalarRange(0.,255.);
	   mapper->Update();
	   actor->SetMapper(mapper);

	   SetCellcenters();

	   GenMediumIDs(0);


}


int VoxToMesh::GetNumberOfContours(){

	return contour->GetNumberOfContours();

}

void VoxToMesh::GetContourValues(double* values){

	contour->GetValues(values);

}



double* VoxToMesh::GetContourValues(){

	return contour->GetValues();


}




void VoxToMesh::GenMediumIDs(int mode){

	if(mode==0){
		GenMediumIDs0();
	}else{
		GenMediumIDs1();
	}

}
void VoxToMesh::GenMediumIDs1(){


	cout << "in GenMediumIDs | ";

	mediumIDs->SetNumberOfComponents(4);

	// vtkPoints* baryPoints = vtkPoints::New();
	vtkDoubleArray* baryPoints = vtkDoubleArray::New();
	baryPoints->SetNumberOfComponents(3);
	baryPoints->SetNumberOfTuples(55);

	vtkPoints* mediumPoints = vtkPoints::New();
	vtkPoints* mediumPointsTemp = vtkPoints::New();
	vtkPoints* cellpoints = vtkPoints::New();
	cellpoints->SetNumberOfPoints(3);
	vtkIdList* ptsList = vtkIdList::New();
	vtkPoints* cartesian = vtkPoints::New();

	double controlBaryPts[3];

	vtkIdType count = 0;
			for(int k=0;k<10;k++){

				double p1 = (double)k/10 + 0.05;

				for(int l=0;l<10-k;l++){

					double p2 = (double)l/10 + 0.05;
					double p3 = 1-p1-p2;

					if(p3<0.05){p3 = 0;}
					// baryPoints->SetNumberOfTuples(count+1);
					baryPoints->SetTuple3(count,p1,p2,p3);

			//		cout << " Set Barycentric Points: (count,p1,p2,p3): (" << count << ", " << baryPoints->GetTuple(count)[0] << ", " << baryPoints->GetTuple(count)[1] << ", " << baryPoints->GetTuple(count)[2] << ")" << endl;


					count++;
				}
			}


			baryPoints->GetTuple(0,controlBaryPts);
	int numBaryPoints = baryPoints->GetNumberOfTuples();


	int numCells = (int)(GetMeshImage()->GetNumberOfCells());

//	cout << "count: " << count << " | Number Of Points: " << numBaryPoints << " | Number Of Cells: " << numCells << endl;


//	cout << "BarPts:_(" << controlBaryPts[0] << ", " << controlBaryPts[1] << ", " << controlBaryPts[2] << ")" << endl;

	for(int i = 0;i<numCells;i++){

					//	cout << "Loop i: " << i << "|";

		int dim = GetMeshImage()->GetCell(i)->GetCellDimension();
		//cell->DeepCopy(this->PolyImage->GetCell(i));

		// int numpo = this->GetMeshImage()->GetCell(i)->GetNumberOfPoints();

		GetMeshImage()->GetCell(i)->GetPoints()->GetPoints(GetMeshImage()->GetCell(i)->GetPointIds(),cellpoints);

		double* n = this->GetNormals()->GetTuple3(i);

		double n_abs = sqrt(pow(n[0],2)+pow(n[1],2)+pow(n[2],2));
		n[0]/=n_abs;
		n[1]/=n_abs;
		n[2]/=n_abs;

		 int cellnumpo = cellpoints->GetNumberOfPoints();

		 	 	 	 //	 cout << "Number of CellPoints: " << cellnumpo << endl;
		 double cps[3];
		 for(int cp = 0; cp<cellnumpo;cp++){
			GetMeshImage()->GetCell(i)->GetPoints()->GetPoint(cp,cps);
			 //cellpoints->GetPoint(cp,cps);

					//	cout << "inLoop: cp: " << cp << " || ";
					//	cout << "CellPoints: P" << cp << "(";
					//	cout << cps[0] << ", " << cps[1] << ", " << cps[2] << ")" << endl;
		}

		 	 	 	// 	 cout << "|| Normalvector: n(" << n[0] << ", " << n[1] << ", " << n[2] << ")" << endl;



		BarycentricToCartesian(GetMeshImage()->GetCell(i)->GetPoints(),baryPoints,cartesian);
		MediumPoints(cartesian,n, mediumPointsTemp);
		AddPointsToPoints(mediumPointsTemp,mediumPoints);



	}

						cout << "after foorLoop: MediumPoints" << endl;


		int medSize = mediumPoints->GetNumberOfPoints();
		mediumIDs->SetNumberOfTuples(medSize);

	for(int i=0;i < medSize;i++){

			double p[3];
			mediumPoints->GetPoint(i,p);

			/*
			 * if(nextMedium !=currentMedium) n*1; else n*2;
			*/

			int p1,p2,p3;
			p1 = (int)p[0];
			p2 = (int)p[1];
			p3 = (int)p[2];

			double medium = GetVoxImage()->GetScalarComponentAsDouble(p1,p2,p3,0); // component 0, because grey scale image

			cout << "p" << i << ": " << p1 << ", " << p2 << ", " << p3 << ", " << medium << endl;


			mediumIDs->SetTuple4(i,p[0],p[1],p[2],medium);

	}






}

void VoxToMesh::GenMediumIDs0(){
	// nur fürlayermodell: 1 point/cell
	// cout << "in GenMediumIDs | ";

		mediumIDs->SetNumberOfComponents(2);
		mediumIDs->SetNumberOfTuples(GetMeshImage()->GetNumberOfCells());



		int numCells = (int)(GetMeshImage()->GetNumberOfCells());


		for(int i = 0;i<numCells;i++){


			double centers[3];
			GetCellcenters()->GetOutput()->GetPoints()->GetPoint(i,centers);
			//this->normalcenters_array->GetTuple(i,centers);

			double* n = GetNormals()->GetTuple3(i);

			double n_abs = sqrt(pow(n[0],2)+pow(n[1],2)+pow(n[2],2));
			n[0]/=n_abs;
			n[1]/=n_abs;
			n[2]/=n_abs;

			double p[3];

			p[0] = centers[0]+n[0];
			p[1] = centers[1]+n[1];
			p[2] = centers[2]+n[2];


			double mediumOuter = GetVoxImage()->GetScalarComponentAsDouble(p[0],p[1],p[2],0); // component 0, because grey scale image

			p[0] = centers[0]-n[0];
			p[1] = centers[1]-n[1];
			p[2] = centers[2]-n[2];

			double mediumInner = GetVoxImage()->GetScalarComponentAsDouble(p[0],p[1],p[2],0); // component 0, because grey scale image


		//	cout << "p" << i << ": (" << p[0] << ", " << p[1] << ", " << p[2] << "): " << medium << endl;


			mediumIDs->SetTuple2(i,mediumInner,mediumOuter);

		}


}

void VoxToMesh::AddPointsToPoints(vtkPoints* inputToAdd, vtkPoints* output){

	// cout << "in AddPoints" << endl;

	vtkPoints* out = vtkPoints::New();

	out->DeepCopy(output);

	int numPointsOutput = out->GetNumberOfPoints();

	cout << "numPointsOutput: " << numPointsOutput << " | ";

	int numPointsInput = inputToAdd->GetNumberOfPoints();

	 cout << "numPointsInput: " << numPointsInput << endl;

	 out->SetNumberOfPoints(numPointsInput+numPointsOutput);
	//output->SetNumberOfPoints(numPointsInput+numPointsOutput);



	for(int i=0;i<numPointsInput;i++){

		double p[3];
		inputToAdd->GetPoint(i,p);



		out->SetPoint(i+numPointsOutput,p);

		// cout << "P:" << (i+numPointsOutput) << ": (" << out->GetPoint(i+numPointsOutput)[0] << ", " << out->GetPoint(i+numPointsOutput)[1] << ", " << out->GetPoint(i+numPointsOutput)[2] << ")" << endl;


	}

	output->DeepCopy(out);

	int s = output->GetNumberOfPoints();

	for(int i=0;i<s;i++){


			cout << "P:" << i << ": (" << output->GetPoint(i)[0] << ", " << output->GetPoint(i)[1] << ", " << output->GetPoint(i)[2] << ")" << endl;


		}


}

void VoxToMesh::BarycentricToCartesian(vtkPoints* cellPoints, vtkDoubleArray* barypts, vtkPoints* cartesian){
/* convert baryzentric --> cartesian
 *
 * Baryzentrische Koordinatensystem der CellPoints
 *
 * c = c1(x1,y1,0)*e_x(x,y,z) + c1(x1,y1,0)*e_y(x,y,z);
 *
 *
 */

	vtkDoubleArray* barPts = vtkDoubleArray::New();
	barPts->DeepCopy(barypts);
	double bpts[3];
	barPts->GetTuple(0,bpts);
//	cout << "in BaryToCart" << endl;
//	cout << "Barycentric Points: (0,p1,p2,p3): (" << 0 << ", " << bpts[0] << ", " << bpts[1] << ", " << bpts[2] << ")" << endl;
	int bSize = barypts->GetNumberOfTuples();
	double lambda[3];

	double point0[3];
	double point1[3];
	double point2[3];


	cellPoints->GetPoint(0,point0);
	cellPoints->GetPoint(1,point1);
	cellPoints->GetPoint(2,point2);



	/*
	 * calc
	 * e_x = vec(point0,poin1)
	 * e_z = cross(e_x,vec(point0,point2))
	 * e_y = cross(e_x,e_z)
	 */

	// e_x
	double e_x[3];
	Vec(point0,point1,e_x);

	Normalize(e_x);
	// vec(point0,point2)
	double vec1[3];
	Vec(point0,point2,vec1);
	double e_z[3];
	double e_y[3];

	Cross(e_x,vec1,e_z);
	Cross(e_x,e_z,e_y);
/*
	cout << " | point0: (" << point0[0] << ", " << point0[1] << ", " << point0[2] << ")";
	cout << " | point1: (" << point1[0] << ", " << point1[1] << ", " << point1[2] << ")";
	cout << " | point2: (" << point2[0] << ", " << point2[1] << ", " << point2[2] << ")" << endl;
	cout << " | vec: (" << vec1[0] << ", " << vec1[1] << ", " << vec1[2] << ")";
	cout << " | e_x: (" << e_x[0] << ", " << e_x[1] << ", " << e_x[2] << ")";
	cout << " | e_y: (" << e_y[0] << ", " << e_y[1] << ", " << e_y[2] << ")";
	cout << " | e_z: (" << e_z[0] << ", " << e_z[1] << ", " << e_z[2] << ")" << endl;
*/
	double planeC0[2];
	vtkDoubleArray* planeCoords = vtkDoubleArray::New();
	//vtkDoubleArray* planeCoords2 = vtkDoubleArray::New();

	planeCoords->SetNumberOfComponents(2);
	planeCoords->SetNumberOfTuples(2);

	// PointInPlaneCoords(point0,point0,e_x,e_y,planeC0); =0
	PointInPlaneCoords(point0,point1,e_x,e_y,0,planeCoords);
	PointInPlaneCoords(point0,point2,e_x,e_y,1,planeCoords);
/*

	cout << " | planeP1: (" << planeCoords->GetTuple(0)[0] << ", " << planeCoords->GetTuple(0)[1] << ")";
	cout << " | planeP2: (" << planeCoords->GetTuple(1)[0] << ", " << planeCoords->GetTuple(1)[1] << ")" << endl;
*/
	vtkPoints* cartTemp = vtkPoints::New();
	cartTemp->SetNumberOfPoints(bSize);

//	cout << "Barycentric coords: " << endl;
	for(int i=0;i<bSize;i++){

		barypts->GetTuple(i,lambda);

	//	cout << "lamda: (" << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << ")";

		double x = lambda[1]*(planeCoords->GetTuple(0)[0])+lambda[2]*(planeCoords->GetTuple(1)[0]); // lambda[0]*planeC0[0] = 0
		double y = lambda[2]*(planeCoords->GetTuple(1)[1]); // lambda[0]*planeC0[1] = 0 ; planeC1[1] = 0 , da auf e_x

	//	cout << " | (x, y) : (" << x << ", " << y << ")" << endl;
		double coords[3];

		PlaneCoordsTo3D(x,y,point0,e_x,e_y,coords);

		cartTemp->SetPoint(i,coords);


	}

	//cout << "before Cartesian-SetPoint" << endl;
	cartesian->DeepCopy(cartTemp);
	//cout << "after Cartesian-SetPoint" << endl;

}

void VoxToMesh::Cross(double* a, double* b, double* c){



	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

	//Einheitsvector
	double d = sqrt(pow(c[0],2)+pow(c[1],2)+pow(c[2],2));

	c[0]/=d;
	c[1]/=d;
	c[2]/=d;

}

double VoxToMesh::Dot(double* a, double* b){

	double c;
	c = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
return c;
}

void VoxToMesh::Vec(double* a, double* b, double* vec){



	//for(int i=0;i<b)
	vec[0] = b[0]- a[0];
	vec[1] = b[1]- a[1];
	vec[2] = b[2]- a[2];


}


void VoxToMesh::Normalize(double* v){

	double d = sqrt(pow(v[0],2)+pow(v[1],2)+pow(v[2],2));

		v[0]/=d;
		v[1]/=d;
		v[2]/=d;

}

void VoxToMesh::PointInPlaneCoords(double* p0, double* point, double* ex, double* ey, double* planeCoords){


	double vec[3];
	double pCoords[2];
	Vec(p0,point,vec);

	pCoords[0] = Dot(vec,ex);
	pCoords[1] = Dot(vec,ey);
/*
	cout << "PointInPlaneCoords: vec: " << vec[0] << ", " << vec[1] << ", " << vec[2] << endl;
	cout << "e_x: (" << ex[0] << ", " << ex[1] << ", " << ex[2] << ") |";
	cout << "e_y: (" << ey[0] << ", " << ey[1] << ", " << ey[2] << ")" << endl;
	cout << "pCoords: (" << pCoords[0] << ", " << pCoords[1] << ")" << endl;

*/
	planeCoords = pCoords;


}

void VoxToMesh::PointInPlaneCoords(double* p0, double* point, double* ex, double* ey, int i, vtkDoubleArray* planeCoords){


	double vec[3];
	double pCoords[2];
	Vec(p0,point,vec);

	pCoords[0] = Dot(vec,ex);
	pCoords[1] = Dot(vec,ey);
/*
	cout << "PointInPlaneCoords: vec: (" << vec[0] << ", " << vec[1] << ", " << vec[2] << ") || ";
	cout << "e_x: (" << ex[0] << ", " << ex[1] << ", " << ex[2] << ") | ";
	cout << "e_y: (" << ey[0] << ", " << ey[1] << ", " << ey[2] << ") || ";
	cout << "pCoords: (" << pCoords[0] << ", " << pCoords[1] << ")" << endl;
*/

	planeCoords->SetTuple(i,pCoords);


}

void VoxToMesh::PlaneCoordsTo3D(double x, double y, double* p0, double* ex, double* ey, double* coords){

	coords[0] = p0[0]+x*ex[0]+y*ey[0];
	coords[1] = p0[1]+x*ex[1]+y*ey[1];
	coords[2] = p0[2]+x*ex[2]+y*ey[2];

}

void VoxToMesh::MediumPoints(vtkPoints* cartesian, double* normalVector, vtkPoints* mediumPoints){
 /* add new mediumPoints to existing Points
  * if(nextMedium !=currentMedium) n*1; else n*2;
*/

	//cout << "in MediumPoints | ";
	vtkPoints* cartTemp = vtkPoints::New();
	vtkPoints* medTemp = vtkPoints::New();

	cartTemp->DeepCopy(cartesian);

	int s = cartesian->GetNumberOfPoints();
//	cout << "Cartesian Number of Points: " << s << endl;

/*
	int numtup = mediumPoints->GetData()->GetNumberOfTuples();
	cout << "mediumPoints Number of Points: " << numtup << endl;
	int numberMP = mediumPoints->GetNumberOfPoints();
*/


	medTemp->SetNumberOfPoints(s);
	for(int i = 0; i < s;i++){


		double x[3];
		double c[3];
		cartTemp->GetPoint(i,c);

		x[0] = c[0] + normalVector[0];
		x[1] = c[1] + normalVector[1];
		x[2] = c[2] + normalVector[2];

	//	cout << "cartesian: c(" << c[0] << ", " << c[1] << ", " << c[2] << ") || ";
		//cout << "normalvector: n(" << normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << ")";
	//	cout << " || x: x(" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;


		medTemp->SetPoint(i,x);
		//cout << "after setpoint medTemp" << endl;

	}

	mediumPoints->DeepCopy(medTemp);
	//cout << "after DeepCopy to MediumPoints" << endl;

}

vtkPoints* VoxToMesh::MediumPoints(vtkPoints* cartesian, double* normalVector){
 /* add new mediumPoints to existing Points
  * if(nextMedium !=currentMedium) n*1; else n*2;
*/

	//cout << "in MediumPoints" << endl;
	vtkPoints* cartTemp = vtkPoints::New();
	vtkPoints* medTemp = vtkPoints::New();

	cartTemp->DeepCopy(cartesian);
	//cout << "after deep copy cartesian" << endl;


	int s = cartTemp->GetNumberOfPoints();
	//cout << "Cartesian Number of Points: " << s << endl;

	for(int i = 0; i < s;i++){


		double x[3];
		double c[3];
		cartTemp->GetPoint(i,c);
	//	cout << "cartesian: c(" << c[0] << ", " << c[1] << ", " << c[2] << ")"<< endl;
		x[0] = c[0] + normalVector[0];
		x[1] = c[1] + normalVector[1];
		x[2] = c[2] + normalVector[2];


	//	cout << "normalvector: n(" << normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << ")" << endl;
	//	cout << "x: x(" << x[0] << ", " << x[1] << ", " << x[2] << ")" << endl;


		medTemp->SetPoint(i,x);

	}

	return medTemp;


}


void VoxToMesh::SetDensID(int id){
	// temporary
	this->densID = id;
}

int VoxToMesh::GetDensID(){
	return this->densID;
}

vtkDoubleArray* VoxToMesh::GetDensityDoubleArray(){
	return this->density_doubleArray;
}

void VoxToMesh::SetIterations(int it){

	this->iter = it;


}
void VoxToMesh::SetConvergenceFactor(double conv){

	this->converge = conv;

}
void VoxToMesh::SetRelaxationFactor(double rela){

	this->relax = rela;

}

void VoxToMesh::SetVolumeOpacity(double op){
	this->opacity = op;
}

void VoxToMesh::SetVoxImage(vtkImageData* input){

	this->VoxImage->DeepCopy(input);

}

void VoxToMesh::SetDensityValue(double density) {

	this->density_double = density;


}
void VoxToMesh::SetDensityValue(vtkDoubleArray* input){
	density_doubleArray->DeepCopy(input);

}
void VoxToMesh::SetDensityValue(std::vector<float> &densities){

	int densSize = densities.size();
 	this->density_doubleArray->Resize(densSize);
	for(int i=0;i<densSize;i++){
		this->density_doubleArray->SetValue(i,double(densities[i]));

	}

}


void VoxToMesh::SetTargetReduction(double red) {

	this->reductionfactor = red;
}

void VoxToMesh::SetTargetReduction(vtkDoubleArray* input) {

	this->reductionfactor_doubleArray->DeepCopy(input);
}

void VoxToMesh::SetFeatureAngle(double angle) {
	this->featureAngle = angle;
}

void VoxToMesh::SetArrowScale(double scale) {
	this->arrowScale = scale;
}

vtkDataArray* VoxToMesh::GetNormals(){
	return this->normals->GetOutput()->GetCellData()->GetNormals();
}

vtkPolyData* VoxToMesh::GetNormalsOutput(){
	return this->normals->GetOutput();
}

vtkPolyData* VoxToMesh::GetMeshImage(){
	return this->PolyImage;
}

vtkActor* VoxToMesh::GetActor() {
	return this->actor;
}


void VoxToMesh::SetCellcenters(){


	GetCellcenters()->VertexCellsOn();
	GetCellcenters()->SetInputData(this->GetNormalsOutput());
	GetCellcenters()->Update();

//		   this->normalcenters_array->SetNumberOfTuples(GetCellcenters()->GetOutput()->GetNumberOfPoints());
//
//
//		     for(int i = 0;i<GetCellcenters()->GetOutput()->GetNumberOfPoints();i++){
//
//
//		    	  //	 vtkDataArray* n = this->GetNormals();
//		    		 //	double n[3];
//		    	//  	 	n->GetTuple(i,t);
//
//		    	double p[3];
//		    	GetCellcenters()->GetOutput()->GetPoints()->GetPoint(i,p);
//
//		  	   this->normalcenters_array->SetTupleValue(i,p);
//
//		  	   //this->normalcenters_array->GetTupleValue(i,n);
//
//		  	 }


}

vtkCellCenters* VoxToMesh::GetCellcenters(){

	return this->Cellcenter;

}

void VoxToMesh::GetMediumIDs(double** medID){

	int t = this->mediumIDs->GetNumberOfTuples();
	int comp = this->mediumIDs->GetNumberOfComponents();
	double tup[comp];


	medID = new double*[t];

	for(int i=0;i<t;i++){

		medID[i] = new double[comp];

		this->mediumIDs->GetTuple(i,tup);

		for(int j=0;j<comp;j++){

			medID[i][j] = tup[j];

		}

	}


}
vtkDoubleArray* VoxToMesh::GetMediumIDs(){
	return this->mediumIDs;

}

vtkSmartPointer<vtkActor> VoxToMesh::GetNormalArrowsActor(){


//	 vtkCellCenters *Cellcenter = vtkCellCenters::New();
//	 Cellcenter->VertexCellsOn();
//	 Cellcenter->SetInputData(this->GetNormalsOutput());
//	 Cellcenter->Update();
//
//	   this->normalcenters_array->SetNumberOfTuples(Cellcenter->GetOutput()->GetNumberOfPoints());
//
//
//	     for(int i = 0;i<Cellcenter->GetOutput()->GetNumberOfPoints();i++){
//
//
//	    	  //	 vtkDataArray* n = this->GetNormals();
//	    		 //	double n[3];
//	    	//  	 	n->GetTuple(i,t);
//
//	    	double p[3];
//	  	   Cellcenter->GetOutput()->GetPoints()->GetPoint(i,p);
//
//	  	   this->normalcenters_array->SetTupleValue(i,p);
//
//	  	   //this->normalcenters_array->GetTupleValue(i,n);
//
//	  	 }

	   vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
	   arrowSource->SetTipRadius(0.1);

	    vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();

	    // glyph3D->SetInputData(dummy_center->GetOutput());
	    glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
	    glyph3D->SetVectorModeToUseNormal();
	    // glyph3D->SetInputData(vm->GetMeshImage());
	    glyph3D->SetInputData(GetCellcenters()->GetOutput());
	    glyph3D->SetScaleFactor(arrowScale);// scalefactor
	    glyph3D->Update();


	    vtkSmartPointer<vtkPolyDataMapper> mapperArrow = vtkSmartPointer<vtkPolyDataMapper>::New();
	     mapperArrow->SetInputConnection(glyph3D->GetOutputPort());
	     // mapperArrow->AddInputDataObject(center->GetOutput());

	     this->actorArrow->SetMapper(mapperArrow);


return this->actorArrow;
}

vtkVolume* VoxToMesh::GetVoxelVolume(){

	 vtkSmartPointer<vtkVolumeProperty> volprop = vtkSmartPointer<vtkVolumeProperty>::New();


	   vtkPiecewiseFunction *opacity = vtkPiecewiseFunction::New();

	   for(int i=0;i<this->GetDensityDoubleArray()->GetSize();i++){
		   opacity->AddPoint(this->GetDensityDoubleArray()->GetValue(i), this->opacity);
	   }

	   volprop->SetScalarOpacity(opacity);


	vtkSmartVolumeMapper *VoxMapper = vtkSmartVolumeMapper::New();
	   VoxMapper->SetBlendModeToComposite();
	   VoxMapper->SetInputData(this->GetVoxImage());



	   VoxVolume->SetMapper(VoxMapper);
	   VoxVolume->SetProperty(volprop);
	   VoxVolume->Update();


return VoxVolume;
}

vtkDoubleArray* VoxToMesh::GetNormalCenters(){

	vtkDoubleArray* centers = vtkDoubleArray::New();

	centers->SetNumberOfComponents(3);
	centers->SetNumberOfTuples(GetCellcenters()->GetOutput()->GetNumberOfPoints());

	for(int i = 0;i<GetCellcenters()->GetOutput()->GetNumberOfPoints();i++){

		double p[3];
		GetCellcenters()->GetOutput()->GetPoints()->GetPoint(i,p);

		centers->SetTupleValue(i,p);

	 }

	return centers;
}

vtkSmartPointer<vtkImageData> VoxToMesh::GetVoxImage(){

	return this->VoxImage;

}

vtkSmartPointer<vtkImageData> VoxToMesh::GetIDVoxelImage(){

	return this->IDVoxelImage;

}



void VoxToMesh::PrintNormalsCenters(){

cout << "Number of Components: " << GetCellcenters()->GetOutput()->GetPointData()->GetNumberOfComponents() << endl;
cout << "Number of Normals: " << GetCellcenters()->GetOutput()->GetPoints()->GetNumberOfPoints() << endl;
cout << "Number of MediumIDs-Tuples: " << this->mediumIDs->GetNumberOfTuples() << endl;
cout << "Number of MediumIDs-Components: " << this->mediumIDs->GetNumberOfComponents() << endl;

	for(int i = 0;i<GetCellcenters()->GetOutput()->GetPoints()->GetNumberOfPoints();i++){

		double p[3];
		GetCellcenters()->GetOutput()->GetPoints()->GetPoint(i,p);
		double n[3];
		GetNormals()->GetTuple(i,n);

/*
		cout << " Normalvector " << i << ": ( " << n[0] << ", " << n[1] << ", " << n[2] << ")";
		cout << " @P(" << p[0] << ", " << p[1] << ", " << p[2] << ")" << endl;*/

    }




	for(int i=0;i<this->mediumIDs->GetNumberOfTuples();i++){

		double* m = this->mediumIDs->GetTuple2(i);
		double p0[3];
		double p1[3];
		double p2[3];
		this->PolyImage->GetCell(i)->GetPoints()->GetPoint(0,p0);
		this->PolyImage->GetCell(i)->GetPoints()->GetPoint(1,p1);
		this->PolyImage->GetCell(i)->GetPoints()->GetPoint(2,p2);
		cout << " Medium (inside, outside): " << i << "(" << m[0] << ", " << m[1] << ") ";
		cout << "P0(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ") | ";
		cout << "P1(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ") | ";
		cout << "P2(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ")" << endl;


	}


}





