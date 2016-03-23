/*
 * createLayerModel.cpp
 *
 *  Created on: Mar 3, 2015
 *      Author: matthias
 */

#include "createModel.hpp"
#include <cmath>
#include "math.h"
// #include "vtkArray.h"
// #include "vtkObjectFactory.h"
// #include "vtkObject.h"
#include "vtkTransform.h"
#include "vtkImageReslice.h"
#include "vtkPolyData.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkPolyDataNormals.h"
#include "vtkTriangle.h"
#include <vtkPointData.h>
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"

#include "vtkPlaneSource.h"
#include "vtkAppendPolyData.h"
#include "vtkProperty.h"
#include "vtkUnsignedCharArray.h"
#include "vtkArrowSource.h"
#include "vtkGlyph3D.h"


// #include "vtkImageData.h"
// #include "vtkImageCast.h"

// vtkStandardNewMacro(createModel);

createModel::createModel(){

	dimensions_array = vtkIntArray::New();
	densities_array = vtkDoubleArray::New();
	layerdepths_array = vtkIntArray::New();
	rotation_array  = vtkDoubleArray::New();
	modelID_int = 0;
	modeID_int = 0;
//	filename = new char;
	//reader = new ReadXML;
	arrowScale = 10;
	PolyImage = vtkPolyData::New();
	mediumIDs = vtkDoubleArray::New();
	normals = vtkPolyDataNormals::New();
	stripper = vtkStripper::New();
	mapper = vtkPolyDataMapper::New();
	actor = vtkActor::New();

	image = vtkImageData::New();
	Cellcenter = vtkCellCenters::New();
	actorArrow = vtkSmartPointer<vtkActor>::New();
	// this->densities_array->InsertNextValue(0);

//	this->layerdepths_array->InsertNextValue(0);

	//Size of LayerImage
	this->dimensions_array->Resize(3);
	this->dimensions_array->SetValue(0,128); //width
	this->dimensions_array->SetValue(1,128); // height
	this->dimensions_array->SetValue(2,128); // depth

	this->rotation_array->InsertNextValue(0); // yaw
	this->rotation_array->InsertNextValue(0); // pitch
	this->rotation_array->InsertNextValue(0); //roll


}


void createModel::setModel(int model) {

	modelID_int = model;

}

/*
void createModel::setFilename(char* filename){

//	this->filename = filename;
}*/

/*
void createModel::setInputVarsLayer(){


	this->reader->SetFilename(filename);
	this->reader->ParseFile();
	// reader->PrintVars();

	this->setModel(reader->GetModelID());
	this->setSize(reader->GetDimensions());
	this->setInputVarsLayer(reader->GetDensities(),reader->GetLayerdepth(),reader->GetRotations());
	this->setMode(reader->GetModeID());


}

vtkDoubleArray* createModel::GetDensities(){


	return this->reader->GetDensities();


}
double createModel::GetDensities(int densID){


	return this->reader->GetDensities(densID);


}

*/

void createModel::setInputVarsLayer(vtkDoubleArray* input_densities,vtkIntArray* input_layerdepts, vtkDoubleArray* input_rotations){

	// setting Density-Array
	unsigned int densSize = sizeof(input_densities);
	this->setDensities(input_densities,densSize);

	// setting layerdepth-Array
	unsigned int layerdepthSize = sizeof(input_layerdepts);
	this->setLayerDepth(input_layerdepts,layerdepthSize);

	// setting rotation
		this->setRotation(input_rotations);


}

void createModel::setInputVarsLayer(std::vector<float> &densities,std::vector<int> &layerdepths, std::vector<float> &rotations){

	int densSize = densities.size();
	int ldSize = layerdepths.size();
	int rotSize = rotations.size();

	this->densities_array->Resize(densSize);
	this->layerdepths_array->Resize(ldSize);
	this->rotation_array->Resize(rotSize);

	for(int i=0;i<densSize;i++){this->densities_array->SetValue(i,double(densities[i]));}
	for(int i=0;i<ldSize;i++){this->layerdepths_array->SetValue(i,double(layerdepths[i]));}
	for(int i=0;i<rotSize;i++){this->rotation_array->SetValue(i,double(rotations[i]));}

}


void createModel::setDensities(vtkDoubleArray* densities, unsigned size) {

	// this->dense->Resize(size);
	// this->densities_array->SetArray(densities,size,1);
	this->densities_array = densities;

}

void createModel::setLayerDepth(vtkIntArray* layerdepths, unsigned size) {

	// this->ldep->Resize(size);
//	this->ldep->SetArray(layerdepths,size,1);
	this->layerdepths_array = layerdepths;
}

void createModel::setSize(vtkIntArray* dims){

	this->dimensions_array = dims;
}

void createModel::setSize(int width, int height, int depth){
	this->dimensions_array->SetValue(0,width);
	this->dimensions_array->SetValue(1,height);
	this->dimensions_array->SetValue(2,depth);
}
void createModel::setMode(int mode){

	modeID_int = mode;
}

//
//int createModel::getLayer(int x, int y, int z){
//
//	int layer = 0; // equals Air
//	int i = 0; // Idx of Array ldep
//
//
//
//	int s = this->layerdepths_array->GetElementComponentSize();
//	double z_new = z-y*tan(this->rotation_array->GetValue(0)*M_PI/180.0);
//	double sum = 0;
//
//
//	if(z_new>=0) {
//
//
//
//		do{
//				if (i<s){
//
//
//					sum += (this->layerdepths_array->GetValue(i))/cos(this->rotation_array->GetValue(0)*M_PI/180.0);
//					layer = i+1;
//
//				}else {
//					layer = 0;
//					break;
//				}
//				i++;
//
//
//			}while(sum<=z_new);
//
//	} else {
//
//		layer = 0;
//	}
//
//	return layer;
//}

int createModel::getLayer(int x, int y, int z){

	int layer = 0; // equals Air
	int i = 0; // Idx of Array layerdepths_array



	int s = this->layerdepths_array->GetSize();
	// cout << "LD: compSize: " << s << endl;
	// double z_new = z;
	double sum = 0;


	if(z>=0) {



		do{
				if (i<s){


					sum += (this->layerdepths_array->GetValue(i));
					layer = i+1;

				}else {
					layer = 0;
					break;
				}
				i++;


			}while(sum<=z);

	} else {

		layer = 0;
	}

	return layer;
}

unsigned int createModel::getLayerDensity(int x, int y, int z){

	int density = 0; // Density of Air
	int layerID = this->getLayer(x,y,z);

	if (layerID >0){
		// layerID == 0 --> density of Air

		density = this->densities_array->GetValue(layerID-1);
	}



	return density;
}

void createModel::setRotation(vtkDoubleArray* rotation){

	this->rotation_array = rotation;

}


void createModel::setRotation(double rot_yaw, double rot_pitch, double rot_roll){

	this->rotation_array->SetValue(0,rot_yaw); //heading
	this->rotation_array->SetValue(1,rot_pitch); // elevation
	this->rotation_array->SetValue(2,rot_roll); //bank

}


void createModel::SetArrowScale(double scale) {
	this->arrowScale = scale;
}

//void createModel::rotateImage(){
//
//	double* bounds = this->image->GetBounds();
//
//	 // Compute the center of the image
//	  double center[3];
//	  center[0] = (bounds[1] + bounds[0]) / 2.0;
//	  center[1] = (bounds[3] + bounds[2]) / 2.0;
//	  center[2] = (bounds[5] + bounds[4]) / 2.0;
//
//	  cout << center[0] << " " << center[1] << " " << center[2] << endl;
//
//	 vtkTransform *transform = vtkTransform::New();
//	 transform->Scale(2,2,1);
//
//	 transform->Translate(2*center[0],2*center[1], center[2]);
//	 transform->RotateX(this->rotation_array->GetValue(2));
//	 transform->RotateY(this->rotation_array->GetValue(1));
//	 transform->RotateZ(this->rotation_array->GetValue(0));
//	 transform->Translate(-2*center[0], -2*center[1], -center[2]);
//	 transform->Scale(0.5,0.5,1);
//
//	  vtkImageReslice *reslice = vtkImageReslice::New();
//
//	  reslice->SetInputDataObject(this->image);
//
//	  // reslice->SetInputConnection(image->GetOutputPort());
//	  reslice->SetResliceTransform(transform);
//	  reslice->SetInterpolationModeToCubic();
//
//
//	  reslice->SetOutputSpacing(image->GetSpacing()[0], image->GetSpacing()[1], image->GetSpacing()[2]);
//	  reslice->SetOutputOrigin(image->GetOrigin()[0],image->GetOrigin()[1],image->GetOrigin()[2]);
//	  reslice->SetOutputExtent(image->GetExtent());
//
//	  reslice->Update();
//
//	 image->DeepCopy(reslice->GetOutput());
//
//
////
//
////	  reslice->Delete();
////	  transform->Delete();
////	  tempImage->Delete();
//
//
//
//}

void createModel::PrintVariable() {


		cout << this->dimensions_array->GetValue(0) << endl;
		cout << this->dimensions_array->GetValue(1) << endl;
		cout << this->dimensions_array->GetValue(2) << endl;

		cout << "layerdepth:" << this->layerdepths_array->GetElementComponentSize() << endl;

		cout << "layerdepth:";

		for (int i=0;i<this->layerdepths_array->GetElementComponentSize();i++){
			cout << " " << this->layerdepths_array->GetValue(i);
		}
		cout << endl;

		cout << "dense:" << this->densities_array->GetElementComponentSize() << endl;

}

void createModel::Update(){

//	this->setInputVarsLayer();

	int width = this->dimensions_array->GetValue(0);
	int height = this->dimensions_array->GetValue(1);
	int depth = this->dimensions_array->GetValue(2);

	this->image->SetDimensions(width,height,depth);
	this->image->AllocateScalars(VTK_UNSIGNED_CHAR,1);


	switch (this->modelID_int) {
		case 0:

			if(this->modeID_int == 0){
				this->CreateLayerModel(width,height, depth);
			}else{
				this->CreateLayerModel();
			}

			break;
		case 1:
			this->CreateSphereModel(width,height, depth);
			break;
		case 2:
			this->Create3DModel(width,height, depth);
			break;
		default:
			this->CreateLayerModel(width,height, depth);
	}

//this->reader->PrintVars();
}


void createModel::CreateLayerModel(int width, int height,int depth){

	unsigned char value = 0;

	union
	{
	unsigned int word;
	unsigned char bytes[2];
	}density;


		for (int k=0; k<depth; k++) {

			for(int j = 0;j<height;j++) {

				for (int i=0;i<width;i++) {

					unsigned char* pixel = static_cast<unsigned char*>(image->GetScalarPointer(i,j,k));


					density.word = this->getLayerDensity(i,j,k);

					pixel[0] = density.bytes[0];


				}
			}

		}


	//this->rotateImage();


}


/*
void createModel::GenMediumIDs(){
	// nur fürlayermodell: 1 point/cell
	// cout << "in GenMediumIDs | ";

		mediumIDs->SetNumberOfComponents(2);
		mediumIDs->SetNumberOfTuples(PolyImage->GetNumberOfCells());



		int numCells = (int)(PolyImage->GetNumberOfCells());


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
*/
void createModel::CreateLayerModel(){
// create Mesh directly


	double xdim = this->dimensions_array->GetValue(0);
	double ydim = this->dimensions_array->GetValue(1);
	double zcoord = 0.0;
	int numLayers = this->layerdepths_array->GetSize();

	for(int i=0;i<numLayers;i++){zcoord+=this->layerdepths_array->GetValue(i);}

	double m[3] = {xdim/2,ydim/2,(zcoord-(this->layerdepths_array->GetValue(0))/2)};
	cout << "Number of Layers: " << numLayers << endl;
	mediumIDs->SetNumberOfComponents(2);
	mediumIDs->SetNumberOfTuples(2*numLayers);

	 // Create a triangle
	  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	  vtkSmartPointer<vtkCellArray> Vertices = vtkSmartPointer<vtkCellArray>::New();
	  vtkSmartPointer<vtkPolyData> pImage = vtkSmartPointer<vtkPolyData>::New();



	  zcoord = 0.0;
	  //setup points (geometry)
	  for(int i=0;i<numLayers;i++){
		  zcoord+=this->layerdepths_array->GetValue(i);

		  std::vector<double> p0;p0.assign(3,0); p0[0] = 0.0; p0[1] = 0.0; p0[2] = zcoord;
		  std::vector<double> p1;p1.assign(3,0); p1[0] = 0.0; p1[1] = ydim; p1[2] = zcoord;
		  std::vector<double> p2;p2.assign(3,0); p2[0] = xdim; p2[1] = 0.0; p2[2] = zcoord;
		  std::vector<double> p3;p3.assign(3,0); p3[0] = xdim; p3[1] = ydim; p3[2] = zcoord;
		  std::vector<double> p0_Rot,p1_Rot,p2_Rot,p3_Rot;

		  Rotate(p0,p0_Rot,m);
		  Rotate(p1,p1_Rot,m);
		  Rotate(p2,p2_Rot,m);
		  Rotate(p3,p3_Rot,m);

	   points->InsertNextPoint ( p0_Rot[0], p0_Rot[1], p0_Rot[2]);
	   points->InsertNextPoint ( p1_Rot[0], p1_Rot[1], p1_Rot[2]);
	   points->InsertNextPoint ( p2_Rot[0], p2_Rot[1], p2_Rot[2]);
	   points->InsertNextPoint ( p3_Rot[0], p3_Rot[1], p3_Rot[2]);


	   vtkSmartPointer<vtkTriangle> triangle0 = vtkSmartPointer<vtkTriangle>::New();
	   vtkSmartPointer<vtkTriangle> triangle1 = vtkSmartPointer<vtkTriangle>::New();

	   triangle0->GetPointIds()->SetId ( 0, 4*i );
	   triangle0->GetPointIds()->SetId ( 1, 4*i+1 );
	   triangle0->GetPointIds()->SetId ( 2, 4*i+2 );

	   triangle1->GetPointIds()->SetId ( 0, 4*i+1 );
	   triangle1->GetPointIds()->SetId ( 1, 4*i+3 );
	   triangle1->GetPointIds()->SetId ( 2, 4*i+2 );


	   Vertices->InsertNextCell ( triangle0 );
	   Vertices->InsertNextCell ( triangle1 );

	   double mediumOuter;
	   double mediumInner;

	   i==numLayers-1 ? mediumInner=0 : mediumInner=i+1;
	   mediumOuter=i;

	   mediumIDs->SetTuple2(2*i,mediumInner,mediumOuter);
	   mediumIDs->SetTuple2(2*i+1,mediumInner,mediumOuter);



	  }


		PolyImage->SetPoints(points);
		PolyImage->SetPolys(Vertices);

		 vtkSmartPointer<vtkUnsignedCharArray> colors =vtkSmartPointer<vtkUnsignedCharArray>::New();
		  colors->SetNumberOfComponents(3);
		  colors->SetName("Colors");
		  // Define some colors
		    unsigned char red[3] = {255, 0, 0};
		    unsigned char green[3] = {0, 255, 0};
		    unsigned char blue[3] = {0, 0, 255};
		    unsigned char blue1[3] = {0, 0, 255};
		  // Add the three colors we have created to the array
		  colors->InsertNextTupleValue(red);
		  colors->InsertNextTupleValue(green);
		  colors->InsertNextTupleValue(blue);
		  colors->InsertNextTupleValue(blue1);

		 std::cout << "There are " << PolyImage->GetNumberOfCells() << " cells." << std::endl;



		normals->SetInputData(PolyImage);
		normals->ComputeCellNormalsOn();
		normals->ComputePointNormalsOff();
		normals->AutoOrientNormalsOn();
		normals->ConsistencyOn();
		normals->SplittingOn();
		//normals->SetFeatureAngle(this->featureAngle);
		normals->Update();


		 stripper->SetInputConnection(normals->GetOutputPort());
		 stripper->Update();

		 mapper->SetInputConnection(stripper->GetOutputPort());
		 mapper->GetInput()->GetCellData()->SetScalars(colors);
		 //   mapper->SetScalarRange(0.,255.);
		 mapper->Update();
		 vtkProperty* backFaces = vtkProperty::New();
		 backFaces->SetSpecular(0.0);
		 backFaces->SetDiffuse(0.0);
		 backFaces->SetAmbient(1.0);
		 backFaces->SetAmbientColor(1.0000, 0.3882, 0.2784);

		 actor->SetMapper(mapper);
		// actor->GetProperty()->SetColor(1,0,0);
		 actor->GetProperty()->SetLighting(false);
		 actor->GetProperty()->SetOpacity(0.5);
		 actor->GetProperty()->SetEdgeVisibility(true);
		 actor->GetProperty()->BackfaceCullingOff();

		 SetCellcenters();

		 PrintNormalsCenters();

}

vtkDoubleArray* createModel::GetMediumIDs(){
	return this->mediumIDs;

}


vtkPolyData* createModel::GetMeshImage(){
	return this->PolyImage;
}

void createModel::SetCellcenters(){


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


vtkSmartPointer<vtkActor> createModel::GetNormalArrowsActor(){


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


vtkPolyData* createModel::GetNormalsOutput(){
	return this->normals->GetOutput();
}
vtkCellCenters* createModel::GetCellcenters(){

	return this->Cellcenter;

}

vtkActor* createModel::GetActor() {
	return this->actor;
}


void createModel::CreateSphereModel(int width, int height,int depth){

// Not yet implemented

	this->CreateLayerModel(width,height, depth); // as precaution

}

void createModel::Create3DModel(int width, int height,int depth){

// not yet implemented

	this->CreateLayerModel(width,height, depth); // as precaution
}


vtkSmartPointer<vtkImageData> createModel::GetOutput(){

	return this->image;

}

vtkSmartPointer<vtkPolyData> createModel::GetOutputPolyImage(){

	return this->PolyImage;

}

vtkDataArray* createModel::GetNormals(){
	return this->normals->GetOutput()->GetCellData()->GetNormals();
}

void createModel::PrintNormalsCenters(){

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
		cout << "\tP0(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ")\t| ";
		cout << "\tP1(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ")\t| ";
		cout << "\tP2(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ")" << endl;


	}


/*
	int s = GetMeshImage()->GetNumberOfCells();
	for(int i=0;i<s;i++){

		int t = GetMeshImage()->GetCell(i)->GetNumberOfPoints();

		cout << "CellNr." << i << ": ";


			double p0[3];
			double p1[3];
			double p2[3];

			double vec01[3];
			double vec02[3];


			GetMeshImage()->GetCell(i)->GetPoints()->GetPoint(0,p0);
			GetMeshImage()->GetCell(i)->GetPoints()->GetPoint(1,p1);
			GetMeshImage()->GetCell(i)->GetPoints()->GetPoint(2,p2);

			Vec(p0,p1,vec01);
			Vec(p0,p2,vec02);

			cout << "P0" << "(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ") |";
			cout << "P1" << "(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ") |";
			cout << "P2" << "(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ") |";
			cout << "Vec01" << "(" << vec01[0] << ", " << vec01[1] << ", " << vec01[2] << ") |";
			cout << "Vec02" << "(" << vec02[0] << ", " << vec02[1] << ", " << vec02[2] << ") |";




		cout << "|" << endl;


	}
*/


}

void createModel::Rotate(std::vector<double> pIn, std::vector<double> &pOut, double *M){

	double cosYaw = cos(this->rotation_array->GetValue(0));
	double cosPitch = cos(this->rotation_array->GetValue(1));
	double cosRoll = cos(this->rotation_array->GetValue(2));
	double sinYaw = sqrt(1-cosYaw*cosYaw);
	double sinPitch = sqrt(1-cosPitch*cosPitch);
	double sinRoll = sqrt(1-cosRoll*cosRoll);


// Rotate:

	std::vector<double> p;
	p.assign(3,0);
	double R_13[3] = {cosYaw*cosPitch,cosYaw*sinPitch*sinRoll-sinYaw*cosRoll,cosYaw*sinPitch*cosRoll+sinYaw*sinRoll};
	double R_23[3] = {sinYaw*cosPitch, sinYaw*sinPitch*sinRoll+cosYaw*cosRoll,sinYaw*sinPitch*cosRoll-cosYaw*sinRoll};
	double R_33[3] = {-sinPitch,cosPitch*sinRoll,cosPitch*cosRoll};


	p[0] = R_13[0]*(pIn[0]-M[0])+R_13[1]*(pIn[1]-M[1])+R_13[2]*(pIn[2]-M[2])+M[0];
	p[1] = R_23[0]*(pIn[0]-M[0])+R_23[1]*(pIn[1]-M[1])+R_23[2]*(pIn[2]-M[2])+M[1];
	p[2] = R_33[0]*(pIn[0]-M[0])+R_33[1]*(pIn[1]-M[1])+R_33[2]*(pIn[2]-M[2])+M[2];

	pOut.assign(p.begin(),p.end());

}




