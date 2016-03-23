/*
 * ImportExport.cpp
 *
 *  Created on: Mar 22, 2015
 *      Author: matthias
 */

#include "ImportExport.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ReadXML.hpp"

//#include "boost/xpressive/xpressive.hpp"
#include "boost/regex.hpp"

#include "vtkXMLPolyDataReader.h"
#include "vtkXMLDataParser.h"
#include "vtkXMLWriter.h"

#include "vtkOBJReader.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkSmartPointer.h"
#include "vtkCell.h"
#include "vtkCellData.h"
#include "vtkTriangle.h"
#include "vtkProperty.h"
#include "vtkCleanPolyData.h"
#include "vtkArrowSource.h"
#include "vtkGlyph3D.h"


ImportExport::ImportExport() {

	verbose = false;

	points = vtkSmartPointer<vtkPoints>::New();
	mediumIDs = vtkDoubleArray::New();
	mediumIDs->SetNumberOfComponents(2);
	Vertices = vtkSmartPointer<vtkCellArray>::New();
	normals = vtkPolyDataNormals::New();
	ppu = 1;
	unit = "";
	ModelName = "";
	configfile = "";
	modelfile = "";
	saveDir = "";
	emitterFileName = "";
	emitterDistance = 1;
	arrowScale = 1;
	chance = 0.1;
	w_thresh = 0.0001;
	numPhotons = 10;
	numPhotonsinteractions = 20;
	importmesh = vtkPolyData::New();
	//reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	saveFile = false;
	stripper = vtkStripper::New();
	mapper = vtkPolyDataMapper::New();
	actor = vtkActor::New();

	Cellcenter = vtkCellCenters::New();
	actorArrow = vtkSmartPointer<vtkActor>::New();

	surfaceFilename = "";
	transmissionFilename = "";
	fetalFilename = "";
	recordAxzFilename = "";
	recordAyzFilename = "";

	//writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
}

ImportExport::~ImportExport() {
	// TODO Auto-generated destructor stub
	cout << "Destructing ImportExport" << endl;
//	verbose = false;

//		ppu = 1;
//		unit = NULL;
//		ModelName = NULL;
//		configfile = NULL;
//		modelfile = NULL;
//		saveDir = NULL;
//		emitterFileName = NULL;
//		arrowScale = 0;
//		chance = 0;
//		w_thresh = 0;
//		numPhotons = 0;
//		numPhotonsinteractions = 0;
//
//		saveFile = false;

		actor->Delete();
		mapper->Delete();
		stripper->Delete();
		Cellcenter->Delete();
		actorArrow->Delete();
	//	writer->Delete();
		points->Delete();
		mediumIDs->Delete();
		Vertices->Delete();
		normals->Delete();
		importmesh->Delete();

	//	reader->Delete();



}

void ImportExport::SetFilename(std::string filename){

//	reader->SetFileName(filename);

}

void ImportExport::SetConfig(std::string filename){

	configfile.assign(filename.begin(),filename.end());
}


void ImportExport::ParseFile(){
if(verbose){
	cout << "Importing data..." << endl;
}
	 ReadXML *reader = new ReadXML;
	 reader->ParseConfig(configfile);

	 reader->GetRefractionIndices(refractionInd);
	 reader->GetAbsorbtionCoefficients(absorbtionCoeffs);
	 reader->GetAnisotropies(anisotropies);
	 reader->GetScatterCoefficients(scatterCoeffs);
	 reader->GetDensities(densities);
	 reader->GetLayerdepth(layerdepths);
	 reader->GetDimensions(dims);
	 reader->GetRotations(rotations);
	 reader->GetTissueNames(tissueNames);
	 reader->GetViewRotations(view_array);
	 reader->GetSurfaceIDVector(surfaceIDs);
	 reader->GetSurfaceConfigNames(surfaceIDnamesConfig);
	 reader->GetSurfaceNamesVector(surfaceIDnames);

	 reader->GetEmitterPointInPlane(emitterPoint);
	 reader->GetDetectionGrid(detGridDims);

	 reader->GetPictureDimensions(picDims);

	 reader->GetEmitterHeading(emitterHeading);

	 chance = reader->GetChance();
	 w_thresh = reader->GetWeightThreshold();
	 numPhotons = reader->GetNumberOfPhotons();
	 numPhotonsinteractions = reader->GetNumberOfPhotonInteractions();
	 saveFile = reader->GetSaveFile();
	 saveDir = reader->GetSaveDir();
	 emitterFileName = reader->GetEmitterFile();
	 ppu = reader->GetPixelPerUnitFactor();
	 arrowScale = reader->GetArrowScaleFactor();
	 emitterDistance = reader->GetEmitterdistanceToSurface();


	 surfaceFilename = reader->GetSurfaceFilename();
	 transmissionFilename = reader->GetTransmissionFilename();
	 fetalFilename = reader->GetFetalFilename();
	 recordAxzFilename = reader->GetRecordAxzFilename();
	 recordAyzFilename = reader->GetRecordAyzFilename();

	 scaleString = reader->GetScaleString();
	 coordsysString = reader->GetCoordSystemString();
	 colormapString = reader->GetColormapString();

	 fetalID = reader->GetFetalID();

	 this->modelfile = reader->GetModelFilename();
	 this->saveDir = reader->GetSaveDir();
	 this->ImportModel(reader->GetModelFilename());


	 delete reader;

}

vtkPolyData* ImportExport::GetMeshImage(){

	return this->importmesh;
}


vtkSmartPointer<vtkActor> ImportExport::GetNormalArrowsActor(){


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
	   //  mapperArrow->AddInputDataObject(Cellcenter->GetOutput());

	     this->actorArrow->SetMapper(mapperArrow);


return this->actorArrow;
}

vtkDataArray* ImportExport::GetNormals(){
	return this->normals->GetOutput()->GetCellData()->GetNormals();
}

void ImportExport::ExportToVTP(vtkPolyData* data, std::string filename){
	// output file: name.vtp
/*
	writer->SetInputData(data);
	writer->SetFileName(filename);
	writer->Write();
	*/
}

void ImportExport::ImportModel(std::string fname){



	ImportObjModel(fname);
	ImportEmitterModel(emitterFileName);



			 stripper->SetInputConnection(normals->GetOutputPort());
			 stripper->Update();

			 mapper->SetInputConnection(stripper->GetOutputPort());
			// mapper->GetInput()->GetCellData()->SetScalars(colors);
			 //   mapper->SetScalarRange(0.,255.);
			 mapper->Update();

//			 vtkProperty* backFaces = vtkProperty::New();
//			 backFaces->SetSpecular(0.0);
//			 backFaces->SetDiffuse(0.0);
//			 backFaces->SetAmbient(1.0);
//			 backFaces->SetAmbientColor(1.0000, 0.3882, 0.2784);


			 actor->SetMapper(mapper);
			// actor->GetProperty()->SetColor(1,0,0);
			 actor->GetProperty()->SetLighting(false);
			 actor->GetProperty()->SetOpacity(0.5);
			 actor->GetProperty()->SetEdgeVisibility(true);
			 actor->GetProperty()->BackfaceCullingOff();

			 SetCellcenters();



}


void ImportExport::ImportXMLModel(std::string fname){

	vtkXMLDataParser *parser = vtkXMLDataParser::New();
	parser->SetFileName(fname.c_str());

	 if(parser->Parse()==1){

		vtkXMLDataElement *root = parser->GetRootElement(); // <model>

		int numberOfComponents = root->GetNumberOfNestedElements();


		 for(int i=0; i<numberOfComponents;i++){

			 vtkXMLDataElement *dataElement = root->GetNestedElement(i);

			 if(strcmp("mediumsList", dataElement->GetName()) == 0){
				 this->ExtractMediumsList(dataElement);


			 }else if(strcmp("name", dataElement->GetName()) == 0){
				 this->ExtractName(dataElement);

			}else if(strcmp("unit", dataElement->GetName()) == 0){
				this->ExtractUnit(dataElement);

			}else if(strcmp("ppu", dataElement->GetName()) == 0){
				this->ExtractPixelPerUnit(dataElement);

			}else if(strcmp("mesh", dataElement->GetName()) == 0){
				this->ExtractMesh(dataElement);

			}

		 }


	}
	parser->Delete();





}

void ImportExport::ImportObjModel(std::string fname){
	if(verbose){
		cout << "in ImportObjModel: " << endl;
		cout << "\t Filename: " << fname << endl;
	}
	vtkOBJReader *objreader = vtkOBJReader::New();
	objreader->SetFileName(fname.c_str());
	objreader->Update();

	if(verbose){
	cout << "numELements: " << objreader->GetOutput()->GetNumberOfElements(1) << endl;
	}
	importmesh->DeepCopy(objreader->GetOutput());

//	 vtkSmartPointer<vtkUnsignedCharArray> colors =vtkSmartPointer<vtkUnsignedCharArray>::New();
//	 colors->SetNumberOfComponents(3);
//	 colors->SetName("Colors");
//	 // Define some colors
//	 unsigned char red[3] = {255, 0, 0};
//	 unsigned char green[3] = {0, 255, 0};
//	 unsigned char blue[3] = {0, 0, 255};
//	 unsigned char blue1[3] = {0, 0, 255};
//	 // Add the three colors we have created to the array
//	 colors->InsertNextTupleValue(red);
//	 colors->InsertNextTupleValue(green);
//	 colors->InsertNextTupleValue(blue);
//	 colors->InsertNextTupleValue(blue1);

	// importmesh->GetCellData()->GetNormals();


		normals->SetInputData(importmesh);
		normals->ComputeCellNormalsOn();
		//normals->ComputePointNormalsOff();
		// normals->AutoOrientNormalsOn();
		//normals->ConsistencyOff();
		//normals->SplittingOn();
		//normals->SetFeatureAngle(60);
		normals->Update();


	//	 SetCellcenters();

	// MediumIDs per Triangle
	//objreader->GetOutput()->GetAttributes(0)->Print(cout);

	GenMediumIDsFromObjModel(fname);

	objreader->Delete();


}
void ImportExport::ImportEmitterModel(std::string file_name){
	std::string line;
		std::ifstream infile(file_name.c_str(),std::ifstream::in);


		boost::regex emitterCh("^(\\d+\\.?\\d*) (\\d+\\.?\\d*)$");

		 if(infile){

			 while(getline( infile , line ) ) {

				 // regex:
				 boost::smatch whatEmitter;

				 if(verbose){
				 cout << "Line: " << line << endl;
				 }

				 if (boost::regex_match( line, whatEmitter, emitterCh )){

					 angles.push_back(atof(whatEmitter[1].str().c_str())); // angle
					 relRadiantIntensity.push_back(atof(whatEmitter[2].str().c_str())); // relRadiantIntensity

					 if(verbose){
						 cout << "angle | relative Radiant Intensity : " << atof(whatEmitter[1].str().c_str()) << " | ";
						 cout << atof(whatEmitter[2].str().c_str()) << endl;

					 }


				 }

			 }
		 }


		 infile.close( );



}

void ImportExport::ExportModel(std::string fname){
/*
 * <?xml version='1.0' encoding="utf-8"?>
<model>
	<name>LayerModel1</name>
	<unit>cm</unit>
	<ppu>100</ppu>
	<mediumsList>
		<medium>
			<name>layer1</name>
			<absorbtion>1</absorbtion>
			<scatter>90</scatter>
			<refraction>1</refraction>
			<anisotropy>0.75</anisotropy>
		<medium>
	</mediumsList>
	<mesh>
		<pointsList>
			<point><x>0</x><y>0</y><z>10</z></point>

		</pointsList>
		<cellsList>
			<cell>
				<points><p0>0</p0><p1>1</p1><p2>2</p2></points>
				<normal><x>0</x><y>0</y><z>-1</z></normal>
				<mediumIDs><inner>1</inner><outer>0</outer></mediumIDs>
			</cell>
		</cellsList>
	</mesh>
</model>

 */
	 ofstream writer;
	 writer.open (fname.c_str());
	 writer << "<?xml version='1.0' encoding=\"utf-8\"?>\n<model>\n";
	 // name
	 writer << "\t<name>" << ModelName << "</name>\n";
	 // unit
	 writer << "\t<unit>" << unit << "</unit>\n";
	 //ppu
	 writer << "\t<ppu>" << ppu << "</ppu>\n";

	 //mediumsList
	 writer << "\t<mediumsList>\n";

	 int numLayers = mediumIDs->GetSize();

	 for(int i = 0;i<numLayers;i++){

		 writer << "\t\t<medium>\n";
		 writer << "\t\t\t<name>" << tissueNames[i] << "</name>\n";
		 writer << "\t\t\t<absorbtion>" << absorbtionCoeffs[i] << "</absorbtion>\n";
		 writer << "\t\t\t<scatter>" << scatterCoeffs[i] << "</scatter>\n";
		 writer << "\t\t\t<refraction>" << refractionInd[i] << "</refraction>\n";
		 writer << "\t\t\t<anisotropy>" << anisotropies[i] << "</anisotropy>\n";
		 writer << "\t\t</medium>\n";

	 }
	 writer << "\t</mediumsList>\n";
	 //mesh
	 writer << "\t<mesh>\n";

	 writer << "\t\t<pointsList>\n";

	int numPoints = points->GetNumberOfPoints();

		 for(int i = 0;i<numPoints;i++){

			 double p[3];
			 points->GetPoint(i,p);
			 writer << "\t\t\t<point>";
			 writer << "<x>" << p[0] << "</x>";
			 writer << "<y>" << p[1] << "</y>";
			 writer << "<z>" << p[2] << "</z>";
			 writer << "</point>\n";

		 }
	writer << "\t\t</pointsList>\n";

	 writer << "\t\t<cellsList>\n";

			 int numVerts = Vertices->GetSize();

			 for(int i = 0;i<numVerts;i++){

				 vtkIdList *pts = vtkIdList::New();
				 Vertices->GetCell(i,pts);

				 writer << "\t\t\t<cell>";
				 writer << "\t\t\t\t<points>";
				 writer << "<p0>" << pts->GetId(0) << "</p0>";
				 writer << "<p1>" << pts->GetId(1) << "</p1>";
				 writer << "<p2>" << pts->GetId(2) << "</p2></points>\n";

				 writer << "\t\t\t\t<normal>";
				 writer << "<x>" << normal_vector[3*i] << "</x>";
				 writer << "<y>" << normal_vector[3*i+1] << "</y>";
				 writer << "<z>" << normal_vector[3*i+2] << "</z></normal>\n";

				 writer << "\t\t\t\t<mediumIDs>";
				 writer << "<inner>" << mediumIDs->GetTuple2(i)[0] << "</inner>";
				 writer << "<outer>" << mediumIDs->GetTuple2(i)[1] << "</outer></mediumIDs>\n";

				 writer << "</cell>\n";

			 }

	 writer << "\t\t</cellsList>\n";
	 writer << "\t</mesh>\n";

	 writer << "</model>";
	 writer.close();



}

void ImportExport::ExtractMediumsList(vtkXMLDataElement *mediumsList){
/*	<medium>
				<name>outer</name>
				<absorbtion>0</absorbtion>
				<scatter>0</scatter>
				<refraction>1</refraction>
				<anisotropy>0</anisotropy>
	<medium>

*/


	int numberOfComponents = mediumsList->GetNumberOfNestedElements();

	std::vector<float> extractAbsorbtion;
	std::vector<float> extractScatter;
	std::vector<float> extractRefraction;
	std::vector<float> extractAnisotropy;

	extractAbsorbtion.assign(numberOfComponents,0);
	extractScatter.assign(numberOfComponents,0);
	extractRefraction.assign(numberOfComponents,0);
	extractAnisotropy.assign(numberOfComponents,0);


	for(int i=0; i<numberOfComponents;i++){

		vtkXMLDataElement *dataElement = mediumsList->GetNestedElement(i);

		int numEl = dataElement->GetNumberOfNestedElements();

		for(int j = 0; j<numEl;j++){

		vtkXMLDataElement* mediumElement = dataElement->GetNestedElement(j);


			if(strcmp("absorbtion", mediumElement->GetName()) == 0){
				extractAbsorbtion[i] = atof(mediumElement->GetCharacterData());

			}else if(strcmp("name", mediumElement->GetName()) == 0){


			}else if(strcmp("scatter", mediumElement->GetName()) == 0){
				extractScatter[i] = atof(mediumElement->GetCharacterData());

			}else if(strcmp("refraction", mediumElement->GetName()) == 0){
				extractRefraction[i] = atof(mediumElement->GetCharacterData());

			}else if(strcmp("anisotropy", mediumElement->GetName()) == 0){
				extractAnisotropy[i] = atof(mediumElement->GetCharacterData());

			}

		}
	}


	absorbtionCoeffs.assign(extractAbsorbtion.begin(),extractAbsorbtion.end());
	scatterCoeffs.assign(extractScatter.begin(),extractScatter.end());
	refractionInd.assign(extractRefraction.begin(),extractRefraction.end());
	anisotropies.assign(extractAnisotropy.begin(),extractAnisotropy.end());


}
void ImportExport::ExtractName(vtkXMLDataElement *name){
	this->ModelName = name->GetCharacterData();
}
void ImportExport::ExtractUnit(vtkXMLDataElement *ut){
	this->unit = ut->GetCharacterData();

}
void ImportExport::ExtractMesh(vtkXMLDataElement *mesh){

/*
 * <mesh>
		<pointsList> ... </pointsList>
		<cellsList>... </cellsList>
	</mesh>
 */

	int numberMeshComp = mesh->GetNumberOfNestedElements();

	for(int i=0; i<numberMeshComp;i++){

		vtkXMLDataElement *meshComponent = mesh->GetNestedElement(i);
		if(strcmp("pointsList", meshComponent->GetName()) == 0){

			this->ExtractPointsList(meshComponent);

		}else if(strcmp("cellsList", meshComponent->GetName()) == 0){
			this->ExtractCellsList(meshComponent);

		}
	}

}
void ImportExport::ExtractPointsList(vtkXMLDataElement *pointsList){

/*
 * <pointsList>
			<point><x>0</x><y>0</y><z>10</z></point>
		</pointsList>
 */
	int numberOfPoints = pointsList->GetNumberOfNestedElements();

		for(int i=0; i<numberOfPoints;i++){

			vtkXMLDataElement *point = pointsList->GetNestedElement(i);

			double x = atof(point->GetNestedElement(0)->GetCharacterData());
			double y = atof(point->GetNestedElement(1)->GetCharacterData());
			double z = atof(point->GetNestedElement(2)->GetCharacterData());

			points->InsertNextPoint(x,y,z);

		}

}
void ImportExport::ExtractCellsList(vtkXMLDataElement *cellsList){

/*
 * <cellsList>
			<cell>
				<points><p0>0</p0><p1>1</p1><p2>2</p2></points>
				<normal><x>0</x><y>0</y><z>-1</z></normal>
				<mediumIDs><inner>1</inner><outer>0</outer></mediumIDs>
			</cell>
	</cellsList>
 */

	int numberOfCells = cellsList->GetNumberOfNestedElements();

	normal_vector.assign(numberOfCells,0);

	for(int i=0; i<numberOfCells;i++){

		vtkXMLDataElement *cell = cellsList->GetNestedElement(i);

		int numEl = cell->GetNumberOfNestedElements(); // 3: points, normal, mediumIDs

		for(int j = 0; j<numEl;j++){

			vtkXMLDataElement* cellElement = cell->GetNestedElement(j);
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

			if(strcmp("points", cellElement->GetName()) == 0){

				int indP0 = atoi(cellElement->GetNestedElement(0)->GetCharacterData());
				int indP1 = atoi(cellElement->GetNestedElement(1)->GetCharacterData());
				int indP2 = atoi(cellElement->GetNestedElement(2)->GetCharacterData());

				triangle->GetPointIds()->SetId ( 0, indP0 );
				triangle->GetPointIds()->SetId ( 1, indP1 );
				triangle->GetPointIds()->SetId ( 2, indP2 );

				Vertices->InsertNextCell(triangle);


			}else if(strcmp("normal", cellElement->GetName()) == 0){

				normal_vector[3*i] = atof(cellElement->GetNestedElement(0)->GetCharacterData());
				normal_vector[3*i+1] = atof(cellElement->GetNestedElement(1)->GetCharacterData());
				normal_vector[3*i+2] = atof(cellElement->GetNestedElement(2)->GetCharacterData());

			}else if(strcmp("mediumIDs", cellElement->GetName()) == 0){
				int indMedInner = atoi(cellElement->GetNestedElement(0)->GetCharacterData());
				int indMedOuter = atoi(cellElement->GetNestedElement(1)->GetCharacterData());

				mediumIDs->SetTuple2(i,indMedInner,indMedOuter);


			}
		}

	}

}
void ImportExport::ExtractPixelPerUnit(vtkXMLDataElement *ppu){
	this->ppu = atof(ppu->GetCharacterData());
}

// Get
void ImportExport::GetAbsorbtionCoefficients(std::vector<float> &absorbtions){
absorbtions.assign(absorbtionCoeffs.begin(),absorbtionCoeffs.end());
}
void ImportExport::GetScatterCoefficients(std::vector<float> &scatterOut){
scatterOut.assign(scatterCoeffs.begin(),scatterCoeffs.end());
}
void ImportExport::GetAnisotropies(std::vector<float> &anisotrop){
anisotrop.assign(anisotropies.begin(),anisotropies.end());
}
void ImportExport::GetRefractionIndices(std::vector<float> &refractions){
refractions.assign(refractionInd.begin(),refractionInd.end());
}
float ImportExport::GetPixelPerUnitFactor(){
return this->ppu;
}
vtkDoubleArray* ImportExport::GetMediumIDs(){
	return this->mediumIDs;
}
void ImportExport::GetMediumIDs(vtkSmartPointer<vtkDoubleArray> &medIDsOut){
	medIDsOut->DeepCopy(this->mediumIDs);
}
vtkSmartPointer<vtkPoints> ImportExport::GetPoints(){
return this->points;
}
void ImportExport::GetPoints(vtkSmartPointer<vtkPoints> &pointsOut){
	pointsOut->DeepCopy(this->points);
}
vtkSmartPointer<vtkCellArray> ImportExport::GetVertices(){
return this->Vertices;
}
void ImportExport::GetVertices(vtkSmartPointer<vtkCellArray> &vertsOut){
	vertsOut->DeepCopy(this->Vertices);
}
void ImportExport::GetNormalVectors(std::vector<float> &normalVector){

	normalVector.assign(this->normal_vector.begin(),this->normal_vector.end());
}

void ImportExport::GetSurfaceNamesVector(std::vector<std::string> &surfNamesVectorOut){
	surfNamesVectorOut.assign(surfaceIDnames.begin(),surfaceIDnames.end());

}
void ImportExport::GetSurfaceIDVector(std::vector<int> &surfIDVectorOut){
	surfIDVectorOut.assign(surfaceIDs.begin(),surfaceIDs.end());
}

void ImportExport::GetEmitterHeading(std::vector<float> &heading){
	heading.assign(emitterHeading.begin(),emitterHeading.end());
}
int ImportExport::GetFetalID(){
	return fetalID;
}
void ImportExport::GetPictureDimensions(std::vector<int> &dims){

	dims.assign(picDims.begin(),picDims.end());
}

float ImportExport::GetEmitterdistanceToSurface(){

	return emitterDistance;

}
void ImportExport::GetEmitterAngles(std::vector<float> &ang){
	ang.assign(angles.begin(),angles.end());

}
void ImportExport::GetEmitterRelativeRadiantIntensities(std::vector<float> &relRadInt){

	relRadInt.assign(relRadiantIntensity.begin(),relRadiantIntensity.end());
}
void ImportExport::GetEmitterPointInPlane(std::vector<float> &emitterPiP){
	emitterPiP.assign(emitterPoint.begin(),emitterPoint.end());
}
void ImportExport::GetSurfaceConfigNames(std::vector<std::string> &surfIDVectorOut){
	surfIDVectorOut.assign(surfaceIDnamesConfig.begin(),surfaceIDnamesConfig.end());
}


std::string ImportExport::GetUnit(){
	return this->unit;
}
std::string ImportExport::GetModelName(){
	return this->ModelName;
}

void ImportExport::GetDimensions(std::vector<float> &dimensions){

	dimensions.assign(dims.begin(),dims.end());
}
void ImportExport::GetDetectionGrid(std::vector<float> &detGrid){
	detGrid.assign(detGridDims.begin(),detGridDims.end());
}

int ImportExport::GetNumberOfPhotons(){
	return this->numPhotons;
}
int ImportExport::GetNumberOfPhotonInteractions(){
	return this->numPhotonsinteractions;
}
bool ImportExport::GetSaveFile(){
	return this->saveFile;
}
std::string ImportExport::GetSaveDir(){
	return this->saveDir;
}
double ImportExport::GetArrowScaleFactor(){
	return this->arrowScale;
}
double ImportExport::GetChance(){
	return this->chance;
}
double ImportExport::GetWeightThreshold(){
	return this->w_thresh;
}

double ImportExport::GetViewAzimuth(){
	return this->view_array[0];
}
double ImportExport::GetViewElevation(){
	return this->view_array[1];
}
double ImportExport::GetViewRoll(){
	return this->view_array[2];
}

vtkActor* ImportExport::GetActor() {
	return this->actor;
}

std::string ImportExport::GetSurfaceFilename(){
	return this->surfaceFilename;
}
std::string ImportExport::GetFetalFilename(){
	return this->fetalFilename;
}
std::string ImportExport::GetTransmissionFilename(){
	return this->transmissionFilename;
}
std::string ImportExport::GetRecordAxzFilename(){
	return this->recordAxzFilename;
}
std::string ImportExport::GetRecordAyzFilename(){
	return this->recordAyzFilename;
}

std::string ImportExport::GetScaleString(){
	return this->scaleString;
}
std::string ImportExport::GetCoordSystemString(){
	return this->coordsysString;
}
std::string ImportExport::GetColormapString(){
	return this->colormapString;
}

// Set
void ImportExport::SetVerbosity(bool verb){
	verbose = verb;
}
void ImportExport::SetVerbosityTrue(){
	verbose = true;
}
void ImportExport::SetVerbosityFalse(){
	verbose = false;
}
void ImportExport::SetAbsorbtionCoefficients(std::vector<float> &absorbtionsIn){
	absorbtionCoeffs.assign(absorbtionsIn.begin(),absorbtionsIn.end());
}
void ImportExport::SetScatterCoefficients(std::vector<float> &scatterIn){
	scatterCoeffs.assign(scatterIn.begin(),scatterIn.end());
}
void ImportExport::SetAnisotropies(std::vector<float> &anisotropiesIn){
	anisotropies.assign(anisotropiesIn.begin(),anisotropiesIn.end());
}
void ImportExport::SetRefractionIndices(std::vector<float> &refractionsIn){
	refractionInd.assign(refractionsIn.begin(),refractionsIn.end());
}
void ImportExport::SetPixelPerUnitFactor(float ppuIn){
	this->ppu = ppuIn;
}
void ImportExport::SetMediumIDs(vtkDoubleArray* mediumIDsIn){
	mediumIDs->DeepCopy(mediumIDsIn);
}
void ImportExport::SetPoints(vtkSmartPointer<vtkPoints> pointsIn){
	points->DeepCopy(pointsIn);
}
void ImportExport::SetVertices(vtkSmartPointer<vtkCellArray> vertsIn){
	Vertices->DeepCopy(vertsIn);
}
void ImportExport::SetNormalVectors(std::vector<float> &normalsVectorIn){
	normal_vector.assign(normalsVectorIn.begin(),normalsVectorIn.end());
}
void ImportExport::SetModelName(char* modelname){
	this->ModelName = modelname;
}
void ImportExport::SetUnit(char* unitIn){
	this->unit = unitIn;
}

void ImportExport::SetCellcenters(){


	//GetCellcenters()->VertexCellsOn();
	GetCellcenters()->SetInputData(this->GetNormalsOutput());
	GetCellcenters()->Update();

}


vtkCellCenters* ImportExport::GetCellcenters(){

	return this->Cellcenter;

}
vtkPolyData* ImportExport::GetNormalsOutput(){
	return this->normals->GetOutput();
}

void ImportExport::GenMediumIDsFromObjModel(std::string file_name){
	/* MediumIDs per Tri: {(inner,outer),(inner,outer),...}
	 * outerIDs, innerIDs, tissueNames
	 */

	std::string line;
	std::ifstream infile(file_name.c_str(),std::ifstream::in);
//	std::vector<std::string> objectGroups;
//	std::vector<int> objectGroupFaceNum;
	//int numFaces = 0;
	int outerTissueID = 0;
	int currentID = 0;
	int innerTissueID = 0;
	bool isSurface = false;
	std::string surfName;


	mediumIDs->SetNumberOfComponents(2);


	boost::regex objGr("^o (\\w+)");
	boost::regex objFaces("^f (.+)");
	boost::regex vertGroups("^g (\\d+)$");
	boost::regex vertGroupsDigitWord("^g (\\d+)\\_(\\w+)");


	 if(infile){

		 while(getline( infile , line ) ) {

			 // regex:
			 boost::smatch whatGr;
			 boost::smatch whatFaces;
			 boost::smatch whatVertGroups;
			 boost::smatch whatVertGroupsWord;
			 if(verbose){
			 cout << "Line: " << line << endl;
			 }
			 if (boost::regex_match( line, whatGr, objGr )){
				 // tissue

				 int tns = tissueNames.size();
				 for(int i=0;i<tns;i++){
					 if(verbose){
					 cout << "TissueName(" << i << "): " << tissueNames[i] << " | ";
					 cout << "found: " << whatGr[1].str().c_str() << endl;
					 }
					 if(strcmp(whatGr[1].str().c_str(),tissueNames[i].c_str())==0){
						 innerTissueID = i;// i = innerTissueID
						 break;
					 }
				 }


			 }else if (boost::regex_match( line, whatFaces, objFaces )){
				 // Triangles

				 mediumIDs->InsertNextTuple2(innerTissueID,outerTissueID);

				 if(isSurface){
					 currentID = mediumIDs->GetNumberOfTuples()-1;
					 surfaceIDs.push_back(currentID);
					 surfaceIDnames.push_back(surfName);
					 if(verbose){
					 cout << "clear0: " << surfName.c_str() << endl;
					 surfName = "";
					 cout << "clear1: " << surfName.c_str() << endl;
					 }else{
						 surfName = "";
					 }
				 }




			}else if (boost::regex_match( line, whatVertGroups, vertGroups )){
				// outer Tissue ID
				isSurface = false;
				outerTissueID = atoi(whatVertGroups[1].str().c_str()); // ==outerTissueID


			}else if (boost::regex_match( line, whatVertGroupsWord, vertGroupsDigitWord )){

				isSurface = true;
				if(verbose){
				cout << "isSurface: " << whatVertGroupsWord[2].str() << " | outerTissue: " << whatVertGroupsWord[1].str() << endl;
				}
				outerTissueID = atoi(whatVertGroupsWord[1].str().c_str());

				surfName = whatVertGroupsWord[2].str();


			}

		 }
	 }


	 infile.close( );


}


void ImportExport::PrintNormalsCenters(){

cout << "Number of Components: " << GetCellcenters()->GetOutput()->GetPointData()->GetNumberOfComponents() << endl;
cout << "Number of Normals: " << GetCellcenters()->GetOutput()->GetPoints()->GetNumberOfPoints() << endl;
cout << "Number of MediumIDs-Tuples: " << this->mediumIDs->GetNumberOfTuples() << endl;
cout << "Number of MediumIDs-Components: " << this->mediumIDs->GetNumberOfComponents() << endl;

	for(int i = 0;i<GetCellcenters()->GetOutput()->GetPoints()->GetNumberOfPoints();i++){

		double p[3];
		GetCellcenters()->GetOutput()->GetPoints()->GetPoint(i,p);
		double n[3];
		GetNormals()->GetTuple(i,n);


		for(int i=0;i<this->mediumIDs->GetNumberOfTuples();i++){

			double* m = this->mediumIDs->GetTuple2(i);
			double p0[3];
			double p1[3];
			double p2[3];
			this->importmesh->GetCell(i)->GetPoints()->GetPoint(0,p0);
			this->importmesh->GetCell(i)->GetPoints()->GetPoint(1,p1);
			this->importmesh->GetCell(i)->GetPoints()->GetPoint(2,p2);
			cout << " Medium (inside, outside): " << i << "(" << m[0] << ", " << m[1] << ") ";
			cout << "\tP0(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ")\t| ";
			cout << "\tP1(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ")\t| ";
			cout << "\tP2(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ")" << endl;


		}

	}
}


void ImportExport::PrintCells(){
	double numTup = normals->GetOutput()->GetCellData()->GetNormals()->GetNumberOfTuples();

		cout << "NumNormals: " << numTup << endl;
		int numCells = importmesh->GetNumberOfCells();

		for(int i=0;i<numCells;i++){

			double p0[3];
			double p1[3];
			double p2[3];

			importmesh->GetCell(i)->GetPoints()->GetPoint(0,p0);
			importmesh->GetCell(i)->GetPoints()->GetPoint(1,p1);
			importmesh->GetCell(i)->GetPoints()->GetPoint(2,p2);

			cout << "Cell(" << i << "): P0(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ") | ";
			cout << "P1(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ") | ";
			cout << "P2(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ") | " << endl;

		}


}


void ImportExport::PrintNormals(){
	int numCells = importmesh->GetNumberOfCells();
	for(int i = 0;i<numCells;i++){

			double n0[3];
			normals->GetOutput()->GetCellData()->GetNormals()->GetTuple(i,n0);
			cout << "Normal(" << i << "): (" << n0[0] << ", " << n0[1] << ", " << n0[2] << ")" << endl;


		}

}

void ImportExport::PrintImport(){

	cout << "########### Import ##################" << endl;
	PrintMediumIDs();
	PrintNormals();
	PrintCells();

}

void ImportExport::PrintMediumIDs(){

		int numTissues = mediumIDs->GetNumberOfTuples();
		int numcomps = mediumIDs->GetNumberOfComponents();
		cout << "MediumIDs:\n################################### \n";
		cout << "Number of Tuples: " << numTissues << "\n";
		cout << "Number of Components: " << numcomps << "\n";

		for(int i = 0;i<numTissues;i++){
			cout << "Cell(" << i << "): ( ";
			double val[numcomps];
			mediumIDs->GetTupleValue(i,val);
			for(int j = 0;j<numcomps;j++){
				cout << val[j] << " ";

			}
			cout << ")\n";
		}

		cout << "\n################################### \n" << endl;


}
