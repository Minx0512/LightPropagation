/*
 * ReadXML.cpp
 *
 *  Created on: Mar 7, 2015
 *      Author: matthias
 * Model
 * 	--> Size (width x height x depth)
 * 	--> Layer --> densities, layerdepths
 *
 */

#include "ReadXML.hpp"


#include <string.h>
#include <cstring>
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkStringArray.h"
#include "vtkXMLDataParser.h"

ReadXML::ReadXML() {
	verbose = false;
	configParser = vtkXMLDataParser::New();
	ModelXMLParser = vtkXMLDataParser::New();
	// Model properties
	dimensions_array = vtkDoubleArray::New();
	rotation_array = vtkDoubleArray::New();
	view_array = vtkDoubleArray::New();
	reductionfactors_array = vtkDoubleArray::New();

	// Tissue properties
	tissueNames_array = vtkStringArray::New();
	layerdepths_array = vtkIntArray::New(); // Layerdepth in px
	thickness_array = vtkDoubleArray::New();; // Layerthickness in cm
	densities_array = vtkDoubleArray::New(); // densities in grey scale [0,255]
	anisotropies_array = vtkDoubleArray::New();;
	absorbtion_array = vtkDoubleArray::New();;
	scatter_array = vtkDoubleArray::New();;
	refractions_array = vtkDoubleArray::New();;

	modelID_int = 0;
	modeID_int = 0;
	ppunit = 0;
	filename = "";
	modelFileName = "";
	saveDir = "";
	emitterFileName = "";
	numTissues = 0; // number of tissues
	currTissueID = 0; // current tissue id

	numPhotons = 0;
	numPhotonsinteractions = 0;
	chance = 0.1;
	w_thresh = 0.001;

	arrScaleFactor = 1;
	saveFile = 0;

	scaleString = "";
	coordsysString = "";
	colormapString = "";



}


ReadXML::~ReadXML() {
	// TODO Auto-generated destructor stub



}

void ReadXML::ParseModelXMLFile(std::string fname){

		 ModelXMLParser->SetFileName(fname.c_str());

		 if(ModelXMLParser->Parse()==1){

			 vtkXMLDataElement *root = ModelXMLParser->GetRootElement();

			 this->ExtractVars(root);


		 }


	// parser1->Delete();
}

void ReadXML::ParseConfig(std::string fname){


	 configParser->SetFileName(fname.c_str());

	 if(configParser->Parse()==1){

		vtkXMLDataElement *root = configParser->GetRootElement();
		this->ExtractConfigVars(root);

	 }


}


void ReadXML::ExtractConfigVars(vtkXMLDataElement *rt){

	 int numberOfComponents = rt->GetNumberOfNestedElements();

	 for(int i=0; i<numberOfComponents;i++){

		 vtkXMLDataElement *dataElement = rt->GetNestedElement(i);

		if(strcmp("modelfile", dataElement->GetName()) == 0){

			this->ExtractModelFile(dataElement);
			if(verbose){cout << "Extract model filename: " << this->modelFileName << endl;}

		}else if(strcmp("arrowscalefactor", dataElement->GetName()) == 0){

			this->ExtractArrowScaleFactor(dataElement);

		}else if(strcmp("SimulationParams", dataElement->GetName()) == 0){

			this->ExtractSimulationParams(dataElement);


		}else if(strcmp("saveLogFiles", dataElement->GetName()) == 0){

			this->ExtractSaveLogFiles(dataElement);


		}else if(strcmp("Logs", dataElement->GetName()) == 0){

			this->ExtractLogs(dataElement);


		}else if(strcmp("saveDir", dataElement->GetName()) == 0){

			this->ExtractSaveDir(dataElement);


		}else if(strcmp("view", dataElement->GetName()) == 0){

			this->ExtractInitViewRotation(dataElement);


		}else if(strcmp("model", dataElement->GetName()) == 0){

			this->ExtractModel(dataElement);


		}else if(strcmp("Emittermodel", dataElement->GetName()) == 0){
			this->ExtractEmittermodel(dataElement);


		}else if(strcmp("analyse", dataElement->GetName()) == 0){

			this->ExtractAnalyse(dataElement);

		}
	 }

}


void ReadXML::ExtractVars(vtkXMLDataElement *rt){

	 int numberOfComponents = rt->GetNumberOfNestedElements();
	 numTissues = 0;

	 for(int i = 0;i<numberOfComponents;i++){
		 vtkXMLDataElement *dataElement = rt->GetNestedElement(i);

		 if(strcmp("tissue", dataElement->GetName()) == 0){
			 numTissues++;
		 }

	 }

	 this->InitTissueArrays(numTissues);

	 for(int i=0; i<numberOfComponents;i++){

		 vtkXMLDataElement *dataElement = rt->GetNestedElement(i);

		if(strcmp("model", dataElement->GetName()) == 0){

			this->ExtractModelID(dataElement);

		}else if(strcmp("mode", dataElement->GetName()) == 0){

			this->ExtractModeID(dataElement);

		}else if(strcmp("PixelPerUnit", dataElement->GetName()) == 0){

			this->ExtractPixelPerUnitFactor(dataElement);

		}else if(strcmp("dimensions", dataElement->GetName()) == 0){

			this->ExtractDimensions(dataElement);

		}else if(strcmp("rotation", dataElement->GetName()) == 0){

			this->ExtractRotations(dataElement);
		}else if(strcmp("reductionfactors", dataElement->GetName()) == 0){

			this->ExtractReductionfactors(dataElement);
		}else if(strcmp("tissue", dataElement->GetName()) == 0){

			int numEl = dataElement->GetNumberOfNestedElements();

			for(int j = 0; j<numEl;j++){

				vtkXMLDataElement* tissueElement = dataElement->GetNestedElement(j);


				if(strcmp("ID", tissueElement->GetName()) == 0){
					this->ExtractTissueID(tissueElement);

				}else if(strcmp("name", tissueElement->GetName()) == 0){
					this->ExtractTissueName(tissueElement);

				}else if(strcmp("layerdepth", tissueElement->GetName()) == 0){
					this->ExtractLayerdepths(tissueElement);

				}else if(strcmp("thickness", tissueElement->GetName()) == 0){
					this->ExtractThickness(tissueElement);

				}else if(strcmp("density", tissueElement->GetName()) == 0){
					this->ExtractDensities(tissueElement);

				}else if(strcmp("anisotropy", tissueElement->GetName()) == 0){
					this->ExtractAnisotropies(tissueElement);

				}else if(strcmp("mua", tissueElement->GetName()) == 0){
					this->ExtractAbsorbtionCoefficients(tissueElement);

				}else if(strcmp("mus", tissueElement->GetName()) == 0){
					this->ExtractScatterCoefficients(tissueElement);

				}else if(strcmp("refraction", tissueElement->GetName()) == 0){
					this->ExtractRefractionIndices(tissueElement);

				}
			}
		}

	 }

}


void ReadXML::SetFilename(std::string fname){

	this->filename = fname;

}

std::string ReadXML::GetEmitterFile(){
	return this->emitterFileName;
}

// Get Model properties
std::string ReadXML::GetModelFilename(){
	return this->modelFileName;
}
int ReadXML::GetModeID(){
	return this->modeID_int;
}
int ReadXML::GetModelID(){
	return this->modelID_int;
}

float ReadXML::GetPixelPerUnitFactor(){
	return this->ppunit;
}

vtkDoubleArray* ReadXML::GetRotations(){
	return this->rotation_array;
}

void ReadXML::GetRotations(std::vector<float> &rot){
	int s = this->rotation_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->rotation_array->GetValue(i));
	}
	rot.assign(values.begin(),values.end());
}

void ReadXML::GetViewRotations(std::vector<float> &rot){
	int s = this->view_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->view_array->GetValue(i));
	}
	rot.assign(values.begin(),values.end());
}

double ReadXML::GetViewAzimuth(){
	return this->view_array->GetValue(0);
}
double ReadXML::GetViewElevation(){
	return this->view_array->GetValue(1);
}
double ReadXML::GetViewRoll(){
	return this->view_array->GetValue(2);
}

vtkDoubleArray* ReadXML::GetDimensions(){
	return this->dimensions_array;
}
void ReadXML::GetDimensions(std::vector<float> &dims){

	int s = this->dimensions_array->GetSize();
	std::vector<double> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = this->dimensions_array->GetValue(i);
	}
	dims.assign(values.begin(),values.end());
}


void ReadXML::GetPictureDimensions(std::vector<int> &dims){

	dims.assign(picDims.begin(),picDims.end());
}

int ReadXML::GetFetalID(){
	return fetalID;
}

void ReadXML::GetDetectionGrid(std::vector<float> &detGrid){
	detGrid.assign(detGridDims.begin(),detGridDims.end());
}

vtkDoubleArray* ReadXML::GetReductionFactors(){
	return this->reductionfactors_array;
}


double ReadXML::GetReductionFactors(int idx){
	return this->reductionfactors_array->GetValue(idx);
}

// Get Tissue properties

void ReadXML::GetTissueNames(std::vector<std::string> &lNames){

	int s = this->tissueNames_array->GetSize();

	std::string c;
	std::vector<std::string> values;

	for(int i= 0;i<s;i++){
		c = this->tissueNames_array->GetValue(i);
		values.push_back(c);

	}

	lNames.assign(values.begin(),values.end());

}
void ReadXML::GetLayerdepth(std::vector<int> &ldpth){

	int s = this->layerdepths_array->GetSize();
	std::vector<double> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = this->layerdepths_array->GetValue(i);
	}
	ldpth.assign(values.begin(),values.end());
}
void ReadXML::GetThicknesses(std::vector<float> &thicks){
	int s = this->thickness_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->thickness_array->GetValue(i));
	}
	thicks.assign(values.begin(),values.end());
}
void ReadXML::GetDensities(std::vector<float> &denses){
	int s = this->densities_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->densities_array->GetValue(i));
	}
	denses.assign(values.begin(),values.end());
}
double ReadXML::GetDensities(int dens){
	return this->densities_array->GetValue(dens);
}
void ReadXML::GetAnisotropies(std::vector<float> &anisotropies){
	int s = this->anisotropies_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->anisotropies_array->GetValue(i));
	}
	anisotropies.assign(values.begin(),values.end());
}
void ReadXML::GetAbsorbtionCoefficients(std::vector<float> &absorbtions){
	int s = this->absorbtion_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->absorbtion_array->GetValue(i));
	}
	absorbtions.assign(values.begin(),values.end());
}
void ReadXML::GetScatterCoefficients(std::vector<float> &scatter){
	int s = this->scatter_array->GetSize();
	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->scatter_array->GetValue(i));
	}
	scatter.assign(values.begin(),values.end());
}
void ReadXML::GetRefractionIndices(std::vector<float> &refractions) {
	int s = this->refractions_array->GetSize();

	std::vector<float> values;
	values.assign(s,0);
	for(int i = 0;i<s;i++){
		values[i] = float(this->refractions_array->GetValue(i));
	}
	refractions.assign(values.begin(),values.end());
}

int ReadXML::GetNumberOfPhotons(){
	return this->numPhotons;
}

int ReadXML::GetNumberOfPhotonInteractions(){
	return this->numPhotonsinteractions;
}

bool ReadXML::GetSaveFile(){
	return bool(saveFile);
}
std::string ReadXML::GetSaveDir(){
	return this->saveDir;
}

std::string ReadXML::GetSurfaceFilename(){
	return this->surfaceFilename;
}
std::string ReadXML::GetFetalFilename(){
	return this->fetalFilename;
}
std::string ReadXML::GetTransmissionFilename(){
	return this->transmissionFilename;
}
std::string ReadXML::GetRecordAxzFilename(){
	return this->recordAxzFilename;
}
std::string ReadXML::GetRecordAyzFilename(){
	return this->recordAyzFilename;
}

std::string ReadXML::GetScaleString(){
	return this->scaleString;
}
std::string ReadXML::GetCoordSystemString(){
	return this->coordsysString;
}
std::string ReadXML::GetColormapString(){
	return this->colormapString;
}

double ReadXML::GetArrowScaleFactor(){
	return this->arrScaleFactor;
}
double ReadXML::GetChance(){
	return this->chance;
}
double ReadXML::GetWeightThreshold(){
	return this->w_thresh;
}

void ReadXML::GetSurfaceNamesVector(std::vector<std::string> &surfNamesVectorOut){
	surfNamesVectorOut.assign(surfaceIDnames.begin(),surfaceIDnames.end());

}
void ReadXML::GetSurfaceIDVector(std::vector<int> &surfIDVectorOut){
	surfIDVectorOut.assign(surfaceIDs.begin(),surfaceIDs.end());
}
float ReadXML::GetEmitterdistanceToSurface(){

	return emitterDistance;

}
void ReadXML::GetEmitterPointInPlane(std::vector<float> &emitterPiP){
	emitterPiP.assign(emitterPoint.begin(),emitterPoint.end());
}
void ReadXML::GetSurfaceConfigNames(std::vector<std::string> &surfNamesOut) {

	surfNamesOut.assign(surfaceIDnamesConfig.begin(),surfaceIDnamesConfig.end());
}

void ReadXML::GetEmitterHeading(std::vector<float> &heading){
	heading.assign(emitterHeading.begin(),emitterHeading.end());
}

// Extract vars
void ReadXML::ExtractReductionfactors(vtkXMLDataElement *redf){

	int num = redf->GetNumberOfNestedElements();

		this->reductionfactors_array->Resize(num);

		for(int i = 0;i<num;i++){

			double value = atof(redf->GetNestedElement(i)->GetCharacterData());
			this->reductionfactors_array->SetValue(i,value);
		}



}
void ReadXML::ExtractDimensions(vtkXMLDataElement *dims){
	int num = dims->GetNumberOfNestedElements();

	this->dimensions_array->Resize(num);

	for(int i = 0;i<num;i++){

		float value = atof(dims->GetNestedElement(i)->GetCharacterData());
		this->dimensions_array->SetValue(i,value);
	}

}
void ReadXML::ExtractRotations(vtkXMLDataElement *rots){
	int num = rots->GetNumberOfNestedElements();

	double value = 0;
	this->rotation_array->Resize(3);
	this->rotation_array->SetValue(0,0);
	this->rotation_array->SetValue(1,0);
	this->rotation_array->SetValue(2,0);


	for(int i = 0;i<num;i++){

		if(strcmp("yaw",rots->GetNestedElement(i)->GetName())==0) {

			value = atof(rots->GetNestedElement(i)->GetCharacterData());
			this->rotation_array->SetValue(0,value);

		} else if(strcmp("pitch",rots->GetNestedElement(i)->GetName())==0) {

				value = atof(rots->GetNestedElement(i)->GetCharacterData());
				this->rotation_array->SetValue(1,value);

		} else if(strcmp("roll",rots->GetNestedElement(i)->GetName())==0) {

				value = atof(rots->GetNestedElement(i)->GetCharacterData());
				this->rotation_array->SetValue(2,value);

		}
	}
}
void ReadXML::ExtractModelID(vtkXMLDataElement *model){
	this->modelID_int = atoi(model->GetCharacterData());
}
void ReadXML::ExtractModeID(vtkXMLDataElement *mode){
	this->modeID_int = atoi(mode->GetCharacterData());
}
void ReadXML::ExtractPixelPerUnitFactor(vtkXMLDataElement *mode){
	this->ppunit = atof(mode->GetCharacterData());
}
void ReadXML::ExtractModelFile(vtkXMLDataElement *mfile){
	this->modelFileName = mfile->GetCharacterData();

}
void ReadXML::ExtractSaveLogFiles(vtkXMLDataElement *save){
	this->saveFile = atoi(save->GetCharacterData());
}
void ReadXML::ExtractSaveDir(vtkXMLDataElement *save){
	this->saveDir = save->GetCharacterData();
}
void ReadXML::ExtractNumberOfPhotonInteractions(vtkXMLDataElement *numInter){
	this->numPhotonsinteractions = atoi(numInter->GetCharacterData());
}
void ReadXML::ExtractNumberOfPhotons(vtkXMLDataElement *numPhotonsDE){
	this->numPhotons = atoi(numPhotonsDE->GetCharacterData());
}
void ReadXML::ExtractSimulationParams(vtkXMLDataElement *rot){

	int num = rot->GetNumberOfNestedElements();

		for(int i = 0;i<num;i++){
			vtkXMLDataElement *dataElement = rot->GetNestedElement(i);
			if(strcmp("numberOfPhotons", dataElement->GetName()) == 0){

				this->ExtractNumberOfPhotons(dataElement);

			}else if(strcmp("numPhotonsInteractions", dataElement->GetName()) == 0){

				this->ExtractNumberOfPhotonInteractions(dataElement);

			}else if(strcmp("chance", dataElement->GetName()) == 0){
				chance = atof(dataElement->GetCharacterData());


			}else if(strcmp("weightThreshold", dataElement->GetName()) == 0){
				w_thresh = atof(dataElement->GetCharacterData());

			}
		}

}

void ReadXML::ExtractInitViewRotation(vtkXMLDataElement *rot){

	int num = rot->GetNumberOfNestedElements();

		this->view_array->Resize(num);

		for(int i = 0;i<num;i++){

			vtkXMLDataElement *dataElement = rot->GetNestedElement(i);
			double value = atof(rot->GetNestedElement(i)->GetCharacterData());
			if(strcmp("Azimuth", dataElement->GetName()) == 0){

				this->view_array->SetValue(0,value);

			}else if(strcmp("Elevation", dataElement->GetName()) == 0){

				this->view_array->SetValue(1,value);

			}else if(strcmp("Roll", dataElement->GetName()) == 0){

				this->view_array->SetValue(2,value);

			}

		}



}

void ReadXML::ExtractArrowScaleFactor(vtkXMLDataElement *scale){
	arrScaleFactor = atof(scale->GetCharacterData());
}
void ReadXML::ExtractTopBottom(vtkXMLDataElement *topbottom){


}


void ReadXML::ExtractModel(vtkXMLDataElement* model){

	int num = model->GetNumberOfNestedElements();

	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = model->GetNestedElement(i);

		if(strcmp("modelfile", dataElement->GetName()) == 0){
			this->ExtractModelFile(dataElement);

		}else if(strcmp("unit", dataElement->GetName()) == 0){


		}else if(strcmp("PixelPerUnit", dataElement->GetName()) == 0){
			this->ExtractPixelPerUnitFactor(dataElement);

		}else if(strcmp("dimensions", dataElement->GetName()) == 0){
			this->ExtractDimensions(dataElement);

		}else if(strcmp("tissueList", dataElement->GetName()) == 0){
			this->ExtractTissueList(dataElement);

		}
	}


}


void ReadXML::ExtractAnalyse(vtkXMLDataElement* analyse){

	int num = analyse->GetNumberOfNestedElements();

	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = analyse->GetNestedElement(i);

		if(strcmp("boundaries", dataElement->GetName()) == 0){

			ExtractAnalyseBoundaries(dataElement);
		}
	}


}


void ReadXML::ExtractAnalyseBoundaries(vtkXMLDataElement* boundaries){

	int num = boundaries->GetNumberOfNestedElements();

	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = boundaries->GetNestedElement(i);

		if(strcmp("surfName", dataElement->GetName()) == 0){

			ExtractAnalyseSurfaceNames(dataElement);
		}
	}

}

void ReadXML::ExtractAnalyseSurfaceNames(vtkXMLDataElement* names){

	surfaceIDnamesConfig.push_back(names->GetCharacterData());


}
void ReadXML::ExtractEmittermodel(vtkXMLDataElement *emittermodel){
	int num = emittermodel->GetNumberOfNestedElements();

		for(int i = 0;i<num;i++){
			vtkXMLDataElement *dataElement = emittermodel->GetNestedElement(i);

			if(strcmp("file", dataElement->GetName()) == 0){

				ExtractEmitterfile(dataElement);
			}else if(strcmp("position", dataElement->GetName()) == 0){

				emitterPoint.push_back(atof(dataElement->GetNestedElement(0)->GetCharacterData()));
				emitterPoint.push_back(atof(dataElement->GetNestedElement(1)->GetCharacterData()));
				emitterPoint.push_back(atof(dataElement->GetNestedElement(2)->GetCharacterData()));

			}else if(strcmp("distance", dataElement->GetName()) == 0){

				emitterDistance = atof(dataElement->GetCharacterData());

			}else if(strcmp("heading", dataElement->GetName()) == 0){
				emitterHeading.assign(3,0);

				emitterHeading[0] = atof(dataElement->GetNestedElement(0)->GetCharacterData());
				emitterHeading[1] = atof(dataElement->GetNestedElement(1)->GetCharacterData());
				emitterHeading[2] = atof(dataElement->GetNestedElement(2)->GetCharacterData());


			}
		}


}
void ReadXML::ExtractEmitterfile(vtkXMLDataElement *emitterfile){

	this->emitterFileName = emitterfile->GetCharacterData();

}

void ReadXML::ExtractLogs(vtkXMLDataElement *logsroot){
/*
	<scale>Log10</scale>
	<coordsystem>radial</coordsystem>
	<colormap>c2w</colormap>
	*/

	int num = logsroot->GetNumberOfNestedElements();
	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = logsroot->GetNestedElement(i);

		if(strcmp("saveLogFiles", dataElement->GetName()) == 0){

			this->ExtractSaveLogFiles(dataElement);


		}else if(strcmp("fetalID", dataElement->GetName()) == 0){
			fetalID = atoi(dataElement->GetCharacterData());

		}else if(strcmp("DetectionGrid", dataElement->GetName()) == 0){
			this->ExtractLogsDetectionGrid(dataElement);
		}else if(strcmp("filenames", dataElement->GetName()) == 0){

			this->ExtractFilenames(dataElement);

		}else if(strcmp("fileDims", dataElement->GetName()) == 0){

			this->ExtractLogsFileDimensions(dataElement);


		}else if(strcmp("scale", dataElement->GetName()) == 0){

			scaleString = dataElement->GetCharacterData();


		}else if(strcmp("coordsystem", dataElement->GetName()) == 0){

			coordsysString = dataElement->GetCharacterData();


		}else if(strcmp("colormap", dataElement->GetName()) == 0){

			colormapString = dataElement->GetCharacterData();


		}


	}



}

void ReadXML::ExtractLogsFileDimensions(vtkXMLDataElement *fileDimensions){
	picDims.assign(4,1);

	int num = fileDimensions->GetNumberOfNestedElements();
		for(int i = 0;i<num;i++){
			vtkXMLDataElement *dataElement = fileDimensions->GetNestedElement(i);

			if(strcmp("R", dataElement->GetName()) == 0){

				picDims[3] = atoi(dataElement->GetCharacterData());

			}else if(strcmp("x", dataElement->GetName()) == 0){
				picDims[0] = atoi(dataElement->GetCharacterData());
			}else if(strcmp("y", dataElement->GetName()) == 0){
				picDims[1] = atoi(dataElement->GetCharacterData());
			}else if(strcmp("z", dataElement->GetName()) == 0){
				picDims[2] = atoi(dataElement->GetCharacterData());
			}


		}



}

void ReadXML::ExtractLogsDetectionGrid(vtkXMLDataElement *logsDetGrid){
	detGridDims.assign(6,1);

	int num = logsDetGrid->GetNumberOfNestedElements();
	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = logsDetGrid->GetNestedElement(i);

		if(strcmp("x", dataElement->GetName()) == 0){
			detGridDims[0] = atof(dataElement->GetCharacterData());
		}else if(strcmp("y", dataElement->GetName()) == 0){
			detGridDims[1] = atof(dataElement->GetCharacterData());
		}else if(strcmp("z", dataElement->GetName()) == 0){
			detGridDims[2] = atof(dataElement->GetCharacterData());
		}else if(strcmp("r", dataElement->GetName()) == 0){
			detGridDims[3] = atof(dataElement->GetCharacterData());
		}else if(strcmp("rho", dataElement->GetName()) == 0){
			detGridDims[5] = atof(dataElement->GetCharacterData());
		}else if(strcmp("phi", dataElement->GetName()) == 0){
			detGridDims[4] = atof(dataElement->GetCharacterData());
		}
	}


}
void ReadXML::ExtractFilenames(vtkXMLDataElement *filenamesroot){
	int num = filenamesroot->GetNumberOfNestedElements();
	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = filenamesroot->GetNestedElement(i);

		if(strcmp("surface", dataElement->GetName()) == 0){
			surfaceFilename = dataElement->GetCharacterData();
		}else if(strcmp("fetal", dataElement->GetName()) == 0){
			fetalFilename = dataElement->GetCharacterData();
		}else if(strcmp("transmission", dataElement->GetName()) == 0){
			transmissionFilename = dataElement->GetCharacterData();
		}else if(strcmp("absorbtionXZ", dataElement->GetName()) == 0){
			recordAxzFilename = dataElement->GetCharacterData();
		}else if(strcmp("absorbtionYZ", dataElement->GetName()) == 0){
			recordAyzFilename = dataElement->GetCharacterData();
		}

	}

}

// Tissue properties
void ReadXML::InitTissueArrays(int numT){


	this->refractions_array->Resize(numT);
	this->reductionfactors_array->Resize(numT);
	this->absorbtion_array->Resize(numT);
	this->anisotropies_array->Resize(numT);
	this->scatter_array->Resize(numT);
	this->tissueNames_array->Resize(numT);
	this->densities_array->Resize(numT);
	this->layerdepths_array->Resize(numT);
	this->thickness_array->Resize(numT);


}
void ReadXML::ExtractDensities(vtkXMLDataElement *dense){

	double value = atof(dense->GetCharacterData());
	this->densities_array->SetValue(currTissueID,value);

}
void ReadXML::ExtractLayerdepths(vtkXMLDataElement *ldeps){

		int value = atoi(ldeps->GetCharacterData());
		this->layerdepths_array->SetValue(currTissueID,value);
}
void ReadXML::ExtractRefractionIndices(vtkXMLDataElement* refractions){

	double value = atof(refractions->GetCharacterData());
	this->refractions_array->SetValue(currTissueID,value);
}
void ReadXML::ExtractAbsorbtionCoefficients(vtkXMLDataElement*absorbtions){
	double value = atof(absorbtions->GetCharacterData());
	this->absorbtion_array->SetValue(currTissueID,value);

}
void ReadXML::ExtractScatterCoefficients(vtkXMLDataElement* scatter){
	double value = atof(scatter->GetCharacterData());
	this->scatter_array->SetValue(currTissueID,value);
}
void ReadXML::ExtractAnisotropies(vtkXMLDataElement* anisotropies){
	double value = atof(anisotropies->GetCharacterData());
	this->anisotropies_array->SetValue(currTissueID,value);
}
void ReadXML::ExtractTissueID(vtkXMLDataElement* tnr){
	currTissueID = atoi(tnr->GetCharacterData());
}
void ReadXML::ExtractTissueName(vtkXMLDataElement* name){
	vtkStdString value = name->GetCharacterData();
	this->tissueNames_array->SetValue(currTissueID,value);
}
void ReadXML::ExtractThickness(vtkXMLDataElement* thickness){
	double value = atof(thickness->GetCharacterData());
	this->thickness_array->SetValue(currTissueID,value);
}

void ReadXML::ExtractTissueList(vtkXMLDataElement* tissueList){

	int num = tissueList->GetNumberOfNestedElements();
	InitTissueArrays(num);
	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = tissueList->GetNestedElement(i);
		if(strcmp("tissue", dataElement->GetName()) == 0){
			this->ExtractTissue(dataElement);

		}
	}


}

void ReadXML::ExtractTissue(vtkXMLDataElement* tissue){

	int num = tissue->GetNumberOfNestedElements();


	for(int i = 0;i<num;i++){
		vtkXMLDataElement *dataElement = tissue->GetNestedElement(i);

		if(strcmp("ID", dataElement->GetName()) == 0){

			this->ExtractTissueID(dataElement);
		}else if(strcmp("name", dataElement->GetName()) == 0){

				this->ExtractTissueName(dataElement);

		}else if(strcmp("layerdepth", dataElement->GetName()) == 0){

				this->ExtractLayerdepths(dataElement);

		}else if(strcmp("thickness", dataElement->GetName()) == 0){

				this->ExtractThickness(dataElement);

		}else if(strcmp("density", dataElement->GetName()) == 0){

				this->ExtractDensities(dataElement);

		}else if(strcmp("anisotropy", dataElement->GetName()) == 0){

				this->ExtractAnisotropies(dataElement);

		}else if(strcmp("mua", dataElement->GetName()) == 0){

				this->ExtractAbsorbtionCoefficients(dataElement);

		}else if(strcmp("mus", dataElement->GetName()) == 0){

				this->ExtractScatterCoefficients(dataElement);

		}else if(strcmp("refraction", dataElement->GetName()) == 0){

				this->ExtractRefractionIndices(dataElement);

		}
	}
}


void ReadXML::PrintArray(std::string header, vtkDoubleArray* arr){

	int num = arr->GetSize();

		cout << header << ": (" ;

		for(int i = 0;i<num;i++){

			cout << arr->GetValue(i) << " ";

		}

		cout << ")" << endl;
}

void ReadXML::PrintArray(std::string header, vtkIntArray* arr){

	int num = arr->GetSize();

		cout << header << ": (" ;

		for(int i = 0;i<num;i++){

			cout << arr->GetValue(i) << " ";

		}

		cout << ")" << endl;
}

void ReadXML::PrintArray(std::string header, vtkStringArray* arr){

	int num = arr->GetSize();

		cout << header << ": (" ;

		for(int i = 0;i<num;i++){

			cout << arr->GetValue(i) << " ";

		}

		cout << ")" << endl;
}



void ReadXML::PrintVars(){


	// ModelID
		cout << "ModelID: "<< this->modelID_int << endl;
	// ModeID
		cout << "ModeID: "<< this->modeID_int << endl;

	// Dimensions
	PrintArray("Dimensions",this->dimensions_array);

	// Densities
	PrintArray("Densities",this->densities_array);
	// Layerdepths
	PrintArray("Layerdepths",this->layerdepths_array);
	// Rotations
	PrintArray("Rotations (yaw,pitch,roll)",this->rotation_array);
	// Reductionfacors
	PrintArray("Reductionfactors",this->reductionfactors_array);

	PrintArray("Layer names",this->tissueNames_array);
	PrintArray("Layer thickness",this->thickness_array);
	PrintArray("Refraction indices",this->refractions_array);
	PrintArray("Absorbtion coeffs",this->absorbtion_array);


	PrintArray("Scatter coeffs",this->scatter_array);
	PrintArray("Anisotropy coeefs",this->anisotropies_array);





}




void ReadXML::SetModelID(int model_ID){
	this->modelID_int = model_ID;
}
void ReadXML::SetRotations(double* rot){
	int size = sizeof(rot);
	this->densities_array->SetArray(rot,size,0);
}
void ReadXML::SetLayerdepth(int* ldeps){
	int size = sizeof(ldeps);
	this->layerdepths_array->SetArray(ldeps,size,0);
}
void ReadXML::SetDensities(double* dens){
	int size = sizeof(dens);
	this->densities_array->SetArray(dens,size,0);
}
void ReadXML::SetDimensions(double* dims){

	this->dimensions_array->SetArray(dims,3,0);

}
void ReadXML::SetModeID(int mode_ID){
	this->modeID_int = mode_ID;
}

void ReadXML::SetVerbosityTrue(){
	verbose = true;
}
void ReadXML::SetVerbosityFalse(){
	verbose = false;
}

