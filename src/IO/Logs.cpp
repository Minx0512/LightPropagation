/*
 * Logs.cpp
 *
 *  Created on: Jul 26, 2015
 *      Author: matthias
 */

#include "Logs.hpp"
#include "../ImageFilter/ImageFilter.hpp"

#include <vector>
#include <string>

#include <vtkImageActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageProperty.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageMapper3D.h>
#include <vtkImageCast.h>


#include <vtkImageLogarithmicScale.h>
#include <vtkImageShiftScale.h>
#include <vtkSmartPointer.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageData.h>
#include <vtkImageExport.h>
#include <vtkPNGWriter.h>

Logs::Logs() {
	// TODO Auto-generated constructor stub



	imageSFT = vtkSmartPointer<vtkImageData>::New();
	imageRFTd = vtkSmartPointer<vtkImageData>::New();

	pixelPerUnitFactor = 100;

}

Logs::~Logs() {
	// TODO Auto-generated destructor stub
}





void Logs::SetSaveDir(std::string directory){
	saveDir = directory;
}

void Logs::SetDimensions(std::vector<float> &dims){
	dimensions.assign(dims.begin(),dims.end());
}
void Logs::SetGridDimensions(std::vector<float> &gdims){
	gridDims.assign(gdims.begin(),gdims.end());
}
void Logs::SetPixelPerUnitFactor(float factor) {
	this->pixelPerUnitFactor = factor;
}
void Logs::SetNumberOfPhotons(int number){
	this->numberPhotons = number;
}
void Logs::SetRTdxyLogs(vtkSmartPointer<vtkImageData> image){


	imageRFTd = image;


}

void Logs::SetSFTLogs(vtkSmartPointer<vtkImageData> image){

	int ext[6];
	image->GetExtent(ext);

	imageSFT = image;



	for(int x=ext[0];x<=ext[1];x++){
		for(int y=ext[2];y<=ext[3];y++){

			surfRxyLog.insert(surfRxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,0));
			fetalRxyLog.insert(fetalRxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,1));
			transTxyLog.insert(transTxyLog.end(),image->GetScalarComponentAsFloat(x,y,0,2));
		}
	}


}


void Logs::SetRxyLogs(std::vector<float> &surffetaltransRxy){

	int s = surffetaltransRxy.size();


	for(int i = 0;i<s;i+=3){

		surfRxyLog.insert(surfRxyLog.end(),surffetaltransRxy[i]);
		fetalRxyLog.insert(fetalRxyLog.end(),surffetaltransRxy[i+1]);
		transTxyLog.insert(transTxyLog.end(),surffetaltransRxy[i+2]);


	}





}

void Logs::SetAxzLog(std::vector<float> &Axz){
	AxzLog.assign(Axz.begin(),Axz.end());
}

void Logs::SetAyzLog(std::vector<float> &Ayz){
	AyzLog.assign(Ayz.begin(),Ayz.end());
}


void Logs::SetAxzCoords(std::vector<int> &coords){
	AxzCoords.assign(coords.begin(),coords.end());
}

void Logs::SetAyzCoords(std::vector<int> &coords){
	AyzCoords.assign(coords.begin(),coords.end());
}
void Logs::SetScaleString(std::string scale){
	scaleString=scale;




}
void Logs::SetCoordSystemString(std::string coord){
	coordsysString = coord;
}
void Logs::SetColormapString(std::string colormap){
	colormapString = colormap;
}


void Logs::WritePositionLog(std::string saveDir,std::string fname, std::vector<float> &vec){
	SetSaveDir(saveDir);
	WritePositionLog(fname, vec);

}

void Logs::WriteAbsorbtionLog(std::string saveDirIn,std::string fname, std::vector<float> &vec){

	SetSaveDir(saveDir);
	WriteAbsorbtionLog(fname,vec);

}

void Logs::WriteConstInFile(std::string saveDir,std::string fname, std::string constName, float constant){

	SetSaveDir(saveDir);
	WriteConstInFile(fname,constName,constant);

}

/*
void Logs::SaveSurfFetTransArr(std::string fnameSurf,std::string fnameFetal,std::string fnameTrans){
	std::ofstream SurfFile;
	std::ofstream FetalFile;
	std::ofstream TransFile;



	std::string filenameSurf = saveDir;
	std::string filenameFetal = saveDir;
	std::string filenameTrans = saveDir;
	filenameSurf.append(fnameSurf.begin(),fnameSurf.end());
	filenameFetal.append(fnameFetal.begin(),fnameFetal.end());
	filenameTrans.append(fnameTrans.begin(),fnameTrans.end());


	std::cout << "WriteArrInFile Filename: " << fnameSurf << ", " << fnameFetal << "," << fnameTrans << " | Number ofElements: " << vec.size() <<  std::endl;

	SurfFile.open (filenameSurf.c_str());
	FetalFile.open (filenameFetal.c_str());
	TransFile.open (filenameTrans.c_str());

	int s = vec.size();
	int phiBins = int(360/gridDims[4]);


/*
 *  phi -->
 *  r |
 *    v
 */

//std::cout << "s: " << s<< std::endl;
//std::cout << "phiBins: " << phiBins<< std::endl;
/*
	for(int j=0;j<phiBins;j++){
		SurfFile << j <<"\t";
		FetalFile << j <<"\t";
		TransFile << j <<"\t";

	}

	SurfFile << std::endl;
	FetalFile << std::endl;
	TransFile << std::endl;
*/
/*
	for(int j=0;j<phiBins;j++){
		SurfFile << vec[0]<<"\t";
		FetalFile << vec[1]<<"\t";
		TransFile << vec[2]<<"\t";

	}
	SurfFile << std::endl;
	FetalFile << std::endl;
	TransFile << std::endl;

int count = 0;
int rs = 1;
	for(int i=1;i<s;i+=3){
		//std::cout<< " i_1: "<< i;
		SurfFile << vec[i] << "\t";
		FetalFile << vec[i+1] << "\t";
		TransFile << vec[i+2] << "\t";
		count++;


		if(count==phiBins){
			SurfFile << std::endl;
			FetalFile << std::endl;
			TransFile << std::endl;

			//std::cout<< " count: "<< count << " i: " << i << std::endl;
			count = 0;
			rs++;
		}


	}
	std::cout<< "r_s: "<< rs << std::endl;

	SurfFile.close();
	FetalFile.close();
	TransFile.close();



}
*/

void Logs::SaveSurfaceArr(std::string fname, int LOG){

	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;
	filenameSurf.append(fname.begin(),fname.end());

	float max = 0;


	SurfFile.open (filenameSurf.c_str());

	std::vector<float> arr;

	switch(LOG){
		case Rxy_Log_SURF:
			arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
		case Rxy_Log_TRANS:
			arr.assign(transTxyLog.begin(),transTxyLog.end());
			break;
		case Rxy_Log_FETAL:
			arr.assign(fetalRxyLog.begin(),fetalRxyLog.end());
			break;
		default: arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
	}

	std::cout << "\t WriteArrInFile Filename: " << fname << " | Number ofElements: " << arr.size() <<  std::endl;


		int s = arr.size();
		int phiBins = int(360/gridDims[4]);

	/*
	 *  phi -->
	 *  r |
	 *    v
	 */

		for(int j=0;j<phiBins;j++){
			SurfFile << arr[0]<<"\t";

			if(max<arr[0]){max = arr[0];}

		}
		SurfFile << std::endl;


	int count = 0;
	int rs = 1;
		for(int i=1;i<s;i++){
			//std::cout<< " i_1: "<< i;
			SurfFile << arr[i] << "\t";
			if(max<arr[i]){max = arr[i];}
			count++;


			if(count==phiBins){
				SurfFile << std::endl;


				//std::cout<< " count: "<< count << " i: " << i << std::endl;
				count = 0;
				rs++;
			}


		}

		SurfFile.close();





}
void Logs::SaveSurfaceArr(std::string fname,vtkSmartPointer<vtkImageData> image){

	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;
	filenameSurf.append(fname.begin(),fname.end());

	SurfFile.open (filenameSurf.c_str());

	int ext[6];
	image->GetExtent(ext);
	int dims[3];
	image->GetDimensions(dims);

	std::cout << "\t WriteArrInFile Filename: " << fname << " | Number ofElements: " << dims[0]*dims[1] <<  std::endl;

	for(int n=ext[2];n<=ext[3];n++){

		for(int m=ext[0];m<=ext[1];m++){

			SurfFile <<	image->GetScalarComponentAsFloat(m,n,0,0) <<"\t";
		}

		SurfFile << endl;


	}
	SurfFile.close();

}

void Logs::SaveRecordRd(std::string fname,vtkSmartPointer<vtkImageData> image){

	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;

		filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_Rd.txt"));

		SurfFile.open (filenameSurf.c_str());

		int ext[6];
		image->GetExtent(ext);
		int dims[3];
		image->GetDimensions(dims);

		std::cout << "\t SaveRecordRd Filename: " << fname << " | Number ofElements: " << dims[0]*dims[1] <<  std::endl;

		for(int n=ext[2];n<=ext[3];n++){

			for(int m=ext[0];m<=ext[1];m++){

				SurfFile <<	image->GetScalarComponentAsFloat(m,n,0,0)*Scale2D(m,n) <<"\t";
			}

			SurfFile << endl;


		}
		SurfFile.close();





}

void Logs::SaveAbsorbtionLog(std::string fname,std::vector<float> &absobArr, std::vector<int> &coords){

	std::ofstream absorbFile;
	std::string filename = saveDir;

	filename.append(fname.begin(),fname.end());

	std::cout << "WriteArrInFile Filename: " << fname << "| Number ofElements: " << absobArr.size() <<  std::endl;

	absorbFile.open (filename.c_str());

	int s = absobArr.size();

	/*
	 *  int idAxz = dimx/2 + c1*dimx+c0;
	 *  int idAyz = dimy/2 + c1*dimy+c0;
	 */

	for(int i=0;i<s;i++){
		if(absobArr[i]!=0){
			absorbFile << coords[2*i] << "\t" << coords[2*i+1] << "\t" << absobArr[i] << std::endl;
		}

	}


	absorbFile.close();


}

void Logs::CalcReflectancePerRadius(std::string fname, int LOG){
	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;

	filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_perRadius.txt"));

	SurfFile.open (filenameSurf.c_str());

	std::vector<float> arr;

		switch(LOG){
		case Rxy_Log_SURF:
			arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
		case Rxy_Log_TRANS:
			arr.assign(transTxyLog.begin(),transTxyLog.end());
			break;
		case Rxy_Log_FETAL:
			arr.assign(fetalRxyLog.begin(),fetalRxyLog.end());
			break;
		default: arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
		}


	std::cout << "\t Calculae Reflectance per radius: " << fname.substr(0,fname.find_last_of(".")).append("_perRadius.txt") << " | Number ofElements: " << arr.size() <<  std::endl;

	int s = arr.size();
	int phiBins = int(360/gridDims[4]);

	int radiusIdx = 0;
	SurfFile << radiusIdx << "\t" << arr[0]<< endl;
	int count = 0;
	float reflectance = 0;

	radiusIdx++;
	for(int i=1;i<s;i++){

		reflectance+=arr[i];
		count++;

		if(count==phiBins){
			SurfFile << radiusIdx << "\t" << reflectance << endl;
			radiusIdx++;
			reflectance = 0;
			count = 0;

		}

	}

	SurfFile.close();



}

float Logs::ScalePerRadius(int ir){


	float scale = 1;
	// in [cm]
	float dr = gridDims[3];

	scale = 2*M_PI*dr*dr*numberPhotons;
	scale = 1/((ir+0.5)*scale);

	return scale;

}

float Logs::ScalePerRadius1D(int ir){

	double dr = gridDims[3];
	float scale1 = 2.0*M_PI*dr*dr*numberPhotons;
		/* area is 2*PI*[(ir+0.5)*dr]*dr.*/
		/* ir+0.5 to be added. */
	float scale2 = 1.0/((ir+0.5)*scale1);

	return scale2;
}



void Logs::CalcReflectancePerRadius(std::string fname,vtkSmartPointer<vtkImageData> image){
	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;

	filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_perRadius.txt"));

	SurfFile.open (filenameSurf.c_str());
	int ext[6];
	image->GetExtent(ext);
	int dims[3];
	image->GetDimensions(dims);

	cout << "\t Calculate Reflectance per radius: " << fname.substr(0,fname.find_last_of(".")).append("_perRadius.txt");
	cout <<  " | Number ofElements: " << dims[0]*dims[1] <<  std::endl;

//cout << "Extent: m: "<< ext[0] << ", " << ext[1] << " | n: " << ext[2] << ", " << ext[3] << endl;
	float sum = 0;
	float sum2 = 0;
	for(int r=ext[0];r<=ext[1];r++){
			// radiusIdx
		float reflectance = 0;

		for(int a=ext[2];a<=ext[3];a++){
			//alphaIdx
			reflectance+=image->GetScalarComponentAsFloat(r,a,0,0);

		}

		SurfFile << r*gridDims[3] << "\t" << reflectance*ScalePerRadius1D(abs(r)) << endl;
	//	sum+=reflectance;
	//	sum2+=reflectance*ScalePerRadius(r);
	//	cout << "1/ScalePerRadius(" << r << "): " << 1/ScalePerRadius(r) << endl;

	}

//	cout << "Sum: " << sum << " | sum2: " << sum2 << endl;
	SurfFile.close();



}
void Logs::CalcReflectancePerAngle(std::string fname, int LOG){
	std::ofstream file;

	std::string filenameSurf = saveDir;
	filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_perAngle.txt"));
	file.open (filenameSurf.c_str());

	std::vector<float> arr;

	switch(LOG){
		case Rxy_Log_SURF:
			arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
		case Rxy_Log_TRANS:
			arr.assign(transTxyLog.begin(),transTxyLog.end());
			break;
		case Rxy_Log_FETAL:
			arr.assign(fetalRxyLog.begin(),fetalRxyLog.end());
			break;
		default: arr.assign(surfRxyLog.begin(),surfRxyLog.end());
			break;
	}


	std::cout << "\t Calculae Reflectance per angle: " << fname.substr(0,fname.find_last_of(".")).append("_perAngle.txt") << " | Number ofElements: " << arr.size() <<  std::endl;

	int s = arr.size();
	int phiBins = int(360/gridDims[4]);

//	int radiusIdx = 0;

	//int count = 0;

	float reflectance = 0;

	int R = int(floor((s-1)/phiBins));

	for(int phi=1;phi<=phiBins;phi++){

		for(int r = 1;r<=R;r++){

			int i = (r-1)*phiBins+phi;
			reflectance+=arr[i];

		}

		file << phi << "\t" << reflectance << endl;

		reflectance = 0;

	}

	file.close();


}




void Logs::ScaleRdTt(vtkSmartPointer<vtkImageData> imageRdRA, std::vector<float> &RdR, std::vector<float> &RdA) {

	int ext[6];
	imageRdRA->GetExtent(ext);
	int dims[3];
	imageRdRA->GetDimensions(dims);


  for(int ir=ext[0]; ir<=ext[1]; ir++){
    for(int ia=ext[2]; ia<=ext[3]; ia++) {

      float reflectance = imageRdRA->GetScalarComponentAsFloat(ir,ia,0,0);
      reflectance*=Scale2D(ir, ia);
      imageRdRA->SetScalarComponentFromFloat(ir,ia,0,0,reflectance);

    }
  }


  for(int ir=0; ir<RdR.size(); ir++) {
    RdR[ir]*=ScalePerRadius1D(ir);
  }

  for(int ia=0; ia<RdA.size(); ia++) {
    RdA[ia]*=ScalePerAngle1D(ia);
  }

}

float Logs::ScalePerAngle1D(int ia){

	double da = gridDims[5]*M_PI/180; // in [rad]
	float scale1  = 2.0*M_PI*da*numberPhotons;
	  /* solid angle is 2*PI*sin(a)*da. sin(a) to be added. */
	float scale2 = 1.0/(sin((ia+0.5)*da)*scale1);

	return scale2;
}

float Logs::Scale2D(int ir, int ia){

	 double dr = gridDims[3];
	 double da = gridDims[5]*M_PI/180; // in [rad]

	float scale1 = 4.0*M_PI*M_PI*dr*sin(da/2)*dr*numberPhotons;
		/* The factor (ir+0.5)*sin(2a) to be added. */
	float scale2 = 1.0/((ir+0.5)*sin(2.0*(ia+0.5)*da)*scale1);

return scale2;

}

void Logs::CalcReflectancePerAngle(std::string fname,vtkSmartPointer<vtkImageData> image){
	std::ofstream SurfFile;

	std::string filenameSurf = saveDir;

	filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_perAngle.txt"));

	SurfFile.open (filenameSurf.c_str());
	int ext[6];
	image->GetExtent(ext);
	int dims[3];
	image->GetDimensions(dims);

	cout << "\t Calculae Reflectance per angle: " << fname.substr(0,fname.find_last_of(".")).append("_perAngle.txt");
	cout <<  " | Number ofElements: " << dims[0]*dims[1];
	cout << " Extent: [" << ext[0] << ", " << ext[1] << "],["<< ext[2] << ", " << ext[3]<< "]"<< endl;

	for(int n=ext[2];n<=ext[3];n++){
				//alphaIdx
		float reflectance = 0;
		for(int m=ext[0];m<=ext[1];m++){		// radiusIdx


			reflectance+=image->GetScalarComponentAsFloat(m,n,0,0);


		}

		SurfFile << n*gridDims[5] << "\t" << reflectance*ScalePerAngle1D(n) << endl;

	}

	SurfFile.close();



}




void Logs::CreateSFT_XYPicture(std::string fnameS,vtkSmartPointer<vtkImageData> imageData){


	ImageFilter *imFilt = new ImageFilter;
	imFilt->SetColorBarWidth(15);
	imFilt->SetBlackWidth(10);
	imFilt->SetInput(imageData);
	imFilt->SetColorToBar(1,1,1);
	imFilt->SetColormapString(colormapString);
	imFilt->SetCoordinatesystemToCartesian();
	imFilt->SetScaleString(scaleString);

	imFilt->Update();

	//imFilt->PrintLogScale();


	std::string filename = saveDir;

	std::string postStrg = "";

	EvalScaleForFname(imFilt->GetScaleMode());
	EvalColormapForFname(imFilt->GetColormapMode());

	postStrg.append(scaleStringFname).append(colormapStringFname);
	postStrg.append(".png");


	filename.append(fnameS.substr(0,fnameS.find_last_of(".")).append(postStrg));

	std::cout << "\t Write PNG file: " << filename.substr(filename.find_last_of("/")+1,filename.length()) << std::endl;

	/*
		vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
		writer->SetFileName(filename.c_str());
		writer->SetInputData(imFilt->GetOutput());
		writer->Write();
	*/

		WritePNG(filename,imFilt->GetOutput());


}

void Logs::CreateAbsorbtionPicture(std::string fname,vtkSmartPointer<vtkImageData> imageData){
	int ext[6];

	vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();

	imageData->GetExtent(ext);

		image->SetExtent(ext);
		image->SetOrigin(imageData->GetOrigin());
		image->SetSpacing(imageData->GetSpacing());
		image->AllocateScalars(VTK_FLOAT,1);

		float scale = 2*M_PI*gridDims[3]*gridDims[3]*gridDims[2]*numberPhotons;

	for(int r=ext[0];r<=ext[1];r++){
		for(int z=ext[2];z<=ext[3];z++){
			float v = imageData->GetScalarComponentAsFloat(r,z,0,0);
			v /= scale*(abs(r)+0.5);
			image->SetScalarComponentFromFloat(r,z,0,0,v);

		}
	}

	ImageFilter *imFilt = new ImageFilter;
	imFilt->SetColorBarWidth(15);
	imFilt->SetBlackWidth(10);
	imFilt->SetInput(image);
	imFilt->SetColorToBar(1,1,1);
	imFilt->SetColormapString(colormapString);
	imFilt->SetCoordinatesystemToCartesian();
	imFilt->SetScaleString(scaleString);

	imFilt->Update();

	//imFilt->PrintLogScale();


	std::string filename = saveDir;
	std::string postStrg = "_Z";
	postStrg.append(std::to_string(gridDims[2]));

	EvalScaleForFname(imFilt->GetScaleMode());
	EvalCoordsysForFname(imFilt->GetCoordinatesystem());
	EvalColormapForFname(imFilt->GetColormapMode());


	postStrg.append(coordsysStringFname).append(scaleStringFname).append(colormapStringFname);
	postStrg.append(".png");

	filename.append(fname.substr(0,fname.find_last_of(".")).append(postStrg));

	std::cout << "\t Write PNG file: " << filename.substr(filename.find_last_of("/")+1,filename.length()) << std::endl;



		WritePNG(filename,imFilt->GetOutput());


}


void Logs::CreateSFT_XYPicture(std::string fnameS,int LOG, int dimR, int dimPhi){

	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
	imageData->SetExtent(0,dimR-1,0,dimPhi-1,0,0);
	imageData->SetSpacing(1,1,1);
	imageData->SetOrigin(0,0,0);
	imageData->AllocateScalars(VTK_FLOAT,1);
/*
	vtkSmartPointer<vtkImageData> imageDataF = vtkSmartPointer<vtkImageData>::New();
	imageDataF->SetExtent(0,dimR,0,dimPhi,0,1);
	imageDataF->SetSpacing(1,1,1);
	imageDataF->SetOrigin(0,0,0);
	imageDataF->AllocateScalars(VTK_FLOAT,1);


	vtkSmartPointer<vtkImageData> imageDataT = vtkSmartPointer<vtkImageData>::New();
	imageDataT->SetExtent(0,dimR,0,dimPhi,0,1);
	imageDataT->SetSpacing(1,1,1);
	imageDataT->SetOrigin(0,0,0);
	imageDataT->AllocateScalars(VTK_FLOAT,1);
*/
//	cout<< "Begin" <<endl;

	int ext[6];
	imageData->GetExtent(ext);

	for(int n=ext[2];n<=ext[3];n++){
		for(int m=ext[0];m<=ext[1];m++){

			imageData->SetScalarComponentFromFloat(m,n,0,0,0);

		}

	}

	std::vector<float> arr;

	switch(LOG){
	case Rxy_Log_SURF:
		arr.assign(surfRxyLog.begin(),surfRxyLog.end());
		break;
	case Rxy_Log_TRANS:
		arr.assign(transTxyLog.begin(),transTxyLog.end());
		break;
	case Rxy_Log_FETAL:
		arr.assign(fetalRxyLog.begin(),fetalRxyLog.end());
		break;
	default: arr.assign(surfRxyLog.begin(),surfRxyLog.end());
		break;
	}


	int s = arr.size();
	int R = 1;
	int phi = 0;
	float max = 0;

	for(int j=0;j<dimPhi;j++){

		imageData->SetScalarComponentFromFloat(0,j,0,0,arr[0]);
	//	cout << surfRxyLog[0] << " ";
		if(max<arr[0]){max = arr[0];}

	}

//	cout << endl;
	for(int i=1;i<s;i++){


		imageData->SetScalarComponentFromFloat(R,phi,0,0,arr[i]);

		if(max<arr[i]){max = arr[i];}

		phi++;
		if(phi==dimPhi){
			R++;
			phi = 0;
	//		cout << endl;
		}

		if(R>dimR-1){
			break;
		}

	}

	ImageFilter *imFilt = new ImageFilter;
	imFilt->SetColorBarWidth(15);
	imFilt->SetBlackWidth(10);
	imFilt->SetInput(imageData);
	imFilt->SetColorToBar(1,1,1);
	imFilt->SetColormapString(colormapString);
	imFilt->SetCoordsystemString(coordsysString);
	imFilt->SetScaleString(scaleString);

	imFilt->Update();

	//imFilt->PrintLogScale();


	std::string filename = saveDir;

	std::string postStrg = "_R";
	postStrg.append(std::to_string(dimR));


	EvalScaleForFname(imFilt->GetScaleMode());
	EvalCoordsysForFname(imFilt->GetCoordinatesystem());
	EvalColormapForFname(imFilt->GetColormapMode());




	postStrg.append(coordsysStringFname).append(scaleStringFname).append(colormapStringFname);
	postStrg.append(".png");


	filename.append(fnameS.substr(0,fnameS.find_last_of(".")).append(postStrg));

	std::cout << "\t Write PNG file: " << filename.substr(filename.find_last_of("/")+1,filename.length()) << std::endl;

/*
	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(filename.c_str());
	writer->SetInputData(imFilt->GetOutput());
	writer->Write();
*/

	WritePNG(filename,imFilt->GetOutput());


/*

	 // Create an image actor
	  vtkSmartPointer<vtkImageActor> imageActor = vtkSmartPointer<vtkImageActor>::New();
	  imageActor->GetMapper()->SetInputData(imFilt->GetOutput());
	  imageActor->GetProperty()->SetInterpolationTypeToNearest();

	  // Visualize
	  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	  renderer->AddActor(imageActor);
	  renderer->ResetCamera();

	  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	  renderWindow->AddRenderer(renderer);

	  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
	    vtkSmartPointer<vtkRenderWindowInteractor>::New();
	  vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();

	  renderWindowInteractor->SetInteractorStyle(style);

	  renderWindowInteractor->SetRenderWindow(renderWindow);
	  renderWindowInteractor->Initialize();
	  renderWindowInteractor->Start();
/*

cout<<"FilenameS: " << filenameS << endl;
/*
	vtkSmartPointer<vtkPNGWriter> writerFet = vtkSmartPointer<vtkPNGWriter>::New();
	writerFet->SetFileName(filenameFet.c_str());
	writerFet->SetInputData(imageDataF);
	writerFet->Write();

	vtkSmartPointer<vtkPNGWriter> writerT = vtkSmartPointer<vtkPNGWriter>::New();
	writerT->SetFileName(filenameT.c_str());
	writerT->SetInputData(imageDataT);
	writerT->Write();
*/

}



void Logs::CreateAbsorbtionPicture(std::string fname,std::vector<float> &absobArr, std::vector<int> &coords, int dim0, int dim1){



		vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

		imageData->SetExtent(-dim0/2,dim0/2,-dim1+1,0,0,0);
		imageData->SetSpacing(1, 1, 1);
		imageData->SetOrigin(0, 0, 0);
		imageData->AllocateScalars(VTK_FLOAT,1);

		int s = absobArr.size();

	//	cout << "CreateAbsorbtionPicture | s = " << s << "| coords.size: " << coords.size() << endl;

	//	cout << "ImageData before: " << endl;

		int ext[6];
		imageData->GetExtent(ext);

		for(int n=ext[2];n<=ext[3];n++){
			for(int m=ext[0];m<=ext[1];m++){

				imageData->SetScalarComponentFromFloat(m,n,0,0,0);

			}

		}


		for(int i=0;i<s;i++){

			int x = int(coords[2*i]);
			int y = int(-coords[2*i+1]);

			float val = (absobArr[i]>0) ? absobArr[i] : 0;
			if(val>0){
				imageData->SetScalarComponentFromFloat(x,y,0,0,val);
			}

		}

		ImageFilter *imFilt = new ImageFilter;
		imFilt->SetColorBarWidth(15);
		imFilt->SetBlackWidth(10);
		imFilt->SetInput(imageData);
		imFilt->SetColorToBar(1,1,1);

		imFilt->SetScaleString(scaleString);
		imFilt->SetCoordinatesystemToCartesian();
		imFilt->SetColormapString(colormapString);


		imFilt->Update();


		//fname = [].txt
		std::string filename = saveDir;
		std::string postStrg = "_Z";
		postStrg.append(std::to_string(gridDims[2]));


		EvalScaleForFname(imFilt->GetScaleMode());
		EvalCoordsysForFname(imFilt->GetCoordinatesystem());
		EvalColormapForFname(imFilt->GetColormapMode());




		postStrg.append(coordsysStringFname).append(scaleStringFname).append(colormapStringFname);
		postStrg.append(".png");

		filename.append(fname.substr(0,fname.find_last_of(".")).append(postStrg));

	std::cout << "\t Write PNG file: " << filename.substr(filename.find_last_of("/")+1,filename.length()) << std::endl;


		WritePNG(filename,imFilt->GetOutput());

}

void Logs::FluenceZ(std::string fname, vtkSmartPointer<vtkImageData> fluenceImage){

	// fluenceImage(r,z)

	int ext[6];
	fluenceImage->GetExtent(ext);

	std::ofstream file;

	std::string filenameSurf = saveDir;

	filenameSurf.append(fname.substr(0,fname.find_last_of(".")).append("_perZ.txt"));
	file.open (filenameSurf.c_str());

	cout << "\t Calculate Fluence per Z: " << fname.substr(0,fname.find_last_of(".")).append("_perZ.txt");
	cout << " Extent: [" << ext[0] << ", " << ext[1] << "],["<< ext[2] << ", " << ext[3]<< "]"<< endl;

	for(int z=ext[2];z<=ext[3];z++){
		float fluenceZ = 0;
		for(int r=ext[0];r<=ext[1];r++){

		//	float dAlpha = 2*M_PI*(abs(r)+0.5)*gridDims[3]*gridDims[3];
		//	dAlpha *= fluenceImage->GetScalarComponentAsFloat(r,z,0,0);
			fluenceZ+=fluenceImage->GetScalarComponentAsFloat(r,z,0,0);


		}

		file << abs(z)*gridDims[2] << "\t" << fluenceZ/(gridDims[2]*numberPhotons) << endl;

	}


	file.close();

}

void Logs::EvalScaleForFname(int scale){

	switch(scale){
	case IF_LOG10:scaleStringFname = "_Log10";
		break;
	case IF_LIN:scaleStringFname = "_Linear";
		break;
	default:scaleStringFname = "_Log10";
		break;

	}



}
void Logs::EvalCoordsysForFname(int system){

	switch(system){
	case IF_COORDS_CARTESIAN:coordsysStringFname = "_Cartesian";
		break;

	case IF_COORDS_RADIAL:coordsysStringFname = "_Radial";
		break;

	default:coordsysStringFname = "_Cartesian";
		break;

	}

}
void Logs::EvalColormapForFname(int map){


	switch(map){
	case COLORMAP_C2W:colormapStringFname = "_C2W";
		break;
	case COLORMAP_RAINBOW: colormapStringFname = "_RAINBOW";
		break;
	case COLORMAP_C2WL: colormapStringFname = "_C2WL";
		break;
	case COLORMAP_GRAYSCALE: colormapStringFname = "_GRAYSCALE";
		break;
	case COLORMAP_BLACKBODY: colormapStringFname = "_BLACKBODY";
		break;
	default:colormapStringFname = "_C2W";
		break;


	}

}



void Logs::CreateAbsorbtionPicture(std::string fname,int Aplane){

	int dim0,dim1;
//	cout << "CreateAbsorbtionPicture 1" << endl;
	switch(Aplane){
	case Axz_Log:

		dim0 = dimensions[0]/gridDims[0];
		dim1 = dimensions[2]/gridDims[2];

		//cout << "CreateAbsorbtionPicture : dimx, dimz "<<dim0<<", "<<dim1 << endl;

			CreateAbsorbtionPicture(fname,AxzLog,AxzCoords,dim0,dim1);

		break;

	case Ayz_Log:

		dim0 = dimensions[1]/gridDims[1];
		dim1 = dimensions[2]/gridDims[2];
	//		cout << "CreateAbsorbtionPicture : dimy, dimz "<<dim0<<", "<<dim1 << endl;
			CreateAbsorbtionPicture(fname,AyzLog,AyzCoords,dim0,dim1);

		break;

	default:
		break;


	}



}

void Logs::WritePNG(std::string fname,vtkSmartPointer<vtkImageData> imageData){

	vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
	writer->SetFileName(fname.c_str());
	writer->SetInputData(imageData);
	writer->Write();

}

