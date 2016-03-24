#include <iostream>
#include <fstream>
#include <stdlib.h>

/* #include "main.hpp"
//#include "../Model/createModel.hpp"
//#include "../Model/VoxToMesh.hpp"
*/

#include "../Simulation/CalculatePropagation.hpp"
#include "../IO/ImportExport.hpp"
#include "../IO/Logs.hpp"
#include "../Analysis/Analysis.hpp"
#include "../Model/Emitter.hpp"
#include "../Interface/CLI.hpp"

#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkRenderer.h"
#include "vtkInteractorStyleTrackballCamera.h"

/*
#include "vtkVolume.h"

#include "vtkVolumeProperty.h"
#include "vtkSmartPointer.h"
#include "vtkSmartVolumeMapper.h"
#include "vtkPolyDataMapper.h"

#include "vtkFixedPointVolumeRayCastMapper.h"

#include "vtkContourFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataMapper.h"
#include "vtkStripper.h"
#include "vtkImageToPolyDataFilter.h"

#include "vtkVolumeProperty.h"
#include "vtkPiecewiseFunction.h"

#include "vtkDecimatePro.h"
#include "vtkQuadricDecimation.h"
#include "vtkDelaunay3D.h"
#include "vtkVolume16Reader.h"

#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include "vtkCellCenters.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkIntArray.h"

//#include "vtkProperty.h"
//#include "vtkMapper.h"
//#include "vtkStripper.h"
//#include "vtkActor.h"
*/

using namespace std;


int main(int argc, char *argv[]){


	CLI *interface = new CLI;
	if(interface->ParseArgs(argc,argv)==EXIT_SUCCESS){


	 /*if (argc < 8)
	    {

	    return EXIT_FAILURE;
	    }
*/
//	 char* filename = argv[1];

	// bool verbose = atoi(argv[3]);
	// bool graphics = atoi(argv[2]);

	// bool printImport = atoi(argv[4]);
	// bool printEmitter = atoi(argv[5]);
	// bool printTimeSpent = atoi(argv[6]);
	// bool performSimulation = atoi(argv[7]);
	// bool analysis = atoi(argv[8]);
	// int deviceType = atoi(argv[9]);



	 // Import: Mesh from other Sources like COMSOL/ Blender --> wavefront object file
	 // Export: PolyData, XML, COMSOL, etc.


	std::vector<float> refractionInd;
	std::vector<float> absorbtionCoeffs;
	std::vector<float> anisotropies;
	std::vector<float> scatterCoeffs;


	std::vector<float> RandNmbersLog;
	std::vector<float> SurfFetalTransRxy;
	std::vector<float> RecordAxz;
	std::vector<float> RecordAyz;
	std::vector<int> RAxzCoords;
	std::vector<int> RAyzCoords;
	std::vector<float> dims;
	std::vector<float> detGrid;


	ImportExport *importE = new ImportExport;
	importE->SetVerbosity(interface->GetVerbose());
	importE->SetConfig(interface->GetConfigFilename());
	importE->ParseFile();

	importE->GetRefractionIndices(refractionInd);
	importE->GetAbsorbtionCoefficients(absorbtionCoeffs);
	importE->GetAnisotropies(anisotropies);
	importE->GetScatterCoefficients(scatterCoeffs);
	importE->GetDimensions(dims);
	importE->GetDetectionGrid(detGrid);


	if(interface->GetPrintImport()){
		importE->PrintImport();

	}



	// Emitter

	std::vector<float> point;
	std::vector<float> rri;
	std::vector<float> angles;
	std::vector<float> eHeading;

	importE->GetEmitterAngles(angles);
	importE->GetEmitterRelativeRadiantIntensities(rri);
	importE->GetEmitterPointInPlane(point);
	importE->GetEmitterHeading(eHeading);

	//cout << "Emitter" << endl;
	Emitter *emitter = new Emitter;
	emitter->SetVerbosity(interface->GetVerbose());

	emitter->SetNumberOfPhotons(importE->GetNumberOfPhotons());
	emitter->SetPointInPlane(point);
	emitter->SetDistanceToSurface(importE->GetEmitterdistanceToSurface());
	emitter->SetMesh(importE->GetMeshImage());
	emitter->SetNormals(importE->GetNormals());
	emitter->SetRadiationCharacteristics(angles,rri);
	emitter->SetEmitterHeading(eHeading);

	std::vector<float> startPosition;
	std::vector<float> startVector;
	std::vector<int> startVectorOcc;
	std::vector<float> eX,eY,eZ;

	emitter->GetEmittermodel(startPosition,startVector,startVectorOcc);
	emitter->GetBasisVectors(eX,eY,eZ);
	if(interface->GetPrintEmitter()){
		emitter->PrintEmitter();

	}



if(interface->GetPerformSimulation()){

   // Calculate Light Propagation

   //Init
     CalculatePropagation* calcProp = new CalculatePropagation;
     calcProp->SetVerbosity(interface->GetVerbose());
     calcProp->SetShowPercentDone(interface->GetTimeSpent());

     calcProp->SetDeviceType(interface->GetDeviceType());

   // Set tissue properties
   calcProp->SetTissueProperties(refractionInd, absorbtionCoeffs,  anisotropies, scatterCoeffs);
   calcProp->SetFetalTissueID(importE->GetFetalID());
   // Set Model geometry
   calcProp->SetDimensions(dims);
   //calcProp->SetDetectionGridDimensions(1,1,1,1,1,1);
   calcProp->SetDetectionGridDimensions(detGrid);

   calcProp->SetTriangles(importE->GetMeshImage());
   calcProp->SetNormals(importE->GetNormals());
   calcProp->SetInnerAndOuterTissueIndices(importE->GetMediumIDs());
   calcProp->SetPixelPerUnitFactor(importE->GetPixelPerUnitFactor());
     // Set Simulation properties
   calcProp->SetStartTrajectories(startVector,startVectorOcc); // SetStartTrajectory vn -0.707107, 0.000000, 0.707107
   calcProp->SetStartPoint(startPosition); // SetStartPoint
   calcProp->SetPOI(point);
   calcProp->SetBasisVectors(eX,eY,eZ);

   calcProp->SetChanceForRoulette(importE->GetChance()); // SetChanceForRoulette
   calcProp->SetWeightThreshold(importE->GetWeightThreshold()); // SetWeightThreshold
   calcProp->SetNumberOfPhotons(importE->GetNumberOfPhotons()); // abs number of Photons
   calcProp->SetNumberOfPhotonInteractions(importE->GetNumberOfPhotonInteractions());

   calcProp->Init();

 if(interface->GetTimeSpent()){
   calcProp->PrintSimulationParams();
   }
    // Execute

   calcProp->Propagation();

    // Get Logs
 //  calcProp->GetLogs(RandNmbersLog,RecordAxz,RecordAyz,SurfFetalTransRxy);


  // calcProp->GetAbsorbtionLogCoords(RAxzCoords,RAyzCoords);
 //  calcProp->GetAbsorbtionLog(RecordAxz, RecordAyz);

   if(interface->GetTimeSpent()){

	   calcProp->PrintTimeSpent();
   }

if(interface->GetPerformAnalysis()){

   // Analysis

   std::vector<int> surfaceIDs; // IDs der Cells | size: numTris
   std::vector<std::string> surfaceIDnames; // ID aus object file | size: numTris
   std::vector<std::string> surfaceIDnamesConfig;
   importE->GetSurfaceIDVector(surfaceIDs);
   importE->GetSurfaceNamesVector(surfaceIDnames);
   importE->GetSurfaceConfigNames(surfaceIDnamesConfig);


   Analysis *analyse = new Analysis;
   analyse->SetVerbosity(interface->GetVerbose());
   //analyse->SetVerbosityFalse();
/*
  // analyse->SetMesh(importE->GetMeshImage());
  // analyse->SetNormals(importE->GetNormals());
   //analyse->SetRandomArray(positionLog);
   //analyse->SetAbsorbtionLog(absorbtionLog);
   analyse->SetSurfaceIDs(surfaceIDs);
   analyse->SetSurfaceIDNames(surfaceIDnames);
   analyse->SetSurfaceIDNamesConfig(surfaceIDnamesConfig);
   analyse->SetNumberOfPhotons(importE->GetNumberOfPhotons());


   analyse->CalcHistogram(100);
   analyse->CalcChiSquared();
   analyse->CalculateSurfaceSums();
*/
   analyse->SetNumberOfPhotons(importE->GetNumberOfPhotons());
   analyse->SetRTdxyLogs(calcProp->GetSFTImage());

   analyse->CalculateSurfaceSums();

   //analyse->PrintChiSquared();
   analyse->PrintAnalysis();


   // write files
if(importE->GetSaveFile()){

	std::string saveDir = importE->GetSaveDir();

	cout << "Writing Files in " << saveDir << "\n"  << endl;

	Logs *log = new Logs;
	log->SetSaveDir(saveDir);
	log->SetDimensions(dims);
	log->SetGridDimensions(detGrid);
	log->SetNumberOfPhotons(importE->GetNumberOfPhotons());
	log->SetRTdxyLogs(calcProp->GetRTdImage());

	log->SetSFTLogs(calcProp->GetSFTImage());
	log->SetPixelPerUnitFactor(importE->GetPixelPerUnitFactor());
	//log->SetAxzLog(RecordAxz);
	//log->SetAyzLog(RecordAyz);
	//log->SetAxzCoords(RAxzCoords);
	//log->SetAyzCoords(RAyzCoords);
	log->SetScaleString(importE->GetScaleString());
	log->SetColormapString(importE->GetColormapString());
	log->SetCoordSystemString(importE->GetCoordSystemString());

	int dimx = int(dims[0]/detGrid[0]);
	int dimy = int(dims[1]/detGrid[1]);
	int dimz = int(dims[2]/detGrid[2]);


//	std::vector<int> picDims;
//	importE->GetPictureDimensions(picDims);

	//int R = picDims[3];


	//log->WriteArrInFile("RandNmbersLog.txt",RandNmbersLog);
//	std::cout << "Save Absorption Logs..." << endl;
//	log->SaveAbsorbtionLog(importE->GetRecordAxzFilename(),RecordAxz,RAxzCoords);
//	log->SaveAbsorbtionLog(importE->GetRecordAyzFilename(),RecordAyz,RAyzCoords);
//	log->SaveSurfFetTransArr(importE->GetSurfaceFilename(),importE->GetFetalFilename(),importE->GetTransmissionFilename());


	log->SaveSurfaceArr(importE->GetSurfaceFilename(),calcProp->GetSFTImage(Rxy_Log_SURF));

	log->SaveRecordRd(importE->GetSurfaceFilename(),calcProp->GetRTdImage(Rxy_Log_SURF));
//	log->SaveSurfaceArr(importE->GetFetalFilename(),calcProp->GetSFTImage(Rxy_Log_FETAL));
//	log->SaveSurfaceArr(importE->GetTransmissionFilename(),calcProp->GetSFTImage(Rxy_Log_TRANS));

	log->CalcReflectancePerRadius(importE->GetSurfaceFilename(),calcProp->GetRTdImage(Rxy_Log_SURF));
	log->CalcReflectancePerAngle(importE->GetSurfaceFilename(),calcProp->GetRTdImage(Rxy_Log_SURF));

	log->CalcReflectancePerRadius(importE->GetTransmissionFilename(),calcProp->GetRTdImage(Rxy_Log_TRANS));
	log->CalcReflectancePerAngle(importE->GetTransmissionFilename(),calcProp->GetRTdImage(Rxy_Log_TRANS));

	log->CalcReflectancePerRadius(importE->GetFetalFilename(),calcProp->GetRTdImage(Rxy_Log_FETAL));
	log->CalcReflectancePerAngle(importE->GetFetalFilename(),calcProp->GetRTdImage(Rxy_Log_FETAL));

	log->FluenceZ(importE->GetRecordAyzFilename(),calcProp->GetPhirzImage());


std::cout << "Create Pictures..." << endl;

	log->CreateSFT_XYPicture(importE->GetSurfaceFilename(),calcProp->GetSFTImage(Rxy_Log_SURF));
	log->CreateSFT_XYPicture(importE->GetTransmissionFilename(),calcProp->GetSFTImage(Rxy_Log_TRANS));
	log->CreateSFT_XYPicture(importE->GetFetalFilename(),calcProp->GetSFTImage(Rxy_Log_FETAL));



	log->CreateAbsorbtionPicture(importE->GetRecordAxzFilename(),calcProp->GetArzImage());
	log->CreateAbsorbtionPicture(importE->GetRecordAyzFilename(),calcProp->GetPhirzImage());





   cout << "Writing Done." << endl;

   SurfFetalTransRxy.clear();

   RecordAxz.clear();
   RecordAyz.clear();
   RandNmbersLog.clear();
   //hist.clear();
   delete log;
}
}


}

   // Graphics output
if(interface->GetGraphicsOut()){



   vtkSmartPointer<vtkCamera> cam = vtkSmartPointer<vtkCamera>::New();

    cam->SetPosition(5/2, -5*3,5/2);
    cam->SetFocalPoint(0, 0, 0);
     cam->SetViewUp(0, 0, -1);
     cam->Azimuth(importE->GetViewAzimuth()); // yaw
     cam->Elevation(importE->GetViewElevation()); // pitch
     cam->Roll(importE->GetViewRoll()); // roll
     cam->Modified();


  vtkRenderer *ren = vtkRenderer::New();

 // ren->AddViewProp(vm->GetVoxelVolume());
  ren->AddActor(importE->GetActor());
 ren->AddActor(importE->GetNormalArrowsActor());
  ren->SetBackground(0,1,1);
 ren->SetActiveCamera(cam);



  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  renWin->SetSize(600,600);
  //

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  iren->Initialize();
  iren->Start();


  ren->Delete();
  renWin->Delete();
  iren->Delete();


}

return EXIT_SUCCESS;
 //delete analyse;
 //delete calcProp;
 // delete emitter;
 //delete importE;
	} else{
		return EXIT_FAILURE;

	}


}

