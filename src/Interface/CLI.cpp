/*
 * CLI.cpp
 *
 *  Created on: Nov 18, 2015
 *      Author: matthias
 */

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "CLI.hpp"



CLI::CLI() {
	// TODO Auto-generated constructor stub


	//argv = NULL;
	argcount = 0;
	configFilename = "config.xml";
	verbose = false;
	graphics = false;
	printImport = false;
	printEmitter = false;
	printTimeSpent = false;
	performSimulation = false;
	analysis = false;
	deviceType = 0; // CPU
	progName.assign("");



}

CLI::~CLI() {
	// TODO Auto-generated destructor stub
}

int CLI::ParseArgs(int argc, char* argvs[]){

	if(argc<=1){
		PrintHelp();
		return EXIT_FAILURE;
	}else{
	//	std::cout << "argc: " <<argc <<std::endl;
	//std::cout << "argv(0): " << argvs[0] << std::endl;
		for(int i=1;i<argc;i++){
		//std::cout << "argvs("<<i<<"): " << argvs[i] << std::endl;

			if(strcmp("-h", argvs[i]) == 0){
				PrintHelp();

				return EXIT_FAILURE;
			}else if(strcmp("-f", argvs[i]) == 0){
			//std::cout << "config:"<< argvs[i+1] << std::endl;
				SetConfigFile(argvs[i+1]);
				i++;
			}else if(strcmp("-v", argvs[i]) == 0){
				SetVerbose(true);
			}else if(strcmp("-g",argvs[i])==0){
				SetGraphicsOutput(true);
			}else if(strcmp("-I",argvs[i])==0){
			//std::cout << argvs[i] << std::endl;
				SetPrintImport(true);
			}else if(strcmp("-E",argvs[i])==0){
				//std::cout << argvs[i] << std::endl;
				SetPrintEmitter(true);
			}else if(strcmp("-t",argvs[i])==0){
				SetPrintTimeSpent(true);
			}else if(strcmp("-s",argvs[i])==0){
				SetPerformSimulation(true);
			}else if(strcmp("-a",argvs[i])==0){
				SetPerformAnalysis(true);
			}else if(strcmp("-gpu",argvs[i])==0){
				SetDeviceType(1);
			}else if(strcmp("-cpu",argvs[i])==0){
				SetDeviceType(0);
			}else if((strcmp("-device",argvs[i])==0 )|| (strcmp("-d",argvs[i])==0)){

				if(strcmp("gpu",argvs[i+1])==0){
					SetDeviceType(1);
				}else if(strcmp("cpu",argvs[i+1])==0){
					SetDeviceType(0);
				}

			}else{
				PrintHelp();
				return EXIT_FAILURE;
			}

		}
	}
return EXIT_SUCCESS;

}

void CLI::GetConfigFilename(std::string &filename){

	filename.assign(configFilename);

}
const char* CLI::GetConfigFilename(){
	return configFilename.c_str();
}


bool CLI::GetVerbose(){
	return verbose;
}
bool CLI::GetGraphicsOut(){
	return graphics;
}
bool CLI::GetPrintImport(){
	return printImport;
}
bool CLI::GetPrintEmitter(){
	return printEmitter;
}
bool CLI::GetTimeSpent(){
	return printTimeSpent;
}
bool CLI::GetPerformSimulation(){
	return performSimulation;
}
bool CLI::GetPerformAnalysis(){
	return analysis;
}

int CLI::GetDeviceType(){
	return deviceType;
}


void CLI::SetConfigFile(char* filename){
	configFilename.assign(filename);

}

void CLI::SetVerbose(bool verb){
	verbose = verb;
}
void CLI::SetGraphicsOutput(bool gr){
	graphics = gr;
}
void CLI::SetPrintImport(bool pImport){
	printImport = pImport;
}
void CLI::SetPrintEmitter(bool pEmitter){
	printEmitter = pEmitter;

}
void CLI::SetPrintTimeSpent(bool pTimeSpent){
	printTimeSpent = pTimeSpent;

}
void CLI::SetPerformSimulation(bool perfSim){
	performSimulation = perfSim;
}
void CLI::SetPerformAnalysis(bool perfAn){
	analysis = perfAn;

}
void CLI::SetDeviceType(int devType){
// 0: GPU ; 1: CPU

	deviceType = devType;
}


void CLI::PrintAbout(){

	std::cout << "\n################################################################" << std::endl;
	std::cout << "Monte Carlo Simulation of Light propagation in biological Tissue" << std::endl;
	std::cout << "Author: Matthias Minx" << std::endl;
	std::cout << "#################################################################" << std::endl;




}
void CLI::PrintHelp(){

	PrintAbout();
	std::cout << "Usage: \n" << this->progName.c_str() << "\n -[hfvgIEtsad] [options]" << std::endl;
	PrintHelpConfigFile();
	PrintHelpVerbose();
	PrintHelpGraphicsOut();
	PrintHelpPrintImport();
	PrintHelpPrintEmitter();
	PrintHelpTimeSpent();
	PrintHelpPerformSimulation();
	PrintHelpPerformAnalysis();
	PrintHelpDeviceType();

	std::cout << "#################################################################" << std::endl;


}
void CLI::PrintHelpHelp(){
	std::cout << "-h \t this help" << std::endl;
}
void CLI::PrintHelpConfigFile(){
	std::cout << "-f filename \t choose xml-configuration file" << std::endl;

}
void CLI::PrintHelpVerbose(){
	std::cout << "-v \t Print more infos during program run. " << std::endl;
}
void CLI::PrintHelpGraphicsOut(){
	std::cout << "-g \t Show simulation model." << std::endl;
}
void CLI::PrintHelpPrintImport(){
	std::cout << "-I \t Print MediumIDs, Cell normals and cells" << std::endl;
}
void CLI::PrintHelpPrintEmitter(){
	std::cout << "-E \t Print emitter: vStart, pStart and POI " << std::endl;
}
void CLI::PrintHelpTimeSpent(){
	std::cout << "-t \t Print calculation time." << std::endl;
}
void CLI::PrintHelpPerformSimulation(){
	std::cout << "-s \t Perform Simulation" << std::endl;
}
void CLI::PrintHelpPerformAnalysis(){
	std::cout << "-a \t Perform Analysis: calculate Surface sums" << std::endl;
}
void CLI::PrintHelpDeviceType(){
	std::cout << "-d gpu|cpu \t Choose OpenCL-Device"<< std::endl;
	std::cout << "-device gpu|cpu" << std::endl;
	std::cout << "-gpu | -cpu" << std::endl;
}
