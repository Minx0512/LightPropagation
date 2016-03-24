/*
 * CLI.hpp
 *
 *  Created on: Nov 18, 2015
 *      Author: matthias
 */

#ifndef CLI_HPP_
#define CLI_HPP_

#include <stdlib.h>
#include <string.h>

class CLI {

	char *argv[];
	int argcount;
	std::string progName;
	std::string configFilename;
	bool verbose;
	bool graphics;
	bool printImport;
	bool printEmitter;
	bool printTimeSpent;
	bool performSimulation;
	bool analysis;
	int deviceType; // GPU | CPU


public:

	const char* GetConfigFilename();
	void GetConfigFilename(std::string &filename);
	bool GetVerbose();
	bool GetGraphicsOut();
	bool GetPrintImport();
	bool GetPrintEmitter();
	bool GetTimeSpent();
	bool GetPerformSimulation();
	bool GetPerformAnalysis();
	int GetDeviceType();


	void PrintAbout();
	void PrintHelp();
	int ParseArgs(int argc,char* argvs[]);


	CLI();
	virtual ~CLI();

private:

	void SetConfigFile(char* filename);
	void SetVerbose(bool verb);
	void SetGraphicsOutput(bool gr);
	void SetPrintImport(bool pImport);
	void SetPrintEmitter(bool pEmitter);
	void SetPrintTimeSpent(bool pTimeSpent);
	void SetPerformSimulation(bool perfSim);
	void SetPerformAnalysis(bool perfAn);
	void SetDeviceType(int devType);

	void PrintHelpHelp();
	void PrintHelpConfigFile();
	void PrintHelpVerbose();
	void PrintHelpGraphicsOut();
	void PrintHelpPrintImport();
	void PrintHelpPrintEmitter();
	void PrintHelpTimeSpent();
	void PrintHelpPerformSimulation();
	void PrintHelpPerformAnalysis();
	void PrintHelpDeviceType();



};



#endif /* CLI_HPP_ */
