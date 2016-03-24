/*
 * PSNR.hpp
 *
 *  Created on: Mar 9, 2015
 *      Author: matthias
 */

#ifndef PSNR_HPP_
#define PSNR_HPP_

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>


class PSNR {


	int inext,inextp;
	long ma[56];
	int iff;
	long mj,mk;
	int i,ii,k;
	long MBIG;
	long MSEED;
	long MZ;
	double FAC;
	long ran3Seed;

	int numberGenerator;



public:


	void SetNumbergenerator(int i);
	void InitRand3();
	void SetSeed(int seed);
	void SetMBig(long mbig);
	void SetMSeed(long mseed);
	void SetFac();
	void SetMZ(int mz);

	void GetRan3Array(std::vector<float> &arr);
	double GetRan3();
	void GetRan3(double &random);
	void GetRan3(std::vector<float> &randoms,int numRandoms);
	void GetInextInextp(std::vector<int> &vec);
	void GetInextInextp(std::vector<int> &vec,int numWorkGroups);
	void GetMBMZ(std::vector<long> &vec);
	

	PSNR();
	virtual ~PSNR();




private:

	


};

#endif /* PSNR_HPP_ */
