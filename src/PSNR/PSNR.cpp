/*
 * PSNR.cpp
 *
 *  Created on: Mar 9, 2015
 *      Author: matthias
 */

#include "PSNR.hpp"



PSNR::PSNR() {


	numberGenerator = 0;

	inext = 0;
	inextp = 0;

	iff = 0;
	mj = 0;
	mk = 0;
	i = 0;
	ii = 0;
	k = 0;
	MBIG = 1000000000;
	MSEED = 161803398;
	MZ = 0;
	FAC = 1.0/MBIG;
	ran3Seed =-(int)time(NULL)%(1<<15);



}

PSNR::~PSNR() {


}



void PSNR::SetNumbergenerator(int i){

	this->numberGenerator = i;
}

void PSNR::InitRand3(){

	if (ran3Seed < 0 || iff == 0) {
			iff=1;
			mj=MSEED-(ran3Seed < 0 ? -ran3Seed : ran3Seed);
			mj %= MBIG;
			ma[55]=mj;
			mk=1;
			for (i=1;i<=54;i++) {

				ii=(21*i) % 55;
				ma[ii]=mk;
				mk=mj-mk;
				if (mk < MZ) mk += MBIG;
				mj=ma[ii];
			}
			for (k=1;k<=4;k++)
				for (i=1;i<=55;i++) {
					ma[i] -= ma[1+(i+30) % 55];
					if (ma[i] < MZ) ma[i] += MBIG;
				}
			inext=0;
			inextp=31;
			ran3Seed=1;
		}
		if (++inext == 56) inext=1;
		if (++inextp == 56) inextp=1;
		mj=ma[inext]-ma[inextp];
		if (mj < MZ) mj += MBIG;
		ma[inext]=mj;


}


void PSNR::SetSeed(int seed){

	this->ran3Seed = -abs(seed);

}
void PSNR::SetMBig(long mbig){
	this->MBIG = mbig;

}
void PSNR::SetMSeed(long mseed){
	this->MSEED = mseed;

}
void PSNR::SetFac(){
	this->FAC = 1/(this->MBIG);

}
void PSNR::SetMZ(int mz){
	this->MZ = mz;

}

void PSNR::GetRan3Array(std::vector<float> &arr){

	std::vector<float> a;
	a.assign(56,0);

	for(int i = 0;i<56;i++){

		a[i] = float(ma[i]);
	}

	arr.assign(a.begin(),a.end());

}

double PSNR::GetRan3(){

	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;

	return mj*FAC;
}

void PSNR::GetRan3(double &random){

	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;

	random = mj*FAC;
}

void PSNR::GetRan3(std::vector<float> &randoms,int numWorkGroups){

	randoms.assign(56*numWorkGroups,0);

	for(int i=0;i<56*numWorkGroups;i++){

		if (++inext == 56) inext=1;
		if (++inextp == 56) inextp=1;
			mj=ma[inext]-ma[inextp];
		if (mj < MZ) mj += MBIG;
			ma[inext]=mj;

		randoms[i] = float(mj);

	}
}

void PSNR::GetInextInextp(std::vector<int> &vec){

	std::vector<int> v;
	v.assign(3,0);

	v[0] = inext;
	v[1] = inextp;
	v[2] = 0;

	vec.assign(v.begin(),v.end());

}
void PSNR::GetInextInextp(std::vector<int> &vec,int numWorkGroups){

	std::vector<int> v;
	v.assign(3*numWorkGroups,0);

	for(int i=0;i<56*numWorkGroups;i++){

		v.push_back(inext);
		v.push_back(inextp);
		v.push_back(0);

	}
	vec.assign(v.begin(),v.end());

}
void PSNR::GetMBMZ(std::vector<long> &vec){

	std::vector<long> v;
	v.assign(2,0);

		v[0] = MBIG;
		v[1] = MZ;

	vec.assign(v.begin(),v.end());

}
