typedef struct iSFT{
	int i;
	float3 w;
}iSFT;



void atomic_add_global(volatile global float *source, const float operand) {
    union {
        unsigned int intVal;
        float floatVal;
    } newVal;
    union {
        unsigned int intVal;
        float floatVal;
    } prevVal;
 
    do {
        prevVal.floatVal = *source;
        newVal.floatVal = prevVal.floatVal + operand;
    } while (atomic_cmpxchg((volatile global unsigned int *)source, prevVal.intVal, newVal.intVal) != prevVal.intVal);
}

float GetDistanceToSurface(float3 point, float3 P0, float3 e_z){
	
	float d = 0.0f;
	float3 normal = normalize(e_z);
	float3 PP0 = point-P0;
	 d = fabs(dot(normal,PP0));
	//printf(" DistToSurf: %f  e_z: %f , %f , %f | P-P0: %f , %f , %f \n", d,e_z.x, e_z.y, e_z.z,PP0.x, PP0.y,PP0.z);
	return d;	
}

float CalcCartesianCoordPerAxis(float3 point, float3 P0, float3 axis){
	
	float d = 0.0f;
	float3 normal = normalize(axis);
	float3 PP0 = point-P0;
	 d = dot(normal,PP0);
	//printf(" DistToSurf: %f  e_z: %f , %f , %f | P-P0: %f , %f , %f \n", d,e_z.x, e_z.y, e_z.z,PP0.x, PP0.y,PP0.z);
	return d;	
}

float CalcRadius(float3 position, float3 POI,float3 e_Z){
	// POI: point of incident
	// e_X, e_Y, e_Z: Einheits vectors
	// exit coordinates
	// use polar coords --> gridDims[3,4]: r, phi = gridDims.x, gridDims.y
	// 
//	float3 n1 = cross(POI-position,e_Z);	
	float3 n = POI-position;

	float r = sqrt(n.x*n.x +  n.y*n.y+n.z*n.z);

//printf("\t\tRadius | coords: %f,%f,%f  , r = %f || coords1: %f,%f,%f  , r1 = %f  | r2(n1)= %f \n",n.x,n.y,n.z, length(n) , n1.x,n1.y,n1.z,length(n1),r );	
	

return r;

}
float CalcAlphaOut(float3 photonV,float3 surfNormal){

	return  fabs(acos(dot(photonV,surfNormal)));

}
int CalcIdx(float3 position, float3 POI,float3 e_X, float3 e_Y, float3 e_Z, float3 gridDims){
	// POI: point of incident
	// e_X, e_Y, e_Z: Einheits vectors
	// exit coordinates
	// use polar coords --> gridDims[3,4]: r, phi = gridDims.x, gridDims.y
	// 
	float3 n = cross(position-POI,e_Z);	
	float radius = length(n);		
	float3 vec = normalize(cross(e_Z,n));
		
	float cosPhix = dot(vec,e_X);
	float cosPhiy = dot(vec,e_Y);
	
	if(cosPhix>1){cosPhix=1;}
	if(cosPhiy>1){cosPhiy=1;}	
	float phi = 0.0f;
	 if(cosPhiy<0){
		 // 3. oder 4. Quadrant
		 phi = 360-degrees(acos(cosPhix));
	 }else{
		phi = degrees(acos(cosPhix));
	 }
	int r = int(radius/gridDims.x);	
	int phiBin = int(phi/gridDims.y);
	int idx = 0;
	 if(r>0){		 
		 idx = r*360/gridDims.y+phiBin;		 
	}
	
	
	
	return idx;
	
}

float3 CalcCoords(float3 position, float3 POI,float3 e_X, float3 e_Y, float3 e_Z){
	// POI: point of incident
	// e_X, e_Y, e_Z: Einheits vectors
	// exit coordinates
	// use polar coords --> gridDims[3,4]: r, phi
	// 
	float3 coords;
	coords.x = CalcCartesianCoordPerAxis(position, POI, e_X); // x-Coordinate	
	coords.y = CalcCartesianCoordPerAxis(position, POI, e_Y); // y-Coordinate	
	coords.z = CalcCartesianCoordPerAxis(position, POI, e_Z); // z-Coordinate
	
	
	return coords;
	
}



void RecordFSTxy(__global float* SurfFetalTransRxy, float3 SFTxyWeight, float3 coords, float alphaOut, float radius){
	
		
		int photonID = get_global_id(0);
		
	
		
		vstore8((float8)(coords,SFTxyWeight,alphaOut,radius),photonID,SurfFetalTransRxy);
		
	//printf("idx: %i | Fetal %f , %f | Surface %f , %f | Trans %f , %f \n", idx, weight.x, FetalRxy[idx], weight.y, SurfRxy[idx],weight.z, TransTxy[idx]);
	
}

void RecordFSTxy(__global float* SurfFetalTransRxy, float3 SFTxyWeight, int idx){
	
		
		int photonID = get_global_id(0);
		
	
		
		vstore4((float4)(idx,SFTxyWeight),photonID,SurfFetalTransRxy);
		
	//printf("idx: %i | Fetal %f , %f | Surface %f , %f | Trans %f , %f \n", idx, weight.x, FetalRxy[idx], weight.y, SurfRxy[idx],weight.z, TransTxy[idx]);
	
}
void RecordFSTxy(__local float* SurfFetalTransRxy, float3 SFTxyWeight, int idx, int count){
	
		
		int photonID = get_global_id(0);
		
	if(SFTxyWeight.x>0 || SFTxyWeight.z>0){		
		//SurfFetalTransRxy[count] = (float4)(idx,SFTxyWeight);
		vstore4((float4)(idx,SFTxyWeight),count,SurfFetalTransRxy);		
	}
	//printf("idx: %i | Fetal %f , %f | Surface %f , %f | Trans %f , %f \n", idx, weight.x, FetalRxy[idx], weight.y, SurfRxy[idx],weight.z, TransTxy[idx]);
	
}
void RecordFSTxy(__local iSFT* SurfFetalTransRxy, float3 SFTxyWeight, int idx, int count){			
		int photonID = get_global_id(0);		
	if(SFTxyWeight.x>0 || SFTxyWeight.z>0){
		
		SurfFetalTransRxy[count].i = idx;
		SurfFetalTransRxy[count].w = SFTxyWeight;
		//vstore4((float4)(idx,SFTxyWeight),photonID,SurfFetalTransRxy);
		
	}
	//printf("idx: %i | Fetal %f , %f | Surface %f , %f | Trans %f , %f \n", idx, weight.x, FetalRxy[idx], weight.y, SurfRxy[idx],weight.z, TransTxy[idx]);
	
}


void RecordAxz(__global float* Axz, float absorbedWeight, float2 coords,float3 dimsXYZ){
	
	
	int idx = (int) (coords.x*dimsXYZ.x+coords.y);
	
	atomic_add_global(&Axz[idx], absorbedWeight);

//	vstore4((float4)(coords,absorbedWeight,0),photonID,AxzAyz);
	
}



void RecordAxz(__global float* Axz, float absorbedWeight, float2 coords){
	int photonID = get_global_id(0);
	
	vstore4((float4)(coords,absorbedWeight,0),photonID,Axz);
	
}


void RecordAyz(__global float* Ayz, float absorbedWeight, float2 coords,float3 dimsXYZ){
	
	int idx = (int) (coords.x*dimsXYZ.y+coords.y);	
	atomic_add_global(&Ayz[idx], absorbedWeight);
	
}





void RecordAyz(__global float* AxzAyz, float absorbedWeight, float2 coords){
	int photonID = get_global_id(0);
	vstore4((float4)(coords,0,absorbedWeight),photonID,AxzAyz);
	
}

void RecordAbsorbtionXZYZ(__global float* AxzAyz, float absorbedWeight, float3 coords, float3 gridDimsXYZ, int3 dimsXYZ){
	// POI: point of incident
	// Plane xz and yz			
	//float3 pPlaneXZ = POI + e_X;
	//float3 pPlaneYZ = POI + e_Y;
	// Axz: size: dimx*dimz
	
		
	if((fabs(coords.y)<gridDimsXYZ.y/2) || (fabs(coords.x)<gridDimsXYZ.x/2)){
		int cz = int(fabs(coords.z/gridDimsXYZ.z));
					
				if(fabs(coords.y)<gridDimsXYZ.y/2){
					//	RecordAxz(AxzAyz, absorbedWeight, (float2)(coords.x,coords.z),dimsXYZ);
					int idx = (int) ((cz+0.5f)*dimsXYZ.x+coords.x+cz);
//printf("Axz: (cx,cz0,cz1) : dimsXYZ.x : (%f,%f,%i) : %i : %i \n",coords.x,coords.z,cz,dimsXYZ.x, 2*idx);

					atomic_add_global(&AxzAyz[2*idx],absorbedWeight);
					//float2 axz = vload2(idx,AxzAyz)+(float2)(absorbedWeight,0);
				
					//vstore2(axz,idx,AxzAyz);
			
		
				}

			if(fabs(coords.x)<gridDimsXYZ.x/2){
				//	RecordAyz(AxzAyz, absorbedWeight, (float2)(coords.y,coords.z), dimsXYZ);
				int idx = (int) ((cz+0.5f)*dimsXYZ.y+coords.y+cz);	
//printf("Ayz: (cy,cz0,cz1) : dimsXYZ.y : (%f,%f,%i) : %i : %i \n",coords.y,coords.z,cz,dimsXYZ.y, 2*idx+1);

				atomic_add_global(&AxzAyz[2*idx+1],absorbedWeight);
				//float2 ayz = vload2(idx,AxzAyz)+(float2)(0,absorbedWeight);
				//vstore2(ayz,idx,AxzAyz);
			
			}

//printf("absorbedWeight: %f \n",absorbedWeight);
	
	
	}
	
	
//printf("absorbedWeight: %f \n",absorbedWeight);		
	
	
	
}
void RecordAbsorbtionRZ(__global float* Arz, float absorbedWeight, float3 coords,float2 gridDimsRZ, float3 dimXYZ, float mua){
	// POI: point of incident
	// Plane xz and yz			
	//float3 pPlaneXZ = POI + e_X;
	//float3 pPlaneYZ = POI + e_Y;
	// Axz: size: dimr*dimz
	// R (-r,+r)
		//int id = get_global_id(0);
		float radius = sqrt(coords.x*coords.x+coords.y*coords.y)*sign(coords.x)/gridDimsRZ.x;
		int dimR = int(ceil(sqrt(pow(dimXYZ.x/2,2)+pow(dimXYZ.y/2,2))));
	
		int id = int(coords.z*(2*dimR+1)/gridDimsRZ.y + radius+dimR);
		float3 tempAz = vload3(id,Arz);
		vstore3((float3)(1,tempAz.y+absorbedWeight,mua),id,Arz);
//	printf("%f, %f : %i | %f\n", radius,coords.z, id, absorbedWeight);		
	
	
}

float Rspecular(__global const float* nmuasg){
  float r1, r2;
	/* direct reflections from the 1st and 2nd layers. */
  float temp;
  float n1;
	n1 = nmuasg[0];//n0
	//n2 = nmuasg[4]; //n1
	//n3 = nmuasg[8]; //n2
	//mua = nmuasg[5]; //mua1
	//mus = nmuasg[6]; //mus1
	float4 n2muas;// n1, mua1, mus1, g1
	n2muas = vload4(1,nmuasg);
	
  temp =(n1 - n2muas.x)/(n1 + n2muas.x);
  r1 = pown(temp,2);
  
  if(n2muas.y == 0.0f && n2muas.z == 0.0f){
/* glass layer. */
		float n3;
		n3 = nmuasg[8]; //n2
    temp = (n2muas.x - n3)/(n2muas.x + n3);
    r2 = temp*temp;
    r1 = r1 + (1-r1)*(1-r1)*r2/(1-r1*r2);
  }
  
  return (r1);	
}



class PhotonClass {

	
	// float3 e_x,e_y,e_z; // Einheitsvektoren des Vektorraumes eines Photons
	float3 v; // normalized trajectory Vector of Photon
	float3 p; // Position of Photon
	float2 alpha; // alpha_in, alpha_out; //Angles of v to Plane in radians
	//float3 nv; //Normal vector of plane	
	//float	stepUnfinished; //, pathlength;
	float w; // photonweight / survival probability
	// pseudo number generator
	
	// Log
	float sumR, sumT;
	float3 SFTxy;
	// short currentTissueID;
	//	short nextTissueID;
	/*
	long rand3Array[56];
	int inext, inextp;
	long MBIG,MZ;
	float random;
	int numberGenerator;
	
	*/
	float3 evY;
	
	
public:
	float stepsize;
	float steptoboundary;
	float unfinStepSize;
	bool traversedFetalTissue;
	int refTri;
float absorbedWeight;
//float3 eYV;
//float3 eZV;

	PhotonClass(){};
	~PhotonClass(){};

void Init(float3 vstart, const float3 pstart, float Rsp){
		// void Init(/* float* n_t, float* mua_t, float* mus_t, float* g_t,*/__global float3 v_start, __global float3 p_start, /*float chance_survival, float w_th,*/ uint sd, short numGen){
	// this->currentTissueID = 0;
	// tissue properties	
	// this->nextTissueID = 1;
	
	// Photon properties
	// normal trajectory Vector of Photon packets
	
	// this->v = (float3) (v_start[0], v_start[1], v_start[2]);
	 this->v = vstart;
  //printf("\tin Init: vStart = %f, %f, %f\n",v.x,v.y,v.z);
	// Position of Photon packet
	//this->p = vload3(0,pstart);
	
	this->p = pstart;
	//this->p = (float3)(0,0,11.9);
	
	//printf("\t pStart = %f, %f, %f\n",p.x,p.y,p.z);
	
	this->w = 1.0f - Rsp; // photonpacket weight / survival probability
		
	//Einheitsvektoren des Vektorraumes eines Photons
	//this->e_x = (float3)(1,0,0);
	//this->e_y = (float3)(0,1,0);
	//this->e_z = (float3)(0,0,1); 
	//Normal vector of plane
	//this->nv = (float3)(0,0,1);	
	//Angles of v to Plane in radians
	// this->alpha_out = 0; 
	//this->alpha = (float2) 0.0f;
	//this->alpha.x = 0;
	//this->alpha.y = 0;
	SFTxy = (float3) 0.0f;
	stepsize = 0.0f;
	unfinStepSize = 0.0f;
	traversedFetalTissue = false;
	// path lengths
	//this->stepUnfinished = 0;
	refTri = -1;
	absorbedWeight = 0.0f;
	evY = (float3)(0,1,0);
	//eYV = (float3)(0,1,0);
	//eZV = (float3)(0,0,1);
				
}

float3 RotTheta(__global const float *vStart,float theta, int vecGroup){
// vStart: {v0,v1,...} w vi = (xi,yi,zi) 

float3 vec;
float3 n; // Rot um theta um vector n = e_x

//n.x = vStart[0];
//n.y = vStart[1];
//n.z = vStart[2];
n = vload3(0,vStart);

//vec.x = vStart[3*vecGroup];
//vec.y = vStart[3*vecGroup+1];
//vec.z = vStart[3*vecGroup+2];
vec = vload3(vecGroup,vStart);
	
float cosTheta = cos(theta);
float sinTheta = sqrt(1-pown(cosTheta,2));

/*	
float3 nDx = dot(n,vec);
float3 nx = cross(n,vec);

float3 nxXn = cross(nx,n);
*/
float3 vNew;
vNew = Rotate(vec, n, cosTheta, sinTheta);

		// = n*nDx+cosTheta*nxXn + sinTheta*nx;
	//	vNew = Rotate(vec, n, cosTheta, sinTheta);
	
	

return vNew;

}


float3 GetvStart(__global const float* vStart,__global const int *vStartOcc, int id){

float3 vStartNew;
//int2 groupThetaID;
int s = vec_step(vStartOcc);		
int sum = 0;
//int vGroupID = 0;
if(s>1){
	
	for(int i = 0;i<s;i++){
		
		//int occ = vStartOcc[i]/360;

		sum+=vStartOcc[i];
		if(id<sum){
		//groupThetaID.x = i;
		
		int numTheta = id-sum+vStartOcc[i];
		float theta = numTheta*360/vStartOcc[i]; // in degrees		
		vStartNew = RotTheta(vStart,theta, i);


		break;
		}
	}

}else{

/*
vStartNew.x = vStart[0];
vStartNew.y = vStart[1];
vStartNew.z = vStart[2];
*/
	vStartNew = vload3(0,vStart);


}
return vStartNew;
}

float GetWeight(){
	
	return this->w;
	
}

void SetWeight(float weight){
	this->w = weight;
//	write_mem_fence(CLK_LOCAL_MEM_FENCE);	
//		write_mem_fence(CLK_GLOBAL_MEM_FENCE);
}

float3 GetSurfFetalTransxy(){return this->SFTxy;}

void SetHeadingVector(float3 vec){this->v = vec;}

float3 GetHeadingVector(){return this->v;}

float3 GetPhotonPosition(){return this->p;}
/*
void ExitVolume(bool partial){
// Photons exit volume	
	
	if(partial){
		// RecordT(w); // Record transmission
		
	}else{
		// RecordT(w); // Record transmission: RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
		this->w = 0;	// Photonpacket dead
	}
	
}

*/



void Transmit(float3 e_x, float dotprod, float n1, float n2){
/*		
* e_x = nv :  Normalenvektor der Ebene	
* dotprod: dot(nv,v)
* n1: Refraction index current medium
* n2: refraction index next medium
* v: Bewegungsvektor des Photons
*/
	//float3 e_x = normalize(nv);
	//float3 e_z = normalize(cross(nv,this->v));
	//float3 e_y = normalize(cross(e_x,e_z));
	float3 e_y = normalize(cross(e_x,cross(e_x,this->v)));
	float sinTransmit = sqrt(1-pown(dotprod,2))*n1/n2;
	float cosTransmit = -sqrt(1-pown(sinTransmit,2));
	

	if(dotprod>0.0f){	
		cosTransmit = fabs(cosTransmit);
	}
//printf("\tTransmit : v_old: %f , %f , %f ", this->v.x, this->v.y, this->v.z);
	this->v = normalize(cosTransmit*e_x - sinTransmit*e_y);
//printf("v_new: %f , %f , %f\n", this->v.x,this->v.y,this->v.z);

}


void Reflect(float3 e_x, float dotprod){
	/*		
		* nv: Normalenvektor der Ebene
		* v: Bewegungsvektor des Photons
		* alpha_out: Brechungswinkel bei ?bergang von Medium n1 zu n2
		* alpha_in: Einfallswinkel
	  */
	// float3 e_z = normalize(cross(e_x,this->v));
	float3 e_y = normalize(cross(e_x,cross(e_x,this->v)));
	//float cosA2 = dot(e_y,this->v);
	float sinRefl = sqrt(1-pown(dotprod,2));
	
//printf("\tReflect : v_old: %f , %f , %f | e_y: %f , %f , %f | cosA1: %f | sinRefl: %f | ", this->v.x, this->v.y, this->v.z, e_y.x,e_y.y, e_y.z, dotprod, sinRefl);
	//printf("\tReflect : v_old: %f , %f , %f ", this->v.x, this->v.y, this->v.z);
	//this->v = normalize(-dotprod*e_x + cosA2*e_y);
	this->v = -normalize(dotprod*e_x + sinRefl*e_y);
//printf("v_new: %f , %f , %f\n", this->v.x,this->v.y,this->v.z);
}


void DCScatter(float g,float randomTheta, float randomPsi ) {	
	

	float cosTheta = 0;
	float g1 = (1-pown(g,2))/(1-g+2*g*randomTheta);
	//g1*=g1;
	
		if(g!=0){			
			
			cosTheta =(1+pown(g,2)-pown(g1,2))/(2*g);
			if(cosTheta <-1){cosTheta = -1;}else if(cosTheta>1){cosTheta = 1;}
			
		}else {			
			cosTheta = 2*randomTheta-1;			
		}

	float sinTheta = sqrt(1.0f-pown(cosTheta,2));

	float psi = 2*M_PI_F*randomPsi;
	float cosPsi = cos(psi);

	float sinPsi = sqrt(1.0f-pown(cosPsi,2));

	if(psi>M_PI_F){
		sinPsi *= -1;
	}

	float3 vn = 0.00f;
 	if(fabs(this->v.z) ==1.00f)  { 	/* normal incident. */
   	 vn.x = sinTheta*cosPsi;
  	 vn.y = sinTheta*sinPsi;
  	 vn.z = cosTheta*sign(this->v.z);	
	  /* SIGN() is faster than division. */
  	}else  {		/* regular incident. */
    	float temp = sqrt(1.0f - this->v.z*this->v.z);

	vn.x = sinTheta*(this->v.x*this->v.z*cosPsi - this->v.y*sinPsi)/temp + this->v.x*cosTheta;
   	vn.y = sinTheta*(this->v.y*this->v.z*cosPsi + this->v.x*sinPsi)/temp + this->v.y*cosTheta;
	vn.z = -sinTheta*cosPsi*temp + this->v.z*cosTheta;

  	}

	this->v = normalize(vn);

}
void VScatter(float g,float randomTheta, float randomPsi ) {	
	//float3 nv = (float3)(0.0f,0.0f,1.0);	


	float3 e_x = normalize(this->v);
	float3 e_rot = 0.0f;
	if(this->v.z!=0.0f){
		float z = -(this->v.x + this->v.y)/(this->v.z);
		e_rot.x = 1;
		e_rot.y = 1;
		e_rot.z = z;

	}else if(this->v.y!=0.0f){
		float y = -(this->v.x)/(this->v.y);
		e_rot.x = 1;
		e_rot.y = y;
		e_rot.z = 0;
	}else if(this->v.x!=0.0f){
		//float x = -(this->v.y + this->v.z)/(this->v.x);
		e_rot.x = 0;
		e_rot.y = 1;
		e_rot.z = 0;
	}	

	e_rot = normalize(e_rot);

	float cosTheta = 0;

	float g1 = (1-pown(g,2))/(1-g+2*g*randomTheta);
		
	if(g!=0){			
			
		cosTheta =(1+pown(g,2)-pown(g1,2))/(2*g);
		if(cosTheta <-1){cosTheta = -1;}else if(cosTheta>1){cosTheta = 1;}
			
	}else {			
		cosTheta = 2*randomTheta-1;			
	}

	float sinTheta = sqrt(1.0f-pown(cosTheta,2));
	float psi = 2*M_PI_F*randomPsi;
	float cosPsi = cos(psi);
	float sinPsi = sqrt(1.0f-pown(cosPsi,2));

	if(psi>M_PI_F){
		sinPsi *= -1;
	}
	//printf("%f ", randomTheta);

	// 1. deflect: rotate theta [0,pi)
	float3 rot_theta = normalize(Rotate(this->v,e_rot,cosTheta, sinTheta));

	// 2. rotate psi [0,2*pi)
	this->v = normalize(Rotate(rot_theta,this->v,cosPsi, sinPsi));

//printf("%f %f\n",randomTheta,randomPsi); // print random numbers used


}


void QScatter(float g,float randomTheta, float randomPsi){	
	//float3 nv = (float3)(0.0f,0.0f,1.0);	

	float cosTheta = 0;
	float g1 = (1-pown(g,2))/(1-g+2*g*randomTheta);
	if(g!=0){

		cosTheta =(1+pown(g,2)-pown(g1,2))/(2*g);
		if(cosTheta <-1){cosTheta = -1;}else if(cosTheta>1){cosTheta = 1;}
			
	}else {			
		cosTheta = 2*randomTheta-1;			
	}


	float theta = acos(cosTheta);
	cosTheta = cos(theta/2);	
	float sinTheta = sqrt(1.0f-pown(cosTheta,2));

	float psi = 2*M_PI_F*randomPsi;

	float cosPsi = cos(psi/2);		
	float sinPsi = sqrt(1.0f-pown(cosPsi,2));
	if(psi/2>M_PI_F){
		sinPsi *= -1;
	}
	
	float3 eZ = (float3)(0,0,1);
	float3 eY = (float3)(0,1,0);
	float3 eX = (float3)(1,0,0);
	
	float4 qY = (float4)(eY*sinTheta,cosTheta);
	float4 qV = (float4)((this->v)*sinPsi,cosPsi);
	float4 qyz = QMult(qY,qV);
	
	this->v = QRotateVec(qyz,this->v);	
	
}

float3 Rotate(float3 vec, float3 axis, float cosAngle, float sinAngle){

	float nDx = dot(axis,vec);
	float3 nXx = cross(axis,vec);
	float3 nXxXn = cross(nXx,axis); 

	float3 nVec = 0.0f;
	
	nVec = axis*nDx + cosAngle*nXxXn + sinAngle*nXx;

return nVec;
}

float4 QuatFromEuler(float psi, float theta){
	
	float4 q;
	float4 qx,qy,qz;
	
	float3 eY = (float3)(0,1,0);
	float3 eZ = this->v;
	float3 eX = normalize(cross(eZ,eY));
	
	qx = (float4)((float3)(eX*sin(theta/2)),cos(theta/2));
	qy = (float4)((float3)(eY*sin(psi/2)),cos(psi/2));
	qz = (float4)((float3)(0,0,0),1);
	/*
	printf("eX: %f, %f, %f \n", eX.x,eX.y,eX.z);
	printf("eY: %f, %f, %f \n", eY.x,eY.y,eY.z);
	printf("eZ: %f, %f, %f \n", eZ.x,eZ.y,eZ.z);
	printf("qx: %f, %f, %f, %f \n", qx.x,qx.y,qx.z,qx.w);
	printf("qy: %f, %f, %f, %f \n", qy.x,qy.y,qy.z,qy.w);
	printf("qz: %f, %f, %f, %f \n", qz.x,qz.y,qz.z,qz.w);
	printf("\t Psi, Theta: %f,  %f \n ",psi*180/M_PI_F,theta*180/M_PI_F);
	*/
	
	float4 qyz=QMult(qy,qz);
	
//	printf("qyz: %f, %f, %f, %f \n", qyz.x,qyz.y,qyz.z,qyz.w);
	
	q = normalize(QMult(qx,qyz));
//	printf("q: %f, %f, %f, %f \n", q.x,q.y,q.z,q.w);

return (float4)(normalize(q));

}

float4 QuatFromEuler2(float yaw, float pitch, float roll){
	//float PIOVER180 = M_PI_F/180;

	float4 q;
	float p = pitch;
	float y = yaw;
	float r = roll; 

//	printf("\t Psi, Theta: %f,  %f \n ",y*180/M_PI_F,p*180/M_PI_F);
	
	float sinp = sin(p);
	float siny = sin(y);
	float sinr = sin(r);
	float cosp = cos(p);
	float cosy = cos(y);
	float cosr = cos(r); 
	q.x = sinr * cosp * cosy - cosr * sinp * siny;
	q.y = cosr * sinp * cosy + sinr * cosp * siny;
	q.z = cosr * cosp * siny - sinr * sinp * cosy;
	q.w = cosr * cosp * cosy + sinr * sinp * siny;


return (float4)(normalize(q));

}

float4 QConjugate(float4 q){


return (float4)(-q.x,-q.y,-q.z,q.w);
}

float4 QMult(float4 q1, float4 q2){
// Multiplying q1 with q2 applies the rotation q2 to q1
float3 v1,v2;
	v1 = (float3)(q1.xyz);
	v2 = (float3)(q2.xyz);
	
	float3 v = q1.w*v2 + q2.w*v1 + cross(v1,v2);
	float w = q1.w*q2.w - dot(v1,v2);

	//printf("QMult: q1:%f, %f, %f,%f | q2: %f, %f, %f,%f : %f,%f, %f,%f \n",q1.x,q1.y,q1.z,q1.w,q2.x,q2.y,q2.z,q2.w, v.x,v.y,v.z,w);
	
	
return (float4)((float3)v,w);
}



float3 QRotateVec(float4 quat, float3 vec){

	
	float3 nVec = normalize(vec);
	float4 nQuat= normalize(quat);
//printf("in RotateVec: nVec:  %f, %f, %f\n", nVec.x,nVec.y,nVec.z);
	float4 cQuat  = QConjugate(nQuat);
	float4 vecQ = (float4)((float3)nVec,0);
	//printf("in RotateVec: vecQ:  %f, %f, %f %f\n", vecQ.x,vecQ.y,vecQ.z, vecQ.w);
	
	float4 rQ = QMult(vecQ,cQuat);
	rQ = QMult(nQuat,rQ);


return normalize((float3)(rQ.xyz));
}





float3 Rotate2(float cosTheta, float sinTheta, float cosPsi, float sinPsi){

	float3 u = normalize(this->v);

//printf("\tu.x: %f , u.y: %f , u.z: %f |  v.x: %f , v.y: %f , v.z: %f \n",u.x,u.y,u.z, this->v.x,this->v.y,this->v.z);

	if(fabs(u.z)>=1){
		u.x = sinTheta*cosPsi;
		u.y = sinTheta*sinPsi;
		u.z = cosTheta*sign(u.z);

	}else{
	float temp = sqrt(1 - u.z*u.z);
    
		//x-direction
		u.x = sinTheta*(u.x * u.z * cosPsi - u.y * sinPsi)/temp + u.x * cosTheta;
		//y-direction
		u.y = sinTheta * (u.y * u.z * cosPsi + u.x * sinPsi)/temp + u.y * cosTheta;
    		//z-direction
    		u.z = -sinTheta * cosPsi * temp + u.z * cosTheta;
	}

return u;

}



void Absorb(float4 muas){	
	if(muas.y > 0.0f && muas.z > 0.0f){
		this->absorbedWeight = (this->w)*(muas.y/(muas.y + muas.z));
	this->w -= (this->absorbedWeight);
	}
	
}
	
void MoveStep(float stepS) {
		
	this->p += stepS*(this->v);
	this->refTri = -1;
		
}

	
float GetFresnel(float n1, float n2, float cosIn){
/*
 * cosIn: angle. 0<a1<90 degrees.
 */
	float r;
	
//printf("\t\t\t%s\n\t\t\t\tn1: %f | n2: %f\n\t\t\t\cosIn: %f\n","in GetFresnel",n1,n2,cosIn);
 if(n1==n2){			  	/** matched boundary. **/
  
    r = 0.0f;
  }else if(fabs(cosIn) == 1.0f){	/** normal incident. **/
    
    r = (n1-n2)/(n1+n2);
    r *= r;
  }else if(cosIn == 0.0f){	/** very slant. **/
   
    r = 1.0f;
  }else {			  		/** general. **/
    float sinIn, sinOut;  /* sine of the incident and transmission angles. */ 
    float cosOut;
    
    sinIn = sqrt(1-cosIn*cosIn);
    sinOut = n1*sinIn/n2;
    //printf("cosIn: %f |sinIn: %f | sinOut: %f\n",cosIn, sinIn, sinOut);
	if(sinOut>=1.0f) {	
	  /* double check for total internal reflection. */
     
     	 r = 1.0f;
   	}else{      
     		cosOut = sqrt(1-pown(sinOut,2));
		
	//printf("cosOut: %f | ",cosOut);
			
		float cap = cosIn*cosOut - sinIn*sinOut; /* c+ = cc - ss. */
	  float cam = cosIn*cosOut + sinIn*sinOut; /* c- = cc + ss. */
	  float sap = sinIn*cosOut + cosIn*sinOut; /* s+ = sc + cs. */
	  float sam = sinIn*cosOut - cosIn*sinOut; /* s- = sc - cs. */
	//printf("cap: %f , cam: %f , sap : %f , sam: %f \n",cap,cam, sap,sam);
			
	      r = 0.5f*pown(sam,2)*(pown(cam,2)+pown(cap,2))/(pown(sap,2)*pown(cam,2)); 
	
			
	}
  }


//printf("in GetFresnel \t r: %f\n",r);
		
	
return r;
}

float GetCosAlphaCrit(float n1,float n2){
// cos_crit = 0 -> 90? 
	float cos_crit = 0.0f;

	if(n1>n2){
		cos_crit = sqrt(1.0f - pown(n2/n1,2));
	}else{
		cos_crit = 0.0f;
	} 
 return cos_crit;
}

int TransmitOrReflect(float4 nMuasGCurr,float4 nMuasGNext,float3 POI, float3 eZ, const float distMin, float3 nv, int currentTissueID, int2 innerOuterID, float random, int triID, float Rsp){
//int TransmitOrReflect(__global float* n, __global float* mua,__global float* mus,float3 POI, float3 eZ, float distMin, bool traversedFetalTissue, int photonID,  float3 nv, int currentTissueID, int outerID, int innerID, float random){
 /* n: Gewebeeigenschaft: Brechungszahlen der Media
	* nv: normal einheits vector
	* random: random number
	*/
//printf("ID: %i\t", get_global_id(0));	
//printf("%s\n\tcurrTissueID: %i\n\touterID: %i\n\tinnerID: %i\n\trandom: %f\n", "in TransmitOrReflect\n",currentTissueID,innerOuterID.y,innerOuterID.x,random);
	//SetAlphaIn(nv);
	int newTissueID; 
	float2 muasNew;
	float r = 0.0f;
	float cos_crit; 
	float cosIn = dot(nv,this->v);
//printf("\tdotprod: %f\n",cosIn); // cos of incident angle
	float absCosIn = fabs(cosIn);
	float n1,n2; //n1: current n | n2: next n
	n1 = nMuasGCurr.x;
	n2 = nMuasGNext.x;
	
	
	//SetAlphaOut_Transmission(n1,n2,this->alpha.x);
	cos_crit = GetCosAlphaCrit(n1,n2);

	/* Get r. */

 	if( absCosIn <= cos_crit){ 
		r=1.0f;

	 }else{
		r = GetFresnel(n1,n2,cosIn);

	}
//printf("\tGetFresnel: %f | random: %f \n", r, random);
	if(random>r){
		//transmit
//printf("\t%s\n","random > r : Transmit");
		//SetAlphaOut_Transmission(n1,n2,this->alpha.x);	
		Transmit(nv, cosIn,n1,n2);
		
		if(cosIn>0){
		 // nv and v same direction
			newTissueID = innerOuterID.y;
//printf("\t\t%s %i\n","dotprod > 0 => newTissueID = outerID:", innerOuterID.y);
		}else{
		// nv and v opposed directions
			newTissueID = innerOuterID.x;
//printf("\t\t%s %i\n","dotprod < 0 => newTissueID = innerID:", innerOuterID.x);
		}
	
		muasNew.x = nMuasGNext.y;
		muasNew.y = nMuasGNext.z;
	}else{ 
//printf("\t%s\n","random < r : Reflect");
		// reflect
		//SetAlphaOut_Reflexion(n1, n2);
		Reflect(nv,cosIn);
		newTissueID = currentTissueID; 	
		muasNew.x = nMuasGCurr.y;
		muasNew.y = nMuasGCurr.z;
//		this->refTri = triID;
	//printf("refTri: %i\n",this->refTri);
//printf("\t\t%s %i\n","newTissueID = currentTissueID:", currentTissueID);
	}
	
	
//printf("Current Tid: %i : mua: %f , mus: %f \n", currentTissueID, muasNew.x,muasNew.y);
	//if(currentTissueID != newTissueID && (mua[newTissueID] == 0 || mus[newTissueID] == 0)){
	if(muasNew.x == 0 || muasNew.y == 0){
		//int gid = get_global_id(0);
	//printf("\tID: %i | Weight: %e |" , get_global_id(0), GetWeight());
		this->SFTxy.y = 0;
		float dist = GetDistanceToSurface(this->p, POI, eZ);
	//printf(" Dist: %e \n", dist);
		if(dist>distMin){
			this->SFTxy.z = this->w*(1-Rsp);		
			
		}else{
		
			if(this->traversedFetalTissue){
				this->SFTxy.y = this->w*(1-Rsp);
			}
			this->SFTxy.x = this->w*(1-Rsp);			
		
		}
		
		do{
			SetWeight(0);
		}while(this->w>0);
		
	//	printf("Weight after: %f \n", this->w);
			
		
	}

//printf("\t\t v_new: (%f, %f, %f)\n",this->v.x,this->v.y,this->v.z);			
	return newTissueID;
}

bool Roulette(float chance, float random){
		
	if(random < chance){
		
		
		this->w /= chance;
		return true;
		
	}else {
		
		this->w = 0;
		return false;
		
	}
	
}

int FindTriangleInPath(float3 pPhoton, float3 vPhoton, int numTriangles,__global float* triangles){
	
// pPhoton: photon position
// vPhoton: photon trajectory
	
	this->steptoboundary = 10000;
	int triangleID = -1;
			
	for(int j=0; j<numTriangles; j++) {
				
				if(j!=this->refTri){
      		int index = j*3;
	      // Read the first point of Triangle KLM
					
		float3 K;				
				K = vload3(index,triangles);
	      // Read the second point of Triangle KLM
    float3 L;				
					L = vload3(index+1,triangles);
	      // Read the third point of Triangle KLM
		float3 M;			
					M = vload3(index+2,triangles);
//printf("Triangle %i: K(%f,%f,%f) , L(%f,%f,%f) , M(%f,%f,%f) , pPhoton(%f, %f, %f) \n",j,K.x,K.y,K.z,L.x,L.y,L.z,M.x,M.y,M.z, pPhoton.x, pPhoton.y,pPhoton.z);
	      // Compute vectors E, F, and G
	      float3 E = K - M;
	      float3 F = L - M;
	      float3 G = pPhoton - M;
	      // Compute the vector containing t and the coordinates k and l
	      float3 tkl;
				//float f = dot(cross(vPhoton, F), E);
			tkl.x	=  dot(cross(G, E), F);
			float w = dot(cross(vPhoton, F), E);
			
//			if(tkl.x/w<(this->steptoboundary)){
					
			tkl.y = dot(cross(vPhoton, F), G);
			tkl.z = dot(cross(G, E), vPhoton);
			
//			}
			tkl/=w;
 //printf("\tin FindTriangle: j = %i\t\tt = %f | k = %f | l = %f\n",j,tkl.x,tkl.y,tkl.z);
	      // Check if t and the intersection point (k, l) are acceptable
	      if(tkl.x < (this->steptoboundary) && tkl.x > 0.0f && tkl.y >= 0.0f && tkl.z >= 0.0f && (tkl.y + tkl.z) < 1) {
					this->steptoboundary = tkl.x;
	         triangleID = j;
	      }
			}
	}
		
	//if(tmin[0]<0.0001f){
		//tmin[0] = 0;
	//}
	
	//mem_fence(CLK_LOCAL_MEM_FENCE);
	return triangleID;
	
}
int FindTriangleInPath2(float3 pPhoton0,float stepSize, float3 vPhoton, int numTriangles,__global float* triangles, __global float* normals){
	
// pPhoton: photon position
// vPhoton: photon trajectory
	
	this->steptoboundary = 10000;
	int triangleID = -1;
	float3 	pPhoton1;
	pPhoton1 = pPhoton0+stepSize*vPhoton;

	float3 nv;
	
int rounds = 0;
int r2 = 0;
	for(int j=0; j<numTriangles; j++) {
				
		if(j!=this->refTri){
      		int index = j*3;
	      // Read the first point of Triangle KLM
					
		
		float3 M;			
		M = vload3(index+2,triangles);
		float3 G = pPhoton0 - M;
		nv = normalize(vload3(j,normals));
		float d0 = dot(nv,(-G));
		float d1 = dot(nv,(M-pPhoton1));
		
		if(fabs(d0)<=stepSize && (d0*d1)<0.0f){
rounds++;
		float3 K;				
		K = vload3(index,triangles);
	      // Read the second point of Triangle KLM
    		float3 L;				
		L = vload3(index+1,triangles);

		//if(dot((K-pPhoton0),vPhoton)>0 || dot((L-pPhoton0),vPhoton)>0 || dot((-G),vPhoton)>0){ // alles nach Ebene E0 in v-Richtung
		//	if(dot((K-pPhoton1),vPhoton)<0 || dot((L-pPhoton1),vPhoton)<0 || dot((M-pPhoton1),vPhoton)<0){ // alles vor Ebene E1 in v-Richtung
r2++;
//printf("Tri: %i  | d0: %f | d1: %f | stepSize: %f | rounds: %i  | r2: %i \n",j,d0, d1,stepSize, rounds,r2);

		
	      // Read the third point of Triangle KLM

	      		// Compute vectors E, F, and G
	      		float3 E = K - M;
	      		float3 F = L - M;
//	      		float3 G = pPhoton0 - M;
	      		// Compute the vector containing t and the coordinates k and l
	      		float3 tkl;
			//float f = dot(cross(vPhoton, F), E);
			
			
			float w = dot(cross(vPhoton, F), E);
			tkl.y = dot(cross(vPhoton, F), G)/w;
			
			if(tkl.y >= 0.0f && tkl.y < 1){
			
				tkl.z = dot(cross(G, E), vPhoton)/w;
			
				if(tkl.z>=0.0f &&  (tkl.y + tkl.z) < 1) {
			tkl.x	=  dot(cross(G, E), F)/w;
			if(tkl.x > 0.0f && tkl.x < (this->steptoboundary)){
//			if(tkl.x/w<(this->steptoboundary)){
					
			
				
//printf("\t\t in ft: j = %i \t\tt = %f | k = %f \n",j,tkl.x,tkl.y);
				
				
			
				
					
//			}
			//tkl/=w;
 //printf("\t in FindTriangle: j = %i\t\tt = %f | k = %f | l = %f\n",j,tkl.x,tkl.y,tkl.z);
	      // Check if t and the intersection point (k, l) are acceptable
	    // 		if(tkl.x < (this->steptoboundary) && tkl.x > 0.0f && tkl.y >= 0.0f && tkl.z >= 0.0f && (tkl.y + tkl.z) < 1) {
						
				
					this->steptoboundary = tkl.x;
	        		triangleID = j;
//printf("\t TriID: %i | stepToBoundary: %f  \n",triangleID,this->steptoboundary);
break;
	      		}
					}
				
			//}
				}
		}
		//}
		}
	}
		
	//if(tmin[0]<0.0001f){
		//tmin[0] = 0;
	//}
	
	//mem_fence(CLK_LOCAL_MEM_FENCE);
//printf("\t \t final TriID: %i with %i rounds\n",triangleID, rounds);
	return triangleID;
	
}

private:

	

};

float DoRand3(__local float* randArr, long MB, long MZ, __local int* inextInextp){
	int2 tempInextInextp = vload2(0,inextInextp);
	//int tempInextp = inextInextp[1];
	
	float random = 0;
		// Reshuffle Rand3-Array
	//	tempInext++;
	//	tempInextp++;
	tempInextInextp++;
	
	//atomic_inc(&inextInextp[0]);
	//atomic_inc(&inextInextp[1]);
	
	if(tempInextInextp.x==56){tempInextInextp.x=1;}
	if(tempInextInextp.y==56){tempInextInextp.y=1;}
	
	//atomic_cmpxchg(&inextInextp[0],56,1);
	//atomic_cmpxchg(&inextInextp[1],56,1);
	
	//	if (tempInext == 56) {tempInext=1;}
  	//if (tempInextp == 56) {tempInextp=1;}
		
		
		random=(randArr[tempInextInextp.x])-(randArr[tempInextInextp.y]);
	
		if (random < 0) {random += MB;}		
		
		atomic_xchg(&randArr[tempInextInextp.x],random);	
	
		vstore3((int3)(tempInextInextp,0),0,inextInextp);
//	atomic_xchg(&inextInextp[0] , tempInext);
//	atomic_xchg(&inextInextp[1] ,tempInextp);
	
		//atomic_xchg(&inextInextp[2],0); // unlock
	
	mem_fence(CLK_LOCAL_MEM_FENCE);
	
	return random/MB;
	
}

float GetRandomNumber(__local int* inextInextp, __local float* rand3Ar){
		
	float rand = -1;
	//	printf(" ID: %i locked: %i \n",get_global_id(0), inextInextp[2]);
	do{
		if(atomic_xchg(&inextInextp[2],1) ==0 ){
			
			rand = DoRand3(rand3Ar, 1000000000, 0, inextInextp);
			
		}		
		
	}while(rand<0.0f);	
//mem_fence(CLK_GLOBAL_MEM_FENCE);
	return rand;
}


__kernel void MCLightPropagation(const int numTriangles, const float factor,
 __global const float* vStart, __global const int* vStartOcc, const float3 pStart,
const float chance, const float w_threshold, int fetalTissueID,
__global float* SurfFetalTransRxyG, __global float* RecordAxzAyz,
__global const float* dimsXYZ, __global const float* e_P0XYZ , __global const float* gridDims,
__global float* const triangles, __global float* const normals,__global const int* mediumIDs,
 __global const float* nMuasG,
__global float* rand3Arr,__global int* inxtInxtpLocked){	    
	
	
	
	/* 
		* __global long* MBigMZero,
	 * __global int* inxtInxtp
	 * Gewebeeigenschaften: jeweils ein Vektor mit L?nge der Anzahl der Gewebe
	 * n: Tissue refraction indices: {n0,n1,n2,...}, Size: NumTissues; specific tissue: n(triN)
	 * mua: ?_a - Absorbtionskoeffizient
	 * mus: ?_s - Streuungskoeffizient
	 * g: anisotropie der Streuung
	 * 	
	 * triangles: Triangles: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3
	 * normals: Normals per Triangle: {(x0,y0,z0),(x1,y1,z1),...}, Size: NumTriangles*3
	 * mediumIDs: index of n and n_outer for every triangle {(n0_in,n0_out),(n1_in,n1_out),...}, Size: NumTriangles*2;
	 */
			
	int id = get_global_id(0);
	int si = get_global_size(0);
	int gid = get_group_id(0);
	int numg = get_num_groups(0);
	//const int gs = get_local_size(0);
		
	__local float rand3Ar[56];
	__local int inxtInxtpLock[3];
	//__local iSFT sftArr[64];
	//__local int count;
	
	__local float3 eX,eY,eZ,POI;
	__local float3 gridDimRPhiRho;
	__local float3 gridDimsXYZ;
	__local float3 DimsXYZ;
	//float3 aDimsXYZ;
	__local float Rsp;
event_t event[2];
event[0] = 0;
event[1] = 0;
event_t event2[2];
event2[0] = 0;
event2[1] = 0;
//event_t event3 = (event_t)0;
//printf("|ID: %i GS: %i ",id,si);

	if(id<si){
		
		Rsp	= Rspecular(nMuasG);
		event[0] = async_work_group_copy(rand3Ar,rand3Arr,56,event[0]);
		event[1] = async_work_group_copy(inxtInxtpLock,inxtInxtpLocked,3,event[1]);

		PhotonClass photon;
	
//		float stepsize;
//		float steptoboundary[1];
//		float unfinStepSize;
		int currTissueID = 0;
		
		//__local	int nextTissueID;
		float3 vStartRot = photon.GetvStart(vStart,vStartOcc,id);
	//	vStartRot = (float3)(0,0,-1);
		
		
		photon.Init(vStartRot, pStart, Rsp);	
//printf("InitWeight: %e | Rsp: %f\n", photon.GetWeight(), Rsp);
//printf("%f, %f, %f \n", photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z);	
		float3 nv;
		nv	= (float3)(0,0,-1);

		bool roulette = true;
		float randx = 0.0f;	
		float randy = 0.0f;
	//	int triangleID[1] = {0};
	
		__local float distMin;
		eX = vload3(1,e_P0XYZ);
		eY = vload3(2,e_P0XYZ);
		eZ = vload3(3,e_P0XYZ);
		POI = vload3(0,e_P0XYZ);
		gridDimRPhiRho = vload3(1,gridDims);
		gridDimsXYZ = vload3(0,gridDims);
		DimsXYZ = vload3(0,dimsXYZ);
		
	distMin	= gridDimsXYZ.z*dimsXYZ[2]/2.0f; // ==r = 11.8;
		//int ind = 0;
	//	int inter = 0;
		int stepNr = 0;
		//printf("distMin: %f | %f\n", distMin, gridDimsXYZ.z);
			
		
	//	float3 tempPos;
	//	float3 tempVec;
	//	tempVec = photon.GetHeadingVector();
	//	tempPos = photon.GetPhotonPosition();
	//	int temptr = 0;
	//	float firstStep = 0.0f;
		
	//	float4 randoms;
		float4 nmuasGCurr;
		float4 nmuasGNext;
		int2 innerOuterID;
	//printf("POI: %f, %f, %f\n", POI.x,POI.y,POI.z);
		//printf("Einheitsvektoren:\n \t e_x: %f , %f ,%f \n \t e_y: %f , %f , %f \n \t e_z: %f , %f , %f\n", eX.x,eX.y,eX.z,eY.x,eY.y,eY.z,eZ.x,eZ.y,eZ.z);
 wait_group_events(2, event);
		 
	barrier(CLK_LOCAL_MEM_FENCE);
//printf("FetalID: %i\n",fetalTissueID);
		do{
			
//printf("Randoms: %f, %f, %f\n", randx,randy,randz);
	// ID next triangle in path + steptoboundary				
			int triID = photon.FindTriangleInPath(photon.GetPhotonPosition(), photon.GetHeadingVector(), numTriangles, triangles);
							
		//	if(stepNr==0){firstStep = steptoboundary[0];}

			//int triID = triangleID[0];
//printf("ID: %i : triID: %i \n", id,triID);
	//	if(stepNr==0){
	//		photon.SetHeadingVector((float3)(0,0,-1));
			
			
	//	}
			
			//do{
		//	randoms.x = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
		//		randoms.x = 0.1;
				//}while(randoms.x<=0.0f);
	//		randoms.y = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
		//randoms.y = 0.62;
		//	randoms.z = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
		//	randoms.z = 0.33;
	//		randoms.w = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
		//	randoms.w = .84;
		
		//	printf("Randoms: %f, %f , %f , %f\n", randoms.x, randoms.y, randoms.z, randoms.w);
			
			
		 nmuasGCurr = vload4(currTissueID,nMuasG);
			
//printf("g(tID=%i): %f \n",currTissueID, nmuasGCurr.w);			
//printf("mua(%i): %f | mus(%i): %f \n", currTissueID, mua[currTissueID],currTissueID, mus[currTissueID]);
			if(nmuasGCurr.y == 0.0f && nmuasGCurr.z == 0.0f){
	//randx = GetRandomNumber(inxtInxtpPoolIDLock, randPool);
			//	photon.stepsize = 5;
				photon.stepsize = photon.steptoboundary;// -log(randx)*factor;
//printf("\t\t Ln. 1014 stepsize: %f ",photon.stepsize);

			}else{
//printf("\t\tLn. 1014: unfinStepSize: %f\n",photon.unfinStepSize);
				if(photon.unfinStepSize==0.0f){
//printf("\t\trandz: %f\n",randz);
					do{
					randx = GetRandomNumber(inxtInxtpLock,rand3Ar);
					}while(randx<=0.0f);
					//printf("rand: %f\n", randx);
	//photon.stepsize = 0.0005f;
					photon.stepsize = -log(randx)*factor/(nmuasGCurr.y+nmuasGCurr.z);
				//	if(photon.stepsize<0){
						//printf("%f ", photon.stepsize);
				//	}
					
					
				}else{
					photon.stepsize = photon.unfinStepSize/(nmuasGCurr.y+nmuasGCurr.z);		
					photon.unfinStepSize = 0.0f;

				}
				
				//if(photon.unfinStepSize<0.0f){
				//	photon.stepsize = -photon.unfinStepSize/(mua[currTissueID]+mus[currTissueID]);		
				//	photon.unfinStepSize = 0.0f;
				//}else
				

			}
			
		//	int triID = photon.FindTriangleInPath2(photon.GetPhotonPosition(),photon.stepsize, photon.GetHeadingVector(), numTriangles, triangles,normals);
		//	if(nmuasGCurr.y == 0.0f && nmuasGCurr.z == 0.0f){
		//		photon.stepsize = photon.steptoboundary;
		//	}

			// EnterTissue() : only 1st loop
			
		//	mem_fence(CLK_LOCAL_MEM_FENCE);
//if(steptoboundary[0]>0.0f && steptoboundary[0]<=stepsize) {
			// dot(nv,this->v)==0
//printf("ID : %i | Step: %i | ", id,stepNr);
//printf("\t triID: %i \n", triID);

//printf("\t\t Ln. 1334 stepsize: %e | STB: %e  | stpS-STB:  %e \n",photon.stepsize, photon.steptoboundary,photon.stepsize-photon.steptoboundary);
//float stpS = photon.stepsize-photon.steptoboundary;	
//printf("Position_old: %f, %f, %f  | w: %f  | tissueID: %i \n ", photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z,photon.GetWeight(),currTissueID);	
//printf(" v: %f, %f, %f  \n ", photon.GetHeadingVector().x,photon.GetHeadingVector().y,photon.GetHeadingVector().z);				
//printf(" nv: %f, %f, %f  \n ", nv.x,nv.y,nv.z);		
			if(photon.steptoboundary<=photon.stepsize) {
				
				
//printf("STB<stSz : %f \n",photon.steptoboundary-photon.stepsize);
		//	if(stpS>=0.00001f) {
	  			// set alpha_in and alpha_out
	  			photon.MoveStep(photon.steptoboundary); // changes position
//if(nmuasGCurr.y == 0.0f && nmuasGCurr.z == 0.0f){
//printf("\t%s : new Position: %f, %f, %f\n", "TransmitOrReflect", photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z);
//}
					randx = GetRandomNumber(inxtInxtpLock,rand3Ar);
				
					//if(triID>=0){
						innerOuterID = vload2(triID,mediumIDs);
						nv = normalize(vload3(triID,normals));
						
				
					//}
				
		//		if(stepNr==0){
		//		nv = (float3)(0,0,1);
		//		}
				
				if(currTissueID == innerOuterID.y){
					nmuasGNext = vload4(innerOuterID.x,nMuasG);
			//	printf("Next Tid: %i : n: %f , mua: %f , mus: %f , g: %f \n", innerOuterID.x, nmuasGNext.x,nmuasGNext.y,nmuasGNext.z,nmuasGNext.w);
				}else{
					nmuasGNext = vload4(innerOuterID.y,nMuasG);
			//		printf("Next Tid: %i : n: %f , mua: %f , mus: %f , g: %f \n", innerOuterID.y, nmuasGNext.x,nmuasGNext.y,nmuasGNext.z,nmuasGNext.w);
				}
				float3 cs = photon.GetPhotonPosition();
				if(fabs(sqrt(cs.x*cs.x+cs.y*cs.y)) < 0.01f){
				//if(stepNr==0){	
				//	printf(" nv_old: %f, %f, %f  \n ", nv.x,nv.y,nv.z);		
					nv = normalize(eZ)*sign(dot(nv,eZ));
				//	printf(" nv_new: %f, %f, %f  \n ", nv.x,nv.y,nv.z);		
				}

				
	  			currTissueID = photon.TransmitOrReflect(nmuasGCurr,nmuasGNext,POI,eZ, distMin, nv, currTissueID, innerOuterID, randx, triID, Rsp); // changes tissue or not and vector trajectory
				//	temptr = 1;
					photon.refTri = triID;
					if(currTissueID == fetalTissueID){
						photon.traversedFetalTissue = true;
					//	printf("fetalTissueID: %i", get_global_id(0));
					}
					
					photon.unfinStepSize = (photon.stepsize - photon.steptoboundary)*(nmuasGCurr.y+nmuasGCurr.z); 			
					nmuasGCurr= vload4(currTissueID,nMuasG);
//printf("%i : StepNr: %i | Weight: %e \n",id, stepNr,photon.GetWeight());

	  		}else{


				//randx = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
				//randy = GetRandomNumber(inxtInxtpPoolIDLock,MBigMZero,rand3Ar);
				photon.MoveStep(photon.stepsize); // changes position

				photon.Absorb(nmuasGCurr);// changes weight	

		
			randx = GetRandomNumber(inxtInxtpLock,rand3Ar); // rand for theta
			randy = GetRandomNumber(inxtInxtpLock,rand3Ar); // rand for psi


		//	printf("%f ", randy);
				//photon.DCScatter(nmuasGCurr.w, randx,randy);// changes vector trajectory (directional cosine
				photon.VScatter(nmuasGCurr.w, randx,randy); // Vector rotation
				//photon.QScatter(nmuasGCurr.w, randx,randy); // rotation via quaternions
			float3 coords = CalcCoords(photon.GetPhotonPosition(),POI,eX, eY, eZ);
			//RecordAbsorbtionXZYZ(RecordAxzAyz,photon.absorbedWeight, coords, gridDimsXYZ,DimsXYZ);
		//	RecordAbsorbtionRZ(RecordAxzAyz, photon.absorbedWeight, coords,nmuasGCurr.y);		
						
			RecordAbsorbtionRZ(RecordAxzAyz, photon.absorbedWeight, coords,(float2)(gridDimRPhiRho.x,gridDimsXYZ.z), DimsXYZ, nmuasGCurr.y);

		}
			
		if(photon.GetWeight()<w_threshold && photon.GetWeight() > 0.0f){
				randx = GetRandomNumber(inxtInxtpLock,rand3Ar);
				//printf("ID: %i | Below Threshold: weight: %f |",id,photon.GetWeight());
				roulette = photon.Roulette(chance,randx);
				//printf(" roulette: %i\n",roulette);
				
		}

	//printf("%f, %f, %f %f\n", photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z, photon.absorbedWeight);		
	//		inter = inter+3;
		
		//	printf("\t current P(t): %e , %e , %e\n",photon.GetPhotonPosition().x, photon.GetPhotonPosition().y, photon.GetPhotonPosition().z);
		//	if(fabs(photon.GetPhotonPosition().y)>500 || fabs(photon.GetPhotonPosition().x)>500 || photon.GetPhotonPosition().z>2.0f || photon.GetPhotonPosition().z< 0.0f){
				
		/*		
		printf("%i : StepNr: %i | T %i | StepS: %e | STB: %e | Weight: %e \n",id, stepNr,temptr, photon.stepsize, photon.steptoboundary[0],photon.GetWeight());
		printf("\t Current Tid: %i : n: %.2f , mua: %.2f , mus: %.2f , g: %.2f \n", currTissueID, nmuasGCurr.x,nmuasGCurr.y,nmuasGCurr.z,nmuasGCurr.w);
		printf("\t P(t-1): %e , %e , %e | P(t): %e , %e , %e\n",tempPos.x, tempPos.y, tempPos.z, photon.GetPhotonPosition().x, photon.GetPhotonPosition().y, photon.GetPhotonPosition().z);
		printf("\t V(t-1): %f , %f , %f | V(t): %f , %f , %f\n\n",tempVec.x,tempVec.y, tempVec.z, photon.GetHeadingVector().x,photon.GetHeadingVector().y,photon.GetHeadingVector().z);
		
		*/
				//photon.SetWeight(0);
				//mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		//	}
	/*
		tempVec = photon.GetHeadingVector();
			tempPos = photon.GetPhotonPosition();
			temptr = 0;
			*/
			//mem_fence(CLK_LOCAL_MEM_FENCE);
			
			// mem_fence(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
//	printf(",%f \n",photon.GetWeight());
// photon.GetWeight() > 0.0f && (
	
			stepNr++;		
			}while(photon.GetWeight() > 0.0f && (photon.GetWeight()>w_threshold || roulette));
		
		//printf("ID: %i | StepNr: %i | Position: %f , %f , %f | \n",id, stepNr-1, photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z);
		//printf("Surface | Fetal | Trans : %f , %f , %f \n", photon.GetSurfFetalTransxy().x, photon.GetSurfFetalTransxy().y,photon.GetSurfFetalTransxy().z);
		//RecordFSTxy(sftArr,lockRecord, photon.GetPhotonPosition(), POI, photon.GetSurfFetalTransxy(),eX, eY,eZ, gridDimRPhoRho);
		event2[0] = async_work_group_copy(rand3Arr,rand3Ar,56,event2[0]);
		event2[1] = async_work_group_copy(inxtInxtpLocked,inxtInxtpLock,3,event2[1]);
			
		//RecordFSTxy(sftArr,count,photon.GetPhotonPosition(), POI, photon.GetSurfFetalTransxy(),eX, eY,eZ, gridDimRPhoRho);
		//int idx = CalcIdx(photon.GetPhotonPosition(), POI,eX, eY, eZ, gridDimRPhiRho);

		float alpha = CalcAlphaOut(photon.GetHeadingVector(),nv);
		float radius = CalcRadius(photon.GetPhotonPosition(),POI,eZ);
		float3 coords = CalcCoords(photon.GetPhotonPosition(),POI,eX, eY, eZ);
	//	RecordFSTxy(SurfFetalTransRxyG,photon.GetSurfFetalTransxy(), idx);
//if(coords.z==0.0f){
//printf("%f  %f \n",radius,photon.GetSurfFetalTransxy().x);	
//printf("Position: %f , %f , %f | Coords: %f, %f, %f  | Radius: %fcm | w: %f \n",photon.GetPhotonPosition().x,photon.GetPhotonPosition().y,photon.GetPhotonPosition().z, coords.x,coords.y,coords.z, radius/factor, photon.GetWeight());	
//}
		RecordFSTxy(SurfFetalTransRxyG, photon.GetSurfFetalTransxy(), coords, alpha, radius);

//	RecordFSTxy(SurfFetalTransRxyG, photon.GetSurfFetalTransxy(), idx);
	wait_group_events(2, event2);
	
			
	}
	
	
barrier(CLK_GLOBAL_MEM_FENCE);
//printf("%s\n","EOF"); 

}

