/* Kernel3.cl - Created with AMD CodeXL */
__kernel void addition(const int num, __global int* A, __global int* B, __global int* C)
{
	
	int id = get_global_id(0);
	
	
	if(id<num){
		
		//if(id==1){
		//	printf("%i : %i \t",id,C[id]);

		//}


		C[id] = A[id] + B[id] + C[id];	
	
	//	if(id==1){
	//		printf("%i \n",C[id]);

	//	}
	}
	
	
	
	
}
