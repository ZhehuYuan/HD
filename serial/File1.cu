#include <iostream>
#include<stdio.h>
// For the CUDA runtime routines (prefixed with "cuda_")
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include "points.h"

#define TEST 1

__global__ void kernel() {
    printf("Hello from GPU\n");
}

__global__ void MyHD(Point* A, Point* B, unsigned int* cMax, int l1, int l2, int* Aloc, long* gpuCount, unsigned int* top) {
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	int subWorkerId = id % 32;
	unsigned int mask = 0xFFFFFFFF;

	unsigned int i = 0;

	while(true){
		if(subWorkerId == 0)i = atomicAdd(top, 1);
		i = __shfl_sync(0xFFFFFFFF, i, 0);
		if (i >= l1)break;
		Point a = A[i];
		unsigned int cmin = INT_MAX;
		int preindex = Aloc[i];
		for (int j = subWorkerId; j < l2; j += 32) {
			if (0 <= preindex - j) {
				int x = a.l - B[preindex - j].l;
				int y = a.w - B[preindex - j].w;
				int z = a.h - B[preindex - j].h;
				unsigned int disleft = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (disleft < *cMax || cmin < *cMax) {
					cmin = 0;
				}
				else if (disleft < cmin) {
					cmin = disleft;
				}
			}
			if (preindex + j < l2) {
				int x = a.l - B[preindex + j].l;
				int y = a.w - B[preindex + j].w;
				int z = a.h - B[preindex + j].h;
				unsigned int disright = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (disright < *cMax || cmin < *cMax) {
					cmin = 0;
				}
				else if (disright < cmin) {
					cmin = disright;
				}
			}
			cmin = __reduce_min_sync(mask, cmin);
			if (cmin == 0)break;
		}
		if (subWorkerId == 0) {
			if (*cMax < cmin) {
				atomicMax(cMax, cmin);
			}
		}
	}
}

__global__ void ZHD(Point* A, Point* B, unsigned int* cMax, int l1, int l2, long* gpuCount) {
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	const int workerId = id / 32;
	const int subWorkerId = id - workerId * 32;
	const int nWorker = blockDim.x * gridDim.x / 32;
	unsigned int mask = 0xFFFFFFFF;

	int jobPerWorker = l1 / nWorker;
	int extraJob = l1 - jobPerWorker * workerId;

	int preindex = 0;
	
	if (workerId < extraJob) {
		for (int i = workerId * (jobPerWorker + 1); i < (workerId + 1) * (jobPerWorker + 1); i++) {
			Point a = A[i];
			unsigned int cmin = INT_MAX;
			int minplace = 0;
			int round = 0;
			for (int j = subWorkerId; j < l2; j += 32) {
				bool modified = false;
				if (0 <= preindex - j) {
					int x = a.l - B[preindex - j].l;
					int y = a.w - B[preindex - j].w;
					int z = a.h - B[preindex - j].h;
					unsigned int disleft = x * x + y * y + z * z;
					if (TEST) {
						gpuCount[id]++;
					}
					if (cmin < *cMax) {
						cmin = 0;
					}
					else if (disleft < *cMax) {
						preindex = preindex - j;
						cmin = 0;
						modified = true;
					}
					else if (disleft < cmin) {
						cmin = disleft;
						minplace = preindex - j;
						modified = true;
					}
				}
				if (preindex + j < l2) {
					int x = a.l - B[preindex + j].l;
					int y = a.w - B[preindex + j].w;
					int z = a.h - B[preindex + j].h;
					unsigned int disright = x * x + y * y + z * z;
					if (TEST) {
						gpuCount[id]++;
					}
					if (cmin < *cMax) {
						cmin = 0;
					}
					else if (disright < *cMax) {
						preindex = preindex + j;
						modified = true;
						cmin = 0;
					}
					else if (disright < cmin) {
						cmin = disright;
						minplace = preindex + j;
						modified = true;
					}
				}
				int tmp = __reduce_min_sync(mask, cmin);
				if (modified) {
					if (tmp == cmin) {
						round = j / 32;
					}
					modified = false;
				}
				cmin = tmp;
				if (cmin == 0)break;
			}
			if (subWorkerId == 0) {
				if (*cMax < cmin) {
					atomicMax(cMax, cmin);
				}
			}
			int tmp = __reduce_max_sync(mask, round);
			if (tmp != round)minplace = 0;
			preindex = __reduce_max_sync(mask, minplace);
		}
	}
	else {
		for (int i = workerId * jobPerWorker + extraJob; i < (workerId + 1) * jobPerWorker + extraJob; i++) {
			Point a = A[i];
			unsigned int cmin = INT_MAX;
			int minplace = 0;
			int round = 0;
			for (int j = subWorkerId; j < l2; j+=32) {
				bool modified = false;
				if (0 <= preindex - j) {
					int x = a.l - B[preindex - j].l;
					int y = a.w - B[preindex - j].w;
					int z = a.h - B[preindex - j].h;
					unsigned int disleft = x * x + y * y + z * z;
					if (TEST) {
						gpuCount[id]++;
					}
					if (cmin < *cMax) {
						cmin = 0;
					}
					else if (disleft < *cMax) {
						preindex = preindex - j;
						cmin = 0;
						modified = true;
					}
					else if (disleft < cmin) {
						cmin = disleft;
						minplace = preindex - j;
						modified = true;
					}
				}
				if (preindex + j < l2) {
					int x = a.l - B[preindex + j].l;
					int y = a.w - B[preindex + j].w;
					int z = a.h - B[preindex + j].h;
					unsigned int disright = x * x + y * y + z * z;
					if (TEST) {
						gpuCount[id]++;
					}
					if (cmin < *cMax) {
						cmin = 0;
					}
					else if (disright < *cMax) {
						preindex = preindex + j;
						modified = true;
						cmin = 0;
					}
					else if (disright < cmin) {
						cmin = disright;
						minplace = preindex + j;
						modified = true;
					}
				}
				int tmp = __reduce_min_sync(mask, cmin);
				if (modified) {
					if (tmp == cmin) {
						round = j / 32;
					}
					modified = false;
				}
				cmin = tmp;
				if (cmin == 0)break;
			}
			if (subWorkerId == 0) {
				if (*cMax < cmin) {
					atomicMax(cMax, cmin);
				}
			}
			int tmp = __reduce_max_sync(mask, round);
			if (tmp != round)minplace = 0;
			preindex = __reduce_max_sync(mask, minplace);
		}
	}
	
}

__global__ void ZHD2(Point* A, Point* B, unsigned int* cMax, int l1, int l2, long* gpuCount, unsigned int* top) {
	const int id = blockDim.x * blockIdx.x + threadIdx.x;
	const int subWorkerId = id % 32;
	unsigned int mask = 0xFFFFFFFF;

	int preindex = 0;
	unsigned int i = 0;

	while (true) {
		if (subWorkerId == 0)i = atomicAdd(top, 1);
		i = __shfl_sync(0xFFFFFFFF, i, 0);
		if (i >= l1)break;
		Point a = A[i];
		unsigned int cmin = INT_MAX;
		int minplace = 0;
		int round = 0;
		for (int j = subWorkerId; j < l2; j += 32) {
			bool modified = false;
			if (0 <= preindex - j) {
				int x = a.l - B[preindex - j].l;
				int y = a.w - B[preindex - j].w;
				int z = a.h - B[preindex - j].h;
				unsigned int disleft = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (cmin < *cMax) {
					cmin = 0;
				}
				else if (disleft < *cMax) {
					preindex = preindex - j;
					cmin = 0;
					modified = true;
				}
				else if (disleft < cmin) {
					cmin = disleft;
					minplace = preindex - j;
					modified = true;
				}
			}
			if (preindex + j < l2) {
				int x = a.l - B[preindex + j].l;
				int y = a.w - B[preindex + j].w;
				int z = a.h - B[preindex + j].h;
				unsigned int disright = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (cmin < *cMax) {
					cmin = 0;
				}
				else if (disright < *cMax) {
					preindex = preindex + j;
					modified = true;
					cmin = 0;
				}
				else if (disright < cmin) {
					cmin = disright;
					minplace = preindex + j;
					modified = true;
				}
			}
			int tmp = __reduce_min_sync(mask, cmin);
			if (modified) {
				if (tmp == cmin) {
					round = j / 32;
				}
				modified = false;
			}
			cmin = tmp;
			if (cmin == 0)break;
		}
		if (subWorkerId == 0) {
			if (*cMax < cmin) {
				atomicMax(cMax, cmin);
			}
		}
		int tmp = __reduce_max_sync(mask, round);
		if (tmp != round)minplace = 0;
		preindex = __reduce_max_sync(mask, minplace);
	}
}

__global__ void EBHD(Point* A, Point* B, unsigned int* cMax, int l1, int l2, long* gpuCount, unsigned int* top) {
	int id = blockDim.x * blockIdx.x + threadIdx.x;
	int subWorkerId = id % 32;
	unsigned int mask = 0xFFFFFFFF;

	unsigned int i = 0;

	while (true) {
		if (subWorkerId == 0)i = atomicAdd(top, 1);
		i = __shfl_sync(0xFFFFFFFF, i, 0);
		if (i >= l1)break;
		Point a = A[i];
		unsigned int cmin = INT_MAX;
		int preindex = 0;
		for (int j = subWorkerId; j < l2; j += 32) {
			if (0 <= preindex - j) {
				int x = a.l - B[preindex - j].l;
				int y = a.w - B[preindex - j].w;
				int z = a.h - B[preindex - j].h;
				unsigned int disleft = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (disleft < *cMax || cmin < *cMax) {
					cmin = 0;
				}
				else if (disleft < cmin) {
					cmin = disleft;
				}
			}
			if (preindex + j < l2) {
				int x = a.l - B[preindex + j].l;
				int y = a.w - B[preindex + j].w;
				int z = a.h - B[preindex + j].h;
				unsigned int disright = x * x + y * y + z * z;
				if (TEST) {
					gpuCount[id]++;
				}
				if (disright < *cMax || cmin < *cMax) {
					cmin = 0;
				}
				else if (disright < cmin) {
					cmin = disright;
				}
			}
			cmin = __reduce_min_sync(mask, cmin);
			if (cmin == 0)break;
		}
		if (subWorkerId == 0) {
			if (*cMax < cmin) {
				atomicMax(cMax, cmin);
			}
		}
	}
}


/*
extern "C" void gpuAlloc(int l1, int l2) {
	int size1 = l1 * sizeof(Point);
	int size2 = l2 * sizeof(Point);
	int size3 = l1 * sizeof(int);

	//getting details of your CUDA device
	cudaDeviceProp props;
	cudaGetDeviceProperties(&props, 0); //change to the proper index of your cuda device
	const int threadsPerBlock = props.maxThreadsPerBlock;
	const int blocksPerGrid = props.multiProcessorCount;

	//allocating the input data in the GPU
	cudaMalloc(&gpuA, size1);
	cudaMalloc(&gpuB, size2);

	cudaMalloc(&cMax, sizeof(unsigned int));

	cudaMalloc(&gpuAloc, size3);

	cudaMalloc(&gpuCount, threadsPerBlock * blocksPerGrid * sizeof(long));
	
	if (TEST)printf("%d, %d\n", threadsPerBlock, blocksPerGrid);
}
*/
extern "C" double HDparallel(Point* A, Point* B, int l1, int l2, int method, int* Aloc) {
	// method 0 = EBHD
	// method 1 = ZHD
	// method 2 = MyHD
	// method 3 = ZHD2

	cudaError_t error;
	int size1 = l1 * sizeof(Point);
	int size2 = l2 * sizeof(Point);
	int size3 = l1 * sizeof(int);

	//getting details of your CUDA device
	cudaDeviceProp props;
	cudaGetDeviceProperties(&props, 0); //change to the proper index of your cuda device
	const int threadsPerBlock = props.maxThreadsPerBlock;
	const int blocksPerGrid = props.multiProcessorCount;


	auto start = std::chrono::high_resolution_clock::now();
	Point* gpuA, * gpuB;
	unsigned int* cMax;
	unsigned int* top;
	unsigned int t = 0;
	unsigned int t2 = 0;
	int* gpuAloc;
	long* gpuCount;
	long* tmpCount;

	cudaMalloc(&gpuA, size1);
	cudaMalloc(&gpuB, size2);

	cudaMalloc(&cMax, sizeof(unsigned int));
	cudaMalloc(&top, sizeof(unsigned int));

	cudaMalloc(&gpuAloc, size3);

	cudaMalloc(&gpuCount, threadsPerBlock * blocksPerGrid * sizeof(long));

	if (TEST)printf("%d, %d\n", threadsPerBlock, blocksPerGrid);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Memory allocation time: " << duration.count() << " ms" << std::endl;

	start = std::chrono::high_resolution_clock::now();
	cudaMemcpy(gpuA, A, size1, cudaMemcpyHostToDevice);
	cudaMemcpy(gpuB, B, size2, cudaMemcpyHostToDevice);
	cudaMemcpy(cMax, &t, sizeof(unsigned int), cudaMemcpyHostToDevice);
	cudaMemcpy(top, &t2, sizeof(unsigned int), cudaMemcpyHostToDevice);
	if (method == 2 || method == 4) {
		cudaMemcpy(gpuAloc, Aloc, size3, cudaMemcpyHostToDevice);
	}
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::cout << "CUDA error 1: " << cudaGetErrorString(error) << std::endl;
	}
	if (TEST) {
		tmpCount = new long[threadsPerBlock * blocksPerGrid](); for (int i = 0; i < threadsPerBlock * blocksPerGrid; i++) {
			tmpCount[i] = 0;
		}
		cudaMemcpy(gpuCount, tmpCount, threadsPerBlock * blocksPerGrid * sizeof(long), cudaMemcpyHostToDevice);
	}
	end = std::chrono::high_resolution_clock::now();
	duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	std::cout << "Memory copy time: " << duration.count() << " ms" << std::endl;

	//computing h(A,B)
	if (method == 0) {
		EBHD << < blocksPerGrid, threadsPerBlock >> > (gpuA, gpuB, cMax, l1, l2, gpuCount, top);
	}else if (method == 1) {
		ZHD << < blocksPerGrid, threadsPerBlock >> > (gpuA, gpuB, cMax, l1, l2, gpuCount);
	}else if (method == 2) {
		MyHD << < blocksPerGrid, threadsPerBlock >> > (gpuA, gpuB, cMax, l1, l2, gpuAloc, gpuCount, top);
	}
	else if (method == 3) {
		ZHD2 << < blocksPerGrid, threadsPerBlock >> > (gpuA, gpuB, cMax, l1, l2, gpuCount, top);
	}

	cudaDeviceSynchronize();
	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::cout << "CUDA error 2: " << cudaGetErrorString(error) << std::endl;
	}

	//copying the result back
	unsigned int distance;
	cudaMemcpy(&distance, cMax, sizeof(unsigned int), cudaMemcpyDeviceToHost);
	if (TEST) {
		cudaMemcpy(tmpCount, gpuCount, threadsPerBlock * blocksPerGrid * sizeof(long), cudaMemcpyDeviceToHost);
		long long tmp = 0;
		long minElem = tmpCount[0];
		long maxElem = tmpCount[0];
		for (int i = 0; i < threadsPerBlock * blocksPerGrid; i++) {
			tmp += tmpCount[i];
			if (tmpCount[i] > maxElem)maxElem = tmpCount[i];
			if (tmpCount[i] < minElem)minElem = tmpCount[i];
		}
		printf("Count: %lld, (%ld, %ld)\n", tmp, maxElem, minElem);
	}
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::cout << "CUDA error 3: " << cudaGetErrorString(error) << std::endl;
	}

	//freeing memory
	cudaFree(gpuA); 
	cudaFree(gpuB);
	cudaFree(cMax);
	cudaFree(gpuAloc);
	cudaFree(gpuCount);

	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		std::cout << "CUDA error 4: " << cudaGetErrorString(error) << std::endl;
	}

	return sqrt((double)distance);
}