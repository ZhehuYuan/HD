// serial.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <set>
#include <limits>
#include <vector>
#include <random>
#include <chrono>
#include <cuda_runtime.h>

#include "points.h"
#include "wkt_loader.h"

#define UINT_MAX 4294967295
#define INT_MAX 2147483647
#define TEST 1

extern "C" double HDparallel(Point * A, Point * B, int l1, int l2, int method, int* Aloc);
//extern "C" void gpuAlloc(int l1, int l2);

using namespace std;


class EBHD
{
public:
    int l1, l2;
    Point* A, * B, * E;
    unsigned int count;

    EBHD(int _l1, int _l2, Point* _A, Point* _B){
        l1 = _l1;
        l2 = _l2;
        B = _B;
        A = _A;
        E = new Point[l1];
        count = 0;
    }

    ~EBHD(){
        delete[] A;
        delete[] B;
        delete[] E;
    }

    void Randomize() {
        for (int i = 0; i < l1; i++) {
            E[i].swap(E[rand() % l1]);
        }
        for (int i = 0; i < l2; i++) {
            B[i].swap(B[rand() % l2]);
        }
    }

    void Excluding() {
        set<Point> setA(A, A + l1);
        set<Point> setB(B, B + l2);
        set<Point> result;
        set_difference(setA.begin(), setA.end(), setB.begin(), setB.end(), inserter(result, result.end()));
        l1 = result.size();
        copy(result.begin(), result.end(), E);
    }

    double EarlyBreakDirectedHD() {
        double cmax = 0;
        Excluding();
        Randomize();
        for (int i = 0; i < l1; i++) {
            double cmin = numeric_limits<double>::max();
            bool breaked = false;
            for (int j = 0; j < l2; j++) {
                double d = E[i].Distance(B[j]);
                if (TEST)count++;
                if (d < cmax) {
                    breaked = true;
                    break;
                }
                if (d < cmin)cmin = d;
            }
            if (cmin > cmax && !breaked)cmax = cmin;
        }

        return cmax;
    }
};

class NBHD
{
public:
    int l1, l2;
    Point* A, * B, * E;
    unsigned int count;

    NBHD(int _l1, int _l2, Point* _A, Point* _B) {
        l1 = _l1;
        l2 = _l2;
        B = _B;
        A = _A;
        E = new Point[l1];
        count = 0;
    }

    ~NBHD() {
        delete[] A;
        delete[] B;
        delete[] E;
    }

    void Randomize() {
        for (int i = 0; i < l1; i++) {
            E[i].swap(E[rand() % l1]);
        }
        for (int i = 0; i < l2; i++) {
            B[i].swap(B[rand() % l2]);
        }
    }

    void Excluding() {
        set<Point> setA(A, A + l1);
        set<Point> setB(B, B + l2);
        set<Point> result;
        set_difference(setA.begin(), setA.end(), setB.begin(), setB.end(), inserter(result, result.end()));
        l1 = result.size();
        copy(result.begin(), result.end(), E);
    }

    double EarlyBreakDirectedHD() {
        double cmax = 0;
        Excluding();
        Randomize();
        for (int i = 0; i < l1; i++) {
            double cmin = numeric_limits<double>::max();
            bool breaked = false;
            for (int j = 0; j < l2; j++) {
                double d = E[i].Distance(B[j]);
                if (TEST)count++;
                if (d < cmax) {
                    breaked = true;
                    break;
                }
                if (d < cmin)cmin = d;
            }
            if (cmin > cmax && !breaked)cmax = cmin;
        }

        return cmax;
    }
};

class ZHD
{
public:
    int l1, l2;
    Point* A, * B, * Az, * Bz;
    PointCMin* Z1, * Z2;
    unsigned int count;

    ZHD(int _l1, int _l2, Point* _A, Point* _B) {
        l1 = _l1;
        l2 = _l2;
        B = _B;
        A = _A;
        Z1 = new PointCMin[l1];
        Z2 = new PointCMin[l2];
        Az = new Point[l1];
        Bz = new Point[l2];
        count = 0;

        ZOrder();
    }

    ZHD() {
        delete[] A;
        delete[] B;
        delete[] Az;
        delete[] Bz;
    }

    int interleaveBits(uint32_t x) {
        uint64_t result = x;
        result = (result | (result << 16)) & 0x0000FFFF0000FFFF;
        result = (result | (result << 8)) & 0x00FF00FF00FF00FF;
        result = (result | (result << 4)) & 0x0F0F0F0F0F0F0F0F;
        result = (result | (result << 2)) & 0x3333333333333333;
        result = (result | (result << 1)) & 0x5555555555555555;
        return (int)result;
    }

    int morton3D(Point p) {
        return (interleaveBits(p.l) << 2) | (interleaveBits(p.w) << 1) | interleaveBits(p.h);
    }

    void ZOrder() {
        for (int i = 0; i < l1; i++) {
            Z1[i] = { A[i], (double)morton3D(A[i]) };
        }
        for (int i = 0; i < l2; i++) {
            Z2[i] = { B[i], (double)morton3D(B[i]) };
        }
        sort(Z1, Z1 + l1);
        sort(Z2, Z2 + l2);
        for (int i = 0; i < l1; i++) {
            Az[i] = Z1[i].p;
        }
        for (int i = 0; i < l2; i++) {
            Bz[i] = Z2[i].p;
        }
    }

    double EarlyBreakDirectedHD() {
        double cmax = 0;
        int preindex = 0;
        // Az and Bz GPU version is on going
        for (int i = 0; i < l1; i++) {
            double cmin = numeric_limits<double>::max();
            int minplace = 0;
            for (int j = 0; j < l2; j++) {
                if (0 <= preindex - j) {
                    double disleft = Az[i].Distance(Bz[preindex - j]);
                    if (TEST)count++;
                    if (disleft < cmax) {
                        preindex = preindex - j;
                        cmin = 0;
                        break;
                    }
                    else if (disleft < cmin) {
                        cmin = disleft;
                        minplace = preindex - j;
                    }
                }
                if (preindex + j < l2) {
                    double disright = Az[i].Distance(Bz[preindex + j]);
                    if (TEST)count++;
                    if (disright < cmax) {
                        preindex = preindex + j;
                        cmin = 0;
                        break;
                    }
                    else if (disright < cmin) {
                        cmin = disright;
                        minplace = preindex + j;
                    }
                }
            }
            if (cmax < cmin) {
                cmax = cmin;
                preindex = minplace;
            }
        }
        return cmax;
    }
};

class MyHD
{
public:
    int x, y, z;
    int l1, l2;
    Point* B;
    vector<PointCMin> A;
    int* lb;
    Point** Bs;

    MyHD(int _l1, int _l2, PointCMin* _A, Point* _B, int _x, int _y, int _z){
        l1 = _l1;
        l2 = _l2;
        A.assign(_A, _A + l1);
        B = _B;
        x = _x;
        y = _y;
        z = _z;

        Bs = new Point * [x + y + z + 1];
        lb = new int[x + y + z + 1]();

        for (int i = 0; i <= x + y + z; i++) {
            Bs[i] = new Point[l2]();
        }
    }

    ~MyHD(){
        for (int i = 0; i <= x+y+z; i++){
            delete[] Bs[i];
        }
        delete[] Bs;
    }

    void CountSort() {
        for (int i = 0; i < l2; i++) {
            int tmp = B[i].l + B[i].w + B[i].h;
            Bs[tmp][lb[tmp]] = B[i];
            lb[tmp]++;
        }
    }

    void RandomizeB() {
        for (int i = 0; i < l2; i++) {
            B[i].swap(B[rand() % l2]);
        }
    } 
    
    void RandomizeA() {
        random_device rd;
        mt19937 g(rd());
        shuffle(A.begin(), A.end(), g);
    }

    double EarlyBreakDirectedHD() {
        // Maybe?
        //RandomizeA();
        //RandomizeB();
        CountSort();
        double cmax = 0;

        for (int i = 0; i < x + y + z + 1; i++) {
            vector<PointCMin>::iterator it = A.begin();
            while(it != A.end()){
                if (i >= sqrt(3) * (*it).cmin) {
                    if (cmax < (*it).cmin) {
                        cmax = (*it).cmin;
                    }
                    if (it != A.end() - 1) {
                        swap(*it, A.back());
                        A.pop_back();
                    }
                    else {
                        A.pop_back();
                        it = A.end();
                    }
                    continue;
                }

                int tmp = (*it).getValue();
                bool breaked = false;
                if (tmp - i >= 0) {
                    for (int j = 0; j < lb[tmp - i]; j++) {
                        double d = (*it).p.Distance(Bs[tmp - i][j]);
                        if (d <= cmax) {
                            if (it != A.end() - 1) {
                                swap(*it, A.back());
                                A.pop_back();
                            }
                            else {
                                A.pop_back();
                                it = A.end();
                            }
                            breaked = true;
                            break;
                        }
                        else if (d < (*it).cmin) {
                            (*it).cmin = d;
                        }
                    }
                }
                if (tmp + i < x + y + z + 1 && !breaked && i != 0) {
                    for (int j = 0; j < lb[tmp + i]; j++) {
                        double d = (*it).p.Distance(Bs[tmp + i][j]);
                        if (d <= cmax) {
                            if (it != A.end() - 1) {
                                swap(*it, A.back());
                                A.pop_back();
                            }
                            else {
                                A.pop_back();
                                it = A.end();
                            }
                            breaked = true;
                            break;
                        }
                        else if (d < (*it).cmin) {
                            (*it).cmin = d;
                        }
                    }
                }
                if(!breaked)++it;
            }
        }

        for (vector<PointCMin>::iterator it = A.begin(); it != A.end(); ++it) {
            if (cmax < (*it).cmin) {
                cmax = (*it).cmin;
            }
        }

        return cmax;
    }
};

class MyHD2 {
public:
    int l1, l2;
    Point* A, * B;
    MyPoint* Az, * Bz;
    int* Aloc;
    unsigned int count;

    MyHD2(int _l1, int _l2, Point* _A, Point* _B) {
        A = _A;
        B = _B;
        l1 = _l1;
        l2 = _l2;
        Az = new MyPoint[l1];
        Bz = new MyPoint[l2];
        Aloc = new int[l1];
        count = 0;

        MySort();
    }

    int interleaveBits(uint32_t x) {
        uint64_t result = x;
        result = (result | (result << 16)) & 0x0000FFFF0000FFFF;
        result = (result | (result << 8)) & 0x00FF00FF00FF00FF;
        result = (result | (result << 4)) & 0x0F0F0F0F0F0F0F0F;
        result = (result | (result << 2)) & 0x3333333333333333;
        result = (result | (result << 1)) & 0x5555555555555555;
        return (int)result;
    }

    int morton3D(Point p) {
        return (interleaveBits(p.l) << 2) | (interleaveBits(p.w) << 1) | interleaveBits(p.h);
    }

    void MySort() {
        for (int i = 0; i < l1; i++) {
            Az[i] = { A[i], morton3D(A[i])};
        }
        for (int i = 0; i < l2; i++) {
            Bz[i] = { B[i], morton3D(B[i])};
        }
        sort(Az, Az + l1);
        sort(Bz, Bz + l2);

        for (int i = 0; i < l1; i++) {
            A[i] = Az[i].p;
        }
        for (int i = 0; i < l2; i++) {
            B[i] = Bz[i].p;
        }
    }

    void AmapB() {
        int j = 0;
        for (int i = 0; i < l2; i++) {
            if (Bz[i].mortonCode >= Az[j].mortonCode) {
                Aloc[j] = i;
                j++;
                i--;
                if (j >= l1)break;
            }
        }
        for (; j < l1; j++) {
            Aloc[j] = l2 - 1;
        }
    }

    double EarlyBreakDirectedHD() {
        double cmax = 0;
        for (int j = 0; j < l1; j++) {
            int loc = Aloc[j];
            double cmin = numeric_limits<double>::max();
            Point tmp = Az[j].p;
            for (int i = 0; i < l2; i++) {
                if (loc + i < l2) {
                    double d = tmp.Distance(Bz[loc + i].p);
                    if (TEST)count++;
                    if (d <= cmax) {
                        cmin = 0;
                        break;
                    }
                    else if(d < cmin){
                        cmin = d;
                    }
                }
                if (loc >= i) {
                    double d = tmp.Distance(Bz[loc - i].p);
                    if(TEST)count++;
                    if (d <= cmax) {
                        cmin = 0;
                        break;
                    }
                    else if (d < cmin) {
                        cmin = d;
                    }
                }
                if (loc < i && loc + i >= l2)break;
            }
            if (cmin > cmax) cmax = cmin;
        }
        return cmax;
    }
};

int main()
{
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);

    //std::cout << "Device name: " << prop.name << std::endl;
    //std::cout << "Max shared memory per block: " << prop.sharedMemPerBlock << " bytes" << std::endl;

    srand(time(nullptr));

    Point* A, * B, * A2, * B2, * A3, * B3, * A4, * B4, * A5, * B5, * A6, * B6, * A7, * B7;
    int l1, l2;
    int* tmp = (int*) nullptr;

    l1 = 500000;
    l2 = 500000;

    A = LoadPoints("D:\\note\\OSU\\Project Hausdorff Distance\\code\\serial\\serial\\hd_datasets\\lakes.bz2.wkt", &l1);
    B = LoadPoints("D:\\note\\OSU\\Project Hausdorff Distance\\code\\serial\\serial\\hd_datasets\\parks_Europe.wkt", &l2);

    int maxl = numeric_limits<int>::min();
    int maxw = numeric_limits<int>::min();
    int minl = numeric_limits<int>::max();
    int minw = numeric_limits<int>::max();

    for (int i = 0; i < l1; i++) {
        Point p = A[i];
        if (p.l < minl) minl = p.l;
        if (p.w < minw) minw = p.w;
        if (p.l > maxl) maxl = p.l;
        if (p.w > maxw) maxw = p.w;
    }

    for (int i = 0; i < l2; i++) {
        Point p = B[i];
        if (p.l < minl) minl = p.l;
        if (p.w < minw) minw = p.w;
        if (p.l > maxl) maxl = p.l;
        if (p.w > maxw) maxw = p.w;
    }
    
    int x = maxl - minl;
    int ratio = 1;
    while (x / ratio > 1000) {
        ratio += 1;
    }
    x /= ratio;

    int y = maxw - minw;
    while (y / ratio > 1000) {
        ratio += 1;
    }
    y /= ratio;

    int z = 1;
    
    printf("%d, %d, %d, (%d, %d, %d, %d), %d, %d\n", x, y, z, maxl, maxw, minl, minw, ratio, ratio);

    std::ofstream fileA("pointsA.csv");
    std::ofstream fileB("pointsB.csv");

    for (int i = 0; i < l1; i++) {
        A[i].l -= minl;
        A[i].l /= ratio;
        A[i].w -= minw;
        A[i].w /= ratio;
        //fileA << A[i].l << "," << A[i].w << "\n";
    }
    for (int i = 0; i < l2; i++) {
        B[i].l -= minl;
        B[i].l /= ratio;
        B[i].w -= minw;
        B[i].w /= ratio;
        //fileB << B[i].l << "," << B[i].w << "\n";
    }
    fileA.close();
    fileB.close();

    A2 = new Point[l1];
    B2 = new Point[l2];

    A3 = new Point[l1];
    B3 = new Point[l2];

    A4 = new Point[l1];
    B4 = new Point[l2];

    A5 = new Point[l1];
    B5 = new Point[l2];

    A6 = new Point[l1];
    B6 = new Point[l2];

    A7 = new Point[l1];
    B7 = new Point[l2];

    /*int x = 1024;
    int y = 1024;
    int z = 1;

    int modify = 50;

    A = new Point[l1];
    B = new Point[l2];
    */
    for (int i = 0; i < l1; i++) {
        //A[i] = { rand() % (x - modify), rand() % (y - modify), rand() % (z - modify) };
        A2[i] = A[i];
        A3[i] = A[i];
        A4[i] = A[i]; 
        A5[i] = A[i];
        A6[i] = A[i];
        A7[i] = A[i];
    }
    for (int i = 0; i < l2; i++) {
        //B[i] = { rand() % (x - modify) + modify, rand() % (y - modify) + modify, rand() % (z - modify) + modify };
        B2[i] = B[i];
        B3[i] = B[i];
        B4[i] = B[i];
        B5[i] = B[i];
        B6[i] = B[i];
        B7[i] = B[i];
    }

    double cmin = numeric_limits<double>::max();


    ZHD* HD = new ZHD(l1, l2, A, B);
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << HD->EarlyBreakDirectedHD() << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "ZHD Execution time: " << duration.count() << " ms // " << HD->count << "count" << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    /*
    ZHD* HD_p = new ZHD(l1, l2, A2, B2);
    start = std::chrono::high_resolution_clock::now();
    std::cout << HDparallel(HD_p->Az, HD_p->Bz, l1, l2, 1, tmp) << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "ZHD_p1 Execution time: " << duration.count() << " ms // " << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    
    ZHD* HD_p2 = new ZHD(l1, l2, A3, B3);
    start = std::chrono::high_resolution_clock::now();
    std::cout << HDparallel(HD_p2->Az, HD_p2->Bz, l1, l2, 3, tmp) << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "ZHD_p2 Execution time: " << duration.count() << " ms // " << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    */
    
    MyHD2* HD2 = new MyHD2(l1, l2, A4, B4);
    start = std::chrono::high_resolution_clock::now();
    end = std::chrono::high_resolution_clock::now();
    HD2->AmapB();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Merge time: " << duration.count() << " ms" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::cout << HD2->EarlyBreakDirectedHD() << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "MyHD Execution time: " << duration.count() << " ms // " << HD2->count << "count" << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    /*
    //gpuAlloc(l1, l2);
    MyHD2* HD2_p = new MyHD2(l1, l2, A5, B5);
    start = std::chrono::high_resolution_clock::now();
    end = std::chrono::high_resolution_clock::now();
    HD2_p->AmapB();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "MyHD_p Merge time: " << duration.count() << " ms" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    std::cout << HDparallel(HD2_p->A, HD2_p->B, l1, l2, 2, HD2_p->Aloc) << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms // " << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    */
    EBHD* HD3 = new EBHD(l1, l2, A, B);
    start = std::chrono::high_resolution_clock::now();
    std::cout << HD3->EarlyBreakDirectedHD() << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "EBHD Execution time: " << duration.count() << " ms // " << HD3->count << "count" << std::endl;

    std::cout << "-------------------------------------------------" << std::endl;
    /*
    EBHD* HD3_p = new EBHD(l1, l2, A, B);
    start = std::chrono::high_resolution_clock::now();
    HD3_p->Excluding();
    HD3_p->Randomize();
    std::cout << HDparallel(HD3_p->E, HD3_p->B, l1, l2, 0, tmp) << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "EBHD_p Execution time: " << duration.count() << " ms" << std::endl;
    */
    

//    113,826,289
// 42,060,145,691
//     80,555,510
// 37,450,969,980
}
