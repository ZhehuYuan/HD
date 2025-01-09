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

#define UINT_MAX 4294967295
#define INT_MAX 2147483647

using namespace std;

struct Point
{
    int l;
    int w;
    int h;

    Point() {
        l = INT_MAX;
        w = INT_MAX;
        h = INT_MAX;
    }

    Point(int a, int b, int c) {
        l = a;
        w = b;
        h = c;
    }

    Point& operator=(const Point& other) {
        if (&other == this)return *this;

        l = other.l;
        w = other.w;
        h = other.h;

        return *this;
    }

    bool operator==(const Point& other) const {
        return l == other.l && w == other.w && h == other.h;
    }

    bool operator<(const Point& other) const {
        int tmp1 = l + w * 10000 + h * 100000000;
        int tmp2 = other.l + other.w * 10000 + other.h * 100000000;
        return tmp1 < tmp2;
    }

    void swap(Point& other) {
        if (&other == this)return;

        std::swap(l, other.l);
        std::swap(w, other.w);
        std::swap(h, other.h);
    }

    double Distance(Point& other) {
        int x = l - other.l;
        int y = w - other.w;
        int z = h - other.h;
        return sqrt(x * x + y * y + z * z);
    }
};

struct PointCMin
{
    Point p;
    double cmin;

    PointCMin() {
        return;
    }

    PointCMin(Point _p, double init) {
        p = _p;
        cmin = init;
    }

    PointCMin& operator=(const PointCMin& other) {
        p = other.p;
        cmin = other.cmin;
        return *this;
    }

    bool operator==(const PointCMin& other) const {
        return cmin == other.cmin;
    }

    bool operator<(const PointCMin& other) const {
        return cmin < other.cmin;
    }
    
    bool operator>(const PointCMin& other) const {
        return cmin > other.cmin;
    }

    int getValue() {
        return p.l + p.w + p.h;
    }
};

struct MyPoint {
    Point p;
    int mortonCode;

    MyPoint() {
        return;
    }

    MyPoint(Point _p, int m) {
        p = _p;
        mortonCode = m;
    }

    bool operator==(const MyPoint& other) const {
        return mortonCode == other.mortonCode;
    }

    bool operator<(const MyPoint& other) const {
        return mortonCode < other.mortonCode;
    }
};

class EBHD
{
public:
    int l1, l2;
    Point* A, * B, * E;

    EBHD(int _l1, int _l2, Point* _A, Point* _B){
        l1 = _l1;
        l2 = _l2;
        B = _B;
        A = _A;
        E = new Point[l1];
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

    ZHD(int _l1, int _l2, Point* _A, Point* _B) {
        l1 = _l1;
        l2 = _l2;
        B = _B;
        A = _A;
        Z1 = new PointCMin[l1];
        Z2 = new PointCMin[l2];
        Az = new Point[l1];
        Bz = new Point[l2];

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

    MyHD2(int _l1, int _l2, Point* _A, Point* _B) {
        A = _A;
        B = _B;
        l1 = _l1;
        l2 = _l2;
        Az = new MyPoint[l1];
        Bz = new MyPoint[l2];
        Aloc = new int[l1];
    
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
    srand(time(nullptr));

    Point* A, * B, * A2, * B2, * B3, * A4, * B4;
    int l1, l2;

    l1 = 1000000;
    l2 = 1000;

    int x = 350;
    int y = 350;
    int z = 350;

    A = new Point[l1];
    B = new Point[l2];

    A2 = new Point[l1];
    B2 = new Point[l2];

    B3 = new Point[l2];

    A4 = new Point[l1];
    B4 = new Point[l2];

    for (int i = 0; i < l1; i++) {
        A[i] = { rand() % x, rand() % y, rand() % z };
        A2[i] = A[i];
        A4[i] = A[i];
    }
    for (int i = 0; i < l2; i++) {
        B[i] = { rand() % x, rand() % y, rand() % z };
        B2[i] = B[i];
        B3[i] = B[i];
        B4[i] = B[i];
    }

    double cmin = numeric_limits<double>::max();
    PointCMin* A3 = new PointCMin[l1];

    for (int i = 0; i < l1; i++) {
        A3[i] = { A[i], cmin };
    }

    ZHD* HD2 = new ZHD(l1, l2, A2, B2);
    auto start = std::chrono::high_resolution_clock::now();
    std::cout << HD2->EarlyBreakDirectedHD() << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;

    MyHD2* HD4 = new MyHD2(l1, l2, A4, B4);
    start = std::chrono::high_resolution_clock::now();
    std::cout << HD4->EarlyBreakDirectedHD() << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    
    EBHD* HD = new EBHD(l1, l2, A, B);
    start = std::chrono::high_resolution_clock::now();
    std::cout << HD->EarlyBreakDirectedHD() << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
    
    MyHD* HD3 = new MyHD(l1, l2, A3, B3, x, y, z);
/*
    start = std::chrono::high_resolution_clock::now();
    std::cout << HD3->EarlyBreakDirectedHD() << std::endl;
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Execution time: " << duration.count() << " ms" << std::endl;
*/
}
