#pragma once

#define UINT_MAX 4294967295
#define INT_MAX 2147483647

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