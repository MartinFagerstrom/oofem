#include "delaunay.h"
#include "flotarry.h"
#include "intarray.h"
#include "alist.h"
#include <math.h>
#include "geometry.h"
#include "node.h"

bool Delaunay :: colinear(FloatArray *p1, FloatArray *p2, FloatArray *p3) {
    double dist = p1->at(1) * ( p2->at(2) - p3->at(2) ) + p2->at(1) * ( p3->at(2) - p1->at(2) ) +
    p3->at(1) * ( p1->at(2) - p2->at(2) );
    // the tolerance probably needs a setter
    if ( dist < 0.0001 && dist > ( -1 ) * 0.0001 ) {
        return true;
    } else {
        return false;
    }
}

void Delaunay :: printTriangles(AList< Triangle > *triangles) {
    for ( int i = 1; i <= triangles->giveSize(); i++ ) {
        triangles->at(i)->printYourself();
    }
}

bool Delaunay :: isInsideCC(FloatArray *p, FloatArray *p1,  FloatArray *p2,  FloatArray *p3) {
    FloatArray *nodesCopy1 = new FloatArray(*p1);
    FloatArray *nodesCopy2 = new FloatArray(*p2);
    FloatArray *nodesCopy3 = new FloatArray(*p3);
    Triangle *tr = new Triangle(nodesCopy1, nodesCopy2, nodesCopy3);
    double r = tr->getRadiusOfCircumCircle();
    FloatArray circumCenter;
    tr->computeCenterOfCircumCircle(circumCenter);
    double distance = circumCenter.distance(p);
    delete tr;
    if ( distance < r ) {
        return true;
    } else {
        return false;
    }
}

void Delaunay :: triangulate(AList< FloatArray > *vertices, AList< Triangle > *triangles) {
    /// 4th order algorithm - four loops, only for testing puropses
    int n = vertices->giveSize();
    int count = 0;
    /// small shift of vertices
    for ( int i = 1; i <= n; i++ ) {
        vertices->at(i)->at(1) += vertices->at(i)->at(1) * 0.000001 * double ( rand() ) / RAND_MAX;
        vertices->at(i)->at(2) += vertices->at(i)->at(2) * 0.000001 * double ( rand() ) / RAND_MAX;
    }

    for ( int i = 1; i <= n; i++ ) {
        for ( int j = i + 1; j <= n; j++ ) {
            for ( int k = j + 1; k <= n; k++ ) {
                bool isTriangle = true;
                if ( colinear( vertices->at(i), vertices->at(j),
                              vertices->at(k) ) ) {
                    isTriangle = false;
                } else {
                    for ( int a = 1; a <= n; a++ ) {
                        if ( a != i && a != j && a != k ) {
                            // checks whether a point a is inside a circumcircle of a triangle ijk
                            if ( isInsideCC( vertices->at(a), vertices->at(i), vertices->at(j),
                                            vertices->at(k) ) ) {
                                isTriangle = false;
                                break;
                            }
                        }
                    }
                }

                if ( isTriangle ) {
                    count++;
                    FloatArray *p1 = new FloatArray();
                    * p1 = * vertices->at(i);
                    FloatArray *p2 = new FloatArray();
                    * p2 = * vertices->at(j);
                    FloatArray *p3 = new FloatArray();
                    * p3 = * vertices->at(k);
                    Triangle *triangle = new Triangle(p1, p2, p3);
                    if ( !triangle->isOrientedAnticlockwise() ) {
                        triangle->changeToAnticlockwise();
                    }

                    triangles->put(count, triangle);
                }
            }
        }
    }
}

