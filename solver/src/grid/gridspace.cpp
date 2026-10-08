/*! \file gridspace.cpp
 * \author Yu Zhou (zhouyu@live.com)
 */
#include "gridspace.h"
#include <stdexcept>

namespace DRL {

namespace GRID {

using std::invalid_argument;

GridSpace::GridSpace(std::size_t d, std::size_t q)
    : D(d) ,Q(initQ(d, q))
    , qAttr(initQAttr(d, q))
    , Qi   (qAttr.i   )
    , Qp   (qAttr.p   )
    , Qnorm(qAttr.norm)
    , Qop  (qAttr.op  )
    , Qhalf(qAttr.half)
    , Qadj (qAttr.adj )
    , Qrel (qAttr.rel )
    , QrelD(qAttr.reld)
    , Qprep(qAttr.perp)
    , Qaxis(qAttr.axis)
    , Qsign(qAttr.sign)
{

}

size_t GridSpace::initQ(size_t d, size_t q)
{
    size_t maxQ {static_cast<size_t>(std::pow(3, d))};
    return q > maxQ || q == 0 ? maxQ : q;
}

GridSpace::QAttrGroup GridSpace::initQAttr(size_t d, size_t q)
{
    const double ns {sqrt(2)};
    const double nc {sqrt(3)};
    const double sq {ns / 2};
    const double cb {nc / 3};
    QAttrGroup qa;

    switch (d) {
    case 1:         // one dimension
        switch (q) {
        case 3:     // D1Q3
            qa.i    = {{0}, {1}, {-1}};   // 0,1,2
            qa.p    = {{0}, {1}, {-1}};
            qa.norm = {0, 1, 1};
            qa.op   = {0, 2, 1};
            qa.half = {1};
            qa.adj  = {{0}, {1}, {2}};
            qa.rel  = {{}, {1}, {2}};         // x+, x-
            qa.reld = {{1,2}};
            qa.perp = {{}, {0}, {0}};
            qa.axis = {0, 0, 0};
            qa.sign = {0, 1, -1};
            break;
        default:
            throw invalid_argument("GridSpace::initQAttr() - Invalid number of"
                                   " directions.");
        }
        break;

    case 2:         // two dimensions
        switch (q) {
        case 5:     // D2Q5
            qa.i    = {{0, 0}, {0, 1}, {1, 0}, { 0,-1}, {-1, 0}}; // 0-4
            qa.p    = {{0, 0}, {0, 1}, {1, 0}, { 0,-1}, {-1, 0}};
            qa.norm = {0, 1, 1, 1, 1};
            qa.op   = {0, 3, 4, 1, 2};
            qa.half = {1, 2};
            qa.adj  = {{0}, {1}, {2}, {3}, {4}};
            qa.rel  = {{}, {1}, {2}, {3}, {4}};                   // x+,y+,x-,y-
            qa.reld = {{2, 4}, {1, 3}};                           // y, x
            qa.perp = {{}, {0, 2, 4}, {0, 1, 3}, {0, 2, 4}, {0, 1, 3}};
            qa.axis = {0, 1, 0,  1,  0};
            qa.sign = {0, 1, 1, -1, -1};
            break;
        case 9:     // D2Q9
            qa.i    = {{0, 0}, {0, 1}, {1, 0}, { 0,-1}, {-1, 0},  // 0-4
                       {1, 1}, {1,-1}, {-1,-1}, {-1, 1}};         // 5-8
            qa.p    = {{0, 0}, {0, 1}, {1, 0}, { 0,-1}, {-1, 0},
                       {sq, sq}, {sq, -sq}, {-sq, -sq}, {-sq, sq}};
            qa.norm = {0, 1, 1, 1, 1, ns, ns, ns, ns};
            qa.op   = {0, 3, 4, 1, 2, 7, 8, 5, 6};
            qa.half = {1, 2, 5, 6};
            qa.adj  = {{0}, {1}, {2}, {3}, {4},
                       {5, 1, 2}, {6, 2, 3}, {7, 3, 4}, {8, 1, 4}};
            qa.rel  = {{}, {1, 5, 8}, {2, 5, 6}, {3, 6, 7}, {4, 7, 8}};
            qa.reld = {{2, 5, 6, 4, 7, 8}, {1, 5, 8, 3, 6, 7}};   // y, x
            qa.perp = {{}, {0, 2, 4}, {0, 1, 3}, {0, 2, 4}, {0, 1, 3}};
            qa.axis = {0, 1, 0,  1,  0};
            qa.sign = {0, 1, 1, -1, -1};
            break;
        default:
            throw invalid_argument("GridSpace::initQAttr() - Invalid number of"
                                   " directions.");
        }
        break;

    case 3:         // three dimensions
        switch (q) {
        case 19:    // D3Q19
            qa.i    = {{0, 0, 0}, { 0, 0, 1}, { 0, 0,-1},               // 0-2
                       {0, 1, 0}, { 0,-1, 0}, { 1, 0, 0}, {-1, 0, 0},   // 3-6
                       {0, 1, 1}, { 0, 1,-1}, { 0,-1, 1}, { 0,-1,-1},   // 7-10
                       {1, 0, 1}, { 1, 0,-1}, {-1, 0, 1}, {-1, 0,-1},   // 11-14
                       {1, 1, 0}, { 1,-1, 0}, {-1, 1, 0}, {-1,-1, 0}};  // 15-18
            qa.p    = {{0, 0, 0}, { 0, 0, 1}, { 0, 0,-1},
                       {0, 1, 0}, { 0,-1, 0}, { 1, 0, 0}, {-1, 0, 0},
                       {0, sq, sq}, { 0, sq,-sq}, { 0,-sq, sq}, { 0,-sq,-sq},
                       {sq, 0, sq}, { sq, 0,-sq}, {-sq, 0, sq}, {-sq, 0,-sq},
                       {sq, sq, 0}, { sq,-sq, 0}, {-sq, sq, 0}, {-sq,-sq, 0}};
            qa.norm = {0, 1, 1, 1, 1, 1, 1,
                       ns, ns, ns, ns, ns, ns, ns, ns, ns, ns, ns, ns};
            qa.op   = { 0,  2,  1,  4,  3,  6,  5,
                       10,  9,  8,  7,
                       14, 13, 12, 11,
                       18, 17, 16, 15};
            qa.half = {1, 3, 5, 7, 8, 11, 12, 15, 16};
            qa.adj  = {{0}, {1}, {2}, {3}, {4}, {5}, {6},
                       { 7, 1, 3}, { 8, 2, 3}, { 9, 1, 4}, {10, 2, 4},
                       {11, 1, 5}, {12, 2 ,5}, {13, 1, 6}, {14, 2, 6},
                       {15, 3, 5}, {16, 4, 5}, {17, 3, 6}, {18, 4, 6}};
            qa.rel  = {{},
                       {1,  7,  9, 11, 13}, {2,  8, 10, 12, 14},  // x+,x-
                       {3,  8,  7, 15, 17}, {4,  9, 10, 16, 18},  // y+,y-
                       {5, 11, 12, 15, 16}, {6, 13, 14, 17, 18}}; // z+,z-
            qa.reld = {{5, 11, 12, 15, 16,  6, 13, 14, 17, 18},   // z
                       {3,  8,  7, 15, 17,  4,  9, 10, 16, 18},   // y
                       {1,  7,  9, 11, 13,  2,  8, 10, 12, 14}};  // x
            qa.perp = {{},
                       {0, 3, 4, 5, 6, 15, 16, 17, 18},           // x+
                       {0, 3, 4, 5, 6, 15, 16, 17, 18},           // x-
                       {0, 1, 2, 5, 6, 11, 12, 13, 14},           // y+
                       {0, 1, 2, 5, 6, 11, 12, 13, 14},           // y-
                       {0, 1, 2, 3, 4,  7,  8,  9, 10},           // z+
                       {0, 1, 2, 3, 4,  7,  8,  9, 10}};          // z-
            qa.axis = {0, 2,  2, 1,  1, 0,  0};
            qa.sign = {0, 1, -1, 1, -1, 1, -1};
            break;
        case 27:    // D3Q27
            qa.i    = {{ 0, 0, 0}, { 0, 0, 1}, { 0, 0,-1},              // 0-2
                       { 0, 1, 0}, { 0,-1, 0}, { 1, 0, 0}, {-1, 0, 0},  // 3-6
                       { 0, 1, 1}, { 0,-1, 1}, { 0,-1,-1}, { 0, 1,-1},  // 7-10
                       { 1, 1, 0}, {-1, 1, 0}, {-1,-1, 0}, { 1,-1, 0},  // 11-14
                       { 1, 0, 1}, {-1, 0, 1}, {-1, 0,-1}, { 1, 0,-1},  // 15-18
                       { 1, 1, 1}, { 1,-1, 1}, { 1,-1,-1}, { 1, 1,-1},  // 19-22
                       {-1, 1, 1}, {-1,-1, 1}, {-1,-1,-1}, {-1, 1,-1}}; // 23-26
            qa.p    = {{ 0, 0, 0}, { 0, 0, 1}, { 0, 0,-1},
                       { 0, 1, 0}, { 0,-1, 0}, { 1, 0, 0}, {-1, 0, 0},
                       { 0, sq, sq}, { 0, -sq,sq}, {  0,-sq,-sq}, {  0, sq,-sq},
                       { sq,sq,  0}, {-sq, sq, 0}, {-sq,-sq,  0}, { sq,-sq,  0},
                       { sq, 0, sq}, {-sq,  0,sq}, {-sq,  0,-sq}, { sq,  0,-sq},
                       { cb,cb, cb}, { cb,-cb,cb}, { cb,-cb,-cb}, { cb, cb,-cb},
                       {-cb,cb, cb}, {-cb,-cb,cb}, {-cb,-cb,-cb}, {-cb,cb,-cb}};
            qa.norm = {0, 1, 1, 1, 1, 1, 1,
                       ns, ns, ns, ns, ns, ns, ns, ns, ns, ns, ns, ns,
                       nc, nc, nc, nc, nc, nc, nc, nc};
            qa.op   = { 0,  2,  1,  4,  3,  6,  5,  9, 10,  7,  8, 13, 14, 11,
                        12, 17, 18, 15, 16, 25, 26, 23, 24, 21, 22, 19, 20};
            qa.half = {1, 3, 5, 7, 10, 11, 12, 15, 16, 19, 22, 23, 26};
            qa.adj  = {{0}, {1}, {2}, {3}, {4}, {5}, {6},
                       { 7, 1, 3}, { 8, 1, 4}, { 9, 2, 4}, {10, 2, 3},
                       {11, 3, 5}, {12, 3 ,6}, {13, 4, 6}, {14, 4, 5},
                       {15, 1, 5}, {16, 1, 6}, {17, 2, 6}, {18, 2, 5},
                       {19,  7, 11, 15}, {20,  8, 14, 15},
                       {21,  9, 14, 18}, {22, 10, 11, 18},
                       {23,  7, 12, 16}, {24,  8, 13, 16},
                       {25,  9, 13, 17}, {26, 10, 12, 17}};
            qa.rel  = {{},
                       {1,  7,  8, 15, 16, 19, 20, 23, 24},    // 1 - x+
                       {2,  9, 10, 17, 18, 21, 22, 25, 26},    // 2 - x-
                       {3,  7, 10, 11, 12, 19, 22, 23, 26},    // 3 - y+
                       {4,  8,  9, 13, 14, 20, 21, 24, 25},    // 4 - y-
                       {5, 11, 14, 15, 18, 19, 20, 21, 22},    // 5 - z+
                       {6, 12, 13, 16, 17, 23, 24, 25, 26}};   // 6 - z-
            qa.rel  = {{5, 11, 14, 15, 18, 19, 20, 21, 22,     // z
                        6, 12, 13, 16, 17, 23, 24, 25, 26},
                       {3,  7, 10, 11, 12, 19, 22, 23, 26,     // y
                        4,  8,  9, 13, 14, 20, 21, 24, 25},
                       {1,  7,  8, 15, 16, 19, 20, 23, 24,     // x
                        2,  9, 10, 17, 18, 21, 22, 25, 26}};
            qa.perp = {{},
                       {0,  3,  4,  5,  6, 11, 12, 13, 14},    // 1 - x+
                       {0,  3,  4,  5,  6, 11, 12, 13, 14},    // 2 - x-
                       {0,  1,  2,  5,  6, 15, 16, 17, 18},    // 3 - y+
                       {0,  1,  2,  5,  6, 15, 16, 17, 18},    // 4 - y-
                       {0,  1,  2,  3,  4,  7,  8,  9, 10},    // 5 - z+
                       {0,  1,  2,  3,  4,  7,  8,  9, 10}};   // 6 - z-
            qa.axis = {0, 2,  2, 1,  1, 0,  0};
            qa.sign = {0, 1, -1, 1, -1, 1, -1};
            break;
        default:
            throw invalid_argument("GridSpace::initQAttr() - Invalid number of"
                                   " directions.");
        }
        break;

    default:
        throw invalid_argument("GridSpace::initQAttr() - Invalid number of"
                               " dimensions.");
    }
    return qa;
}

} // namespace GRID

} // namespace DRL
