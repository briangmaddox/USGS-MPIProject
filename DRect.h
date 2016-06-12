/***************************************************************************
                          DRect.h  -  description
                             -------------------
    begin                : Thu Oct 12 2000
    copyright            : (C) 2000 by Chris Bilderback
    email                : cbilderback@usgs.gov
 ***************************************************************************/


#ifndef DRECT_H
#define DRECT_H
class DRect
{
public:
    //main constructor for the DRect class
    DRect() : left(0), right(0), top(0), bottom(0)
    {}

    double left, right, top, bottom;
};
#endif //DRECT_H
