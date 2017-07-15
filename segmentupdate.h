#ifndef SEGMENTUPDATE_H
#define SEGMENTUPDATE_H

#include<vector>
#include"cell.h"

class cell;

class SegmentUpdate{
public:


    SegmentUpdate( std::vector<cell*> active_cells, size_t countdown );

    std::vector<cell*> active_cells;//points to adresses in segment.Synaps
    size_t timer;


};


#endif // SEGMENTUPDATE_H
