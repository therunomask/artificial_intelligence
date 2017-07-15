#ifndef TOPLAYER_H
#define TOPLAYER_H

#include"layer.h"

class top_layer : public layer{
public:
    top_layer(size_t Number_of_Column_per_layer, size_t Number_of_Cells_per_Column, brain& pBrain);
//pointer to upper level stays NULL
    //Three_CellActivityList is really a two_cellactivityList in this case
    void Three_CellListUpdater(void);

};



#endif // TOPLAYER_H
