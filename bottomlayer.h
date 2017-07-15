#ifndef BOTTOM_LAYER_H
#define BOTTOM_LAYER_H

#include"layer.h"
#include"magicnumbers.h"
#include<vector>
#include<iostream>

class bottom_layer : public layer{
public:
    //feed function pointer to FindBestColumn()!
    bottom_layer(size_t Number_of_Column_per_layer, size_t Number_of_Cells_per_Column,brain& pBrain,std::vector<bool>(*sensoryinput)(size_t time));
    //FindBestcolumns() is the best place to redefine dynamic of lowest layer by
    //model specific behavior!

    //not yet implemented; write 1. constructor, 2. FindBestcolumn, 3. Three_CellListUpdater
    void FindBestColumns();


    std::vector<bool>(*external_input)(size_t time);
    void ConnectedSynapsesUpdate()
    {
        std::cout<<"lowest layer should not use ConnectedSynapsesUpdate! \n";
        std::abort();
    }//not necessary in lowest layer

    void Three_CellListUpdater();
};


#endif // BOTTOM_LAYER_H
