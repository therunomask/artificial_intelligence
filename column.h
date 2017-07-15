#ifndef COLUMN_H
#define COLUMN_H
#include"layer.h"
#include"cell.h"
#include"segment.h"
#include"magicnumbers.h"
#include<vector>


class layer;

class column{//contains connections to input and list of cells
    //that contain lateral connections to predict activation of the column
private:

    double Overlap_Average;
    static const int MinOverlap=static_cast<int>(minimal_overlap_under_consideration);


public:
    column(layer& layer_to_belong_to, size_t Number_of_Cells_per_Column);
    column(const column& dummycolumn);
    constexpr static double Average_Exp= column_geometric_factor; //Take average of overlap
              //correspons to 65% of whole value was determined in the last 100 steps
    std::vector<std::pair<column*,double>> ConnectedSynapses;//list of pointers to "connected" synapses and connection strength

    layer& MotherLayer;
    std::vector<cell> CellList;//List of cells in the colummn

    std::vector<bool> active;//activity of feed forward input
    bool expect;

    double ActivityLog;//running average via geometric series

    double boosting;//boost value increases
                        //activity of columns which are not active enough


    double feed_input(void);//computes activation caused by input
    //connections count only as "connected" or "not connected" no further weights

    double tell_overlap_average(void){//return running average overlap
        return Overlap_Average;
    }
    segment* BestMatchingSegmentInColumnActivateCells(size_t time);


    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);

};



#endif // COLUMN_H
