#ifndef LAYER_H
#define LAYER_H

#include"brain.h"
#include"column.h"
#include"magicnumbers.h"


class brain;
class layer{//write a derived class lowest_layer; output_layer
    //write constructor!!!!
private:


    constexpr static double FractionOfActiveColumns=active_pillers_per_pillar;
    size_t DesiredLocalActivity;
    //             v change to std::vector<column&>

    //initialize to (Num_Columns,100)!
    //running average of activity of cells
    constexpr static double CondsInc=Learning_Increment_spatial;//connectedness increment
    constexpr static double MaximumConnectedness=Maximum_Connectedness;
    constexpr static double AverageOverlapMin = Average_Overlap_lower_boundary;
    //find meaningful values!!



public:

    layer(size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column, brain& pBrain);
    layer(const layer& dummylayer);

    std::vector<column*> ActColumns;//active columns; maybe turn into array of active

    std::vector<column*> TempActColumns;//active columns; maybe turn into array of active

    brain& MotherBrain;
    const size_t Num_Columns;
    layer*  p_lower_level;//pointer to lower layer to receive input
    layer*  p_upper_level;//pointer to layer above current one

    std::vector<column> ColumnList;//columns in the layer

    std::vector<cell*> CellActivityList;//update as fast as possible

    std::vector<std::vector<cell*>> Three_CellActivityList;//update parallel
                                                            //outer vector for timesteps

    std::vector<cell*> CellUpdateList;
    std::vector<cell*> PendingActivity;
    std::vector<cell*> PendingExpectation;

    void virtual FindBestColumns(void);
    void ActiveColumnUpdater(void);
    void virtual ConnectedSynapsesUpdate(void);
    void ActivityLogUpdate(void);
    void Do_SegmentUpdate(void);

    void CellExpectInitiator(void);
    void CellUpdater(void);
    void CellLearnInitiator(void);
    void virtual Three_CellListUpdater(void);
    void forgetting(void);
    std::vector<column *> ColumnSynapseAdding(size_t ColumnsToFind);

    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);

};

#endif // LAYER_H
