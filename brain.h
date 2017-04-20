#ifndef BRAIN_HEADER_H
#define BRAIN_HEADER_H

#include<vector>
#include<deque>
#include<unordered_set>
#include<iostream>
#include <fstream>
#include <mutex>

class brain;
class layer;
class column;
class LowSyn;
class cell;
class segment;
class debughelper;




//change objects:
//SegmentUpdateList -> saves active Cells; which Segment; timer
/*
 *
 *
 *
 * fix first:       problems occur when updating segments, because there may be pointer to segments
 *                      , which don't exist anymore. Probably there still lurk other problems around.
 *
 *fix layer::forgetting() ->uses segment::operator= , which we don't want and made thrownig errorr
 *
 *
 *in brain constructer: initiate stuff like threeCellActivityList randomly
 * instead of empty
 *
 * change boosting, perfect matching should be stable
 *
 *
 * do we want a whole sequence of n predicting segments in order to
 * predict in n timesteps, or do we allow for a gap in the sequence?
 *
 *
 *
 * check if things depend on a fixed segment number
 *
 *program crashes if more than 4 columns are active in the lowest layer
 *
 *
 *think of relevant variables that describe the workflow of our system
 *
 * successrate of columns ( column gets into expected state and is activated accordingly)
 * successrate of cells
 * rate of activation of cells/columns
 * rate of change in the synapses; columns as well as segments
 * rate of change of weights of synapses
 * distribution of weights of synapses
 * rate of blind synapse adding
 * average number of synapses per segment
 * time each function takes
 *
 *
 *
 *
 *

propagate confusion to higher layers

 */



/* Theoretical Aspects to be reconsidered:
 *
 *check multithreadding+flowchart
 *
 * EndOfSequence?
 *
 * Dopaminsystem?
 *
 * dynamically change parameters (e.g. number of columns)
 *
 * Is Inhibition radius a useful concept, or should one model a second inhibiting network?
 *
 *
 */

/*numbers to opmize
 *
 *
 * how to handle number of synapses in each segment?
 * 3* active columns per layer
 *
 * how to handle number of segments for each cell?
 * at most cells_per_column*3*#columns/#active_columns
 *
 * cells per columns?
 * ?3?
 *
 * columns per layer?
 * 50 * ?4?
 *
 * layers?
 * ?5?
 *
 * InitCon= ?0.5?
 * Learnincrement= ?0.07?
 *
 * active pillars per pillar: ?2%?
 *
 * ConThr=?0.3? <=> neccessary connectedness value of synapse of pillar
 *
 */

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//          optimize with respect to these numbers
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//general dimensions of the system
#define layers_per_brain                            4//8 with multithreadding
#define active_pillers_per_pillar                   0.02
#define cells_per_column                            3
#define pillars_per_layer                           4 / active_pillers_per_pillar
#define segments_per_cell                           cells_per_column/active_pillers_per_pillar/7
#define synapses_per_segment                        3*active_pillers_per_pillar*pillars_per_layer
//end of general dimensions

//synapse
#define initial_connectedness                       0.5
#define Minimal_sum_of_synapseweights_for_activity  0.3*cells_per_column*active_pillers_per_pillar*pillars_per_layer
#define Learning_Increment                          0.07

//column
#define column_geometric_factor                     99/100.0
#define minimal_overlap_under_consideration         0.3 * pillars_per_layer*active_pillers_per_pillar
#define Initial_Activity_log                        2.5 //(1/(1-column_geometric_factor) ^(1/active_pillers_per_pillar)))
#define Initial_Overlap_Average                     Initial_Activity_log* active_pillers_per_pillar*pillars_per_layer

//layer
#define Learning_Increment_spatial                  0.01
#define Learning_Decrement_spatial                  0.01
#define Average_Overlap_lower_boundary              0.01
#define Homogenous_Overlap_Increment                0.1
#define Forgetfulness                               0.01

//magic boosting function
#define maximum_boosting                            3.0


class segment{
private:
    static constexpr double InitCon=initial_connectedness;
            //cells in the segment

public:
    segment(cell& Cell_to_belong_to, size_t TempActivationCountdown);
    segment(const segment& dummysegment);
    segment operator=(const segment& dummysegment);
    ~segment(void){
        who_am_I();
        std::cout<<"and I will now be destroyed\n";
        throw std::invalid_argument("segment was deleted\n");
    }

    cell& MotherCell;                            //3*activeCollumns per layer
    static constexpr double MinSynapseWeightActivity=Minimal_sum_of_synapseweights_for_activity;
    std::deque< std::pair <cell*,double>> Synapse;//pointer to adresses of cells in the segment
    //and connectedness values of the synapses in the segment
    size_t ActivationCountdown;//1 if prediction's due to this segment
            //should activate the column in question in the next
            //timestep rather than predicting prediction of activation

    static constexpr bool PositiveLearning= true;
    static constexpr double LearnIncrement= Learning_Increment;//possibly change to harder
                     // punishment for larger segments

    void AddCell( cell*  const newcell){//adds new cell to segment with
            //standard connectedness close to disconnection
        double con=InitCon;
        Synapse.push_back(std::pair <cell*,double>(newcell,con));
        //CellCon.push_back(InitCon);
    }



    void BlindSynapseAdding(size_t t);
    std::vector<cell *> GetActiveCells();


    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);
};

class SegmentUpdate{
public:
    ~SegmentUpdate(void){
        //std::cout<<"deleting SegmentUpdate\n";
    }

    SegmentUpdate(segment* SegmentAddress, std::vector<cell*> active_cells );
    segment* SegmentAddress;
    std::vector<cell*> active_cells;//points to adresses in segment.Synaps
    size_t timer;

    void AdaptingSynapses(bool positive);//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
};

class cell{

public:
    cell(column& Column_to_belong_to);
    cell(const cell& dummycell);

    column& MotherColumn;
    std::deque<segment> SegList;//list of segments (net of horizontal
            // connections) of this cell
    std::vector<std::vector<segment*>> ActiveSegments;//list of active segments of
    //    ^ maybe not use std::vector                 // last timesteps
    std::vector<bool> active;//saves activity of last few timesteps
    std::vector<bool> expect;//predicts input due to past experience and dendrite information
    std::vector<SegmentUpdate> SegmentUpdateList;

    void UpdateActiveSegments(void);
    segment* BestSegmentInCell(size_t t);
    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);
};

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

    bool active;//activity of feed forward input
    bool expect;

    double ActivityLog;//running average via geometric series

    double boosting;//boost value increases
                        //activity of columns which are not active enough


    double feed_input(void);//computes activation caused by input
    //connections count only as "connected" or "not connected" no further weights

    double tell_overlap_average(void){//return running average overlap
        return Overlap_Average;
    }
    segment* BestMatchingSegmentInColumnActivateCells(void);


    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);

};

class layer{//write a derived class lowest_layer; output_layer
    //write constructor!!!!
private:


    constexpr static double FractionOfActiveColumns=active_pillers_per_pillar;
    size_t DesiredLocalActivity;
    //             v change to std::vector<column&>

    //initialize to (Num_Columns,100)!
    //running average of activity of cells
    constexpr static double CondsInc=Learning_Increment_spatial;//connectedness increment
    constexpr static double CondsDec=Learning_Decrement_spatial;//connectedness decrement
    constexpr static double AverageOverlapMin = Average_Overlap_lower_boundary;
    constexpr static double SpecialOverlapIncrement= Homogenous_Overlap_Increment;
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
    double ActivityLogUpdateFindMaxActivity(void);
    void virtual BoostingUpdate_StrenthenWeak(double MaxActivity);
    void Do_SegmentUpdate(void);

    void CellExpectInitiator(void);
    void CellUpdater(void);
    void CellLearnInitiator(void);
    void virtual Three_CellListUpdater(void);
    void forgetting(void);

    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);

};


class top_layer : public layer{
public:
    top_layer(size_t Number_of_Column_per_layer, size_t Number_of_Cells_per_Column, brain& pBrain);
//pointer to upper level stays NULL
    //Three_CellActivityList is really a two_cellactivityList in this case
    void Three_CellListUpdater(void);

};

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
    void BoostingUpdate_StrenthenWeak()
    {
        std::cout<<"lowest layer should not use BoostingUpdate_StrenthenWeak! \n";
        std::abort();
    }//not necessary in lowest layer
    void Three_CellListUpdater();
};

class debughelper{
private:
public:
    debughelper(void);
    ~debughelper(void);
    std::vector<std::vector<double>> success_column;
    std::vector<std::vector<double>> success_cell;
    std::vector<std::vector<double>> activation_column;
    std::vector<std::vector<double>> activation_cell;
    std::vector<std::vector<double>> avg_synapses_per_segment;
    std::ofstream Log;

    template<typename T>
    std::ofstream& operator<< ( T message);
    //std::ofstream& operator<< ( double bla);
    void tell(std::vector<std::vector<double> > *dummyvec);

};


class brain{
private:
    const size_t NumLevels;
    std::vector<layer> ListOfLevels;



public:
    bottom_layer LowestLayer;
    top_layer    HighestLayer;
    std::vector<layer*> AllLevels;
    size_t time;
    debughelper Martin_Luther;
    size_t max_activation_counter;
    bool max_activation_counter_change;
    bool multithreadding=false;


    brain(size_t Number_of_Levels, size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column,std::vector<bool>(*sensoryinput)(size_t time));
    friend void BrainConstructionHelper(brain& Init_brain, size_t Number_of_Levels, size_t Number_of_Column_per_Layer,size_t Number_of_Cells_per_Column/*,std::vector<bool>(*sensoryinput)(size_t time)*/);

    void update(void);

    void inventory(void);
    void pointerCheck(void);
    /*updates:
     * layer::actcolumns
     * layer::SegmentUpdateList
     * layer::CellActivityList
     * layer::Three_CellActivityList
     * column::Overlap_Average
     * column::ConnectedSynapses
     * column::active
     * column::ActivityLog
     * column::boosting
     * cell::ActiveSegments
     * cell::active
     * cell::expect
     * cell::learn
     * segment::active
     * segment::Synapse
     * segment::EndOfSeq
     * */

};

void UpdateInitialiser(layer* player);
void ThreadUpdater(layer* DummyLayer);

#endif // BRAIN_HEADER_H
