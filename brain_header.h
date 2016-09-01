#ifndef BRAIN_HEADER_H
#define BRAIN_HEADER_H

#include<vector>
#include<deque>

class brain;
class layer;
class column;
class LowSyn;
class cell;
class segment;



class segment{
private:
    std::deque<const cell*> CellAddr;//pointer to adresses of cells in the segment
    std::deque<double> CellCon;//connectedness values of the synapses in the segment
    //length of those vectors must be the same!!
    const static double InitCon;//initial connectedness for new
            //cells in the segment
    bool EndOfSeq;//true if predictions due to this segment
            //should activate the column in question in the next
            //timestep rather than predicting prediction of activation
public:
    void AddCell(const cell*  const newcell){//adds new cell to segment with
            //standard connectedness close to disconnection
        CellAddr.push_back(newcell);
        CellCon.push_back(InitCon);
    }
    void UpdateCon(std::vector<const cell*> winners );//increase connectedness
        //of winners, decrease connectedness
        //of all other cells in the segment, disconnect them if
        //connection too weak.
};
class cell{
    std::vector<segment> SegList;//list of segments (net of horizontal
            // connections) of this cell
    bool active;//representing input
    bool expect;//predicts input due to past experience and dendrite information
};

class column{//contains connections to input and list of cells
    //that contain lateral connections to predict activation of the column
private:
    std::vector<cell> CellList;//List of cells in the colummn
    constexpr static double Average_Exp= 99/100.0; //Take average of overlap
              //over last 20 values
    double Overlap_Average=100;//should be active every 5th?! time
                //pick meaningful average value as initial value
    static const int MinOverlap=123456789;//find meaningful value

    std::vector<size_t> ConnectedSynapses;//list of indices of "connected" synapses

public:
    std::vector<double> connectedness;//connectedness

    double feed_input(const std::vector<bool> &input);//computes activation caused by input
    //connections count only as "connected" or "not connected" no further weights
    double boosting;//boost value increases
                        //activity of columns which are not active enough
    double tell_overlap_average(void){//return running average overlap
        return Overlap_Average;
    }

};


class layer{//write a derived class lowest_layer
private:
     layer* const p_lower_level;//pointer to lower layer to receive input
    const size_t Num_Columns;
    std::vector<size_t> ActColumns;//active columns; maybe turn into array of active
                                //columns of the last few time steps
    std::vector<double> ActivityLog;//activity log
    //initialize to (Num_Columns,100)!
    //running average of activity of cells
    const double AverageExp=0.99;//parameter in running average;
    // find meaningful  value  ^  !
    std::vector<column> ColumnList;//columns in the layer
    float ConThr;//threshhold for synapses in order to be "connected"

    //maybe reconsider connectedness mechanism
    const double CondsInc=0.01;//connectedness increment
    const double CondsDec=0.01;//connectedness decrement
    constexpr static double AverageOverlapMin = 0.01;
    constexpr static double SpecialOverlapIncrement= 0.1;
    //find meaningful values!!

    //stuff for inhibition
    size_t DesiredLocalActivity;//desired local activity;
    //less than ColumnList.size()!!
    //a measure of how many columns should be active on average
    //probably unnecessary consider changing inhibition
    size_t InhiRad;//inhibition radius used to update neighbors
    std::vector<std::vector<size_t>> neighbors;//list of neighbors for each column, consider changing inhibition
    size_t MinOverlap;//minimum overlap
    //end of inhibition stuff

    void Update_Synapses(void);

public:
     std::vector<bool> activation_learning( void );//triggers prediction of lower level
                            //computes current activation-
                            //distribution of culumns. triggers same function of
                            //lower layer to use as input
                        //somehow propagate expected activation of lower level to it!!
    virtual std::vector<bool> current_prediction( void );//triggers activation of same level
        //is "virtual" because lowest_layer needs to redefine this function
    
    //reconsider private/public position of memberfunctions!
    double overlap(std::vector<bool> input,size_t column);//maybe returning a vector and omitting second arguent is smarter
    //computes overlap of a particular column with given input

    void updateBoost(std::vector<bool> input);
    //checks if activity log agrees with the wanted average activity

    void updateLow(std::vector<bool> input);
    //updates all the synapses of all the columns to the lower level
        //due to input from lower level
        //all synapses in active columns get reinforced if they were
        //active right before activation of the corresponding column
        //all the nonactive synapses in active columns are decremented
        //---- the synapses which belong to cells whose overlap with input
        //is lower than proper, are uniformly reinforced as well


    void updateSegPred();//needs activation pattern of current timestep
            //updates the segments of the cells of all columns
            //randomly adds cells to segments with
            //nonpredicting synapses and deletes these
            // due to prediction ->
            //case 1: correct prediction:
            //add random new segments and reinforce
            //correctly predicting segments, degrade other segments
            //case 2: unpredicted activation:
            //1. find cell within active column whose segments
            // fit activation pattern of the layer of last timestep best
            //2. designate this segment to be predicting activation and reinforce it
};


class brain{
private:
    const size_t NumLevels;
    const std::vector<layer> ListLevels;

};

#endif // BRAIN_HEADER_H
