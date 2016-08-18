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



class LowSyn{//list of synapses to lower level of each column
    std::vector<double> conds;//connectedness
    std::vector<size_t> ConLow;//connections to lower layer

};

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
    bool active;
    std::vector<std::vector<LowSyn>> PotSyn;//list of potential synapses
    std::vector<std::vector<LowSyn*>> ConSyn;//list of pointers to "connected" synapses

public:
    bool feed_input(std::vector<bool> input);//returnes updated activation
};


class layer{//write a derived class lowest_layer
private:
    const layer* p_lower_level;//pointer to lower layer to receive input
    const size_t Num_Columns;
    std::vector<size_t> ActColumns;//active columns; maybe turn into array of active
                                //columns of the last few time steps
    std::vector<std::vector<bool>> ActLog;//activity log
    //records the activity of all columns over the last few steps
    std::vector<column> ColumnList;//columns in the layer
    float ConThr;//threshhold for synapses in order to be "connected"

    //maybe reconsider connectedness mechanism
    double CondsInc;//connectedness increment
    double CondsDec;//connectedness decrement

    //stuff for inhibition
    size_t DesLocAct;//desired local activity
    //a measure of how many columns should be active on average
    //probably unnecessary consider changing inhibition
    size_t InhiRad;//inhibition radius used to update neighbors
    std::vector<std::vector<size_t>> neighbors;//list of neighbors for each column, consider changing inhibition
    size_t MinOverlap;//minimum overlap
    //end of inhibition stuff

    std::vector<double> boost();//vector of boost values increases
                        //activity of columns which are not active enough
public:
    std::vector<bool> current_activation( void );//triggers prediction of lower level
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
