#ifndef BRAIN_HEADER_H
#define BRAIN_HEADER_H

#include<vector>
#include<deque>
#include<unordered_set>
#include<iostream>
#include <fstream>
#include <mutex>
#include"bottomlayer.h"
#include"layer.h"
#include"toplayer.h"
#include"debughelper.h"
#include"magicnumbers.h"


class brain;
class layer;
class column;
class cell;
class segment;
class debughelper;
class SegmentUpdate;





/*
 *
 * change initialization of cell synapses so that no two of the same segment point to the same cell
 *
 * disabled influence of expectation on activity of columns in feed_input()
 *
 * 2. expectation over time changes boosting
 * 3. impose L1-norm on connections between columns
 * 4. propagate information about wrong expectation to higher levels
 *
 *
 * do we want a whole sequence of n predicting segments in order to
 * predict in n timesteps, or do we allow for a gap in the sequence?
 *
 * move inventory and pointercheck to debughelper
 *
 *add multithreadding also to forgetting, threecelllistUpdater loop
 *
 * program crashes if more than 4 columns are active in the lowest layer
 * find out how to properly write custom destructors
 *
 * think of relevant variables that describe the workflow of our system
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
 * change Forgetfulness for each layer seperately
 * replace deque with list or forward_list
 *

propagate confusion to higher layers

 */



/* Theoretical Aspects to be reconsidered:
 *
 *check multithreadding+flowchart
 *
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







class brain{
private:
    const size_t NumLevels;




public:
    std::vector<layer> ListOfLevels;
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


};

void UpdateInitialiser(layer* player);
void ThreadUpdater(layer* DummyLayer);

#endif // BRAIN_HEADER_H
