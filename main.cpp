#include <iostream>
#include <vector>

class synapse{
    double conds;//connectedness
    std::vector<size_t> ConLow;//connections to lower layer
    std::vector<size_t> ConHor;//connections within this layer
    std::vector<size_t> ConUp;//connections to upper layer

};

class layer{
private:
    const size_t Num_Columns;
    std::vector<size_t> ActColumns;//active columns; maybe turn into array of active
                                //columns of the last few time steps
    std::vector<std::vector<bool>> ActLog;//activity log
    //records the activity of all columns over the last few steps
    std::vector<std::vector<synapse>> PotSyn;//list of potential synapses
    float ConThr;//threshhold for synapses in order to be "connected"
    std::vector<std::vector<synapse>> ConSyn;//list of "connected" synapses

    //maybe reconsider connectedness mechanism
    double CondsInc;//connectedness increment
    double CondsDec;//connectedness decrement

    //stuff for inhibition
    size_t DesLocAct;//desired local activity
    //a measure of how many columns should be active on average
    //probably unnecessary consider changing inhibition
    size_t InhiRad;//inhibition radius used to update neighbors
    std::vector<std::vector<size_t>> DistMat;//distance matrix, used to update neighbors
    std::vector<std::vector<size_t>> neighbors;//list of neighbors for each column, consider changing inhibition
    size_t MinOverlap;//minimum overlap
    //end of inhibition stuff

    std::vector<double> boost();//vector of boost values increases
                        //activity of columns which are not active enough




public:
    //reconsider private/public position of memberfunctions!
    double overlap(std::vector<bool> input,size_t column);//maybe returning a vector and omitting second arguent is smarter
    //computes overlap of a particular column with given input

    void updateBoost(std::vector<bool> input);
    //checks if activity log agrees with the wanted average activity

    void updateSyn(std::vector<bool> input);
    //updates all the synapses of all the cells in all the columns
        //all synapses in active columns get reinforced if they were
        //active right before activation of the corresponding column
        //all the nonactive synapses in active columns are decremented
        //---- the synapses which in cells whose overlap with input
        //is lower then proper, are uniformly reinforced as well

};

class brain{
private:
    const size_t NumLevels;
    const std::vector<layer> ListLevels;

};

int main(int argc, char *argv[])
{
    std::vector<int> l{1,2};
    for(auto &k : l){
        std::cout<<k<<"vector entry"<<std::endl;
    }
    std::cout << "Hello World!" << std::endl;
    return 0;
}
