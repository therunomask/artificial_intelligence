#include <iostream>
#include<vector>

class synapse{
    double conds;//connectedness
    std::vector<size_t> ConIn;//connections to inputlayer
};

class level{
private:
    const size_t Num_Columns;
    std::vector<size_t> ActColumns;//active columns
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



public:
    //reconsider private/public position of memberfunctions!
    double overlap(std::vector<bool> input,size_t column);//maybe returning a vector and omitting second arguent is smarter
    //computes overlap of a particular column with given input
    std::vector<double> boost();
};

class brain{
private:
    const size_t NumLevels;
    const std::vector<level> ListLevels;

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
