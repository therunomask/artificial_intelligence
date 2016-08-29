#include "brain_header.h"
#include <vector>
#include <deque>
#include <queue>


void segment::UpdateCon(std::vector<const cell*> winners ){//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
}


  std::vector<bool> layer::current_activation(void){
    std::vector<bool>input( p_lower_level->current_prediction());
        //trigger prediction evaluation of lower levels while optaining input
        //
        //missing: propagate expectation to lower level!!
        //
    std::vector<double> overlap;//reserve space for length of ColumnList
    for(auto &pillar : ColumnList){
        //compute feed input of columns in the layer
        overlap.push_back(pillar.feed_input(input));
    }
    //finding #DesiredLocalActivity highest overlapping columns

    std::priority_queue<std::pair<double, int>,std::vector<std::pair<double,int>>,std::greater<std::pair<double,int>>> winner;
    //priority queue orderes its elements automatically
    //long definition to make priority queue order its elements in increasing order
    for (int i = 0; i < DesiredLocalActivity; ++i) {
        winner.push(std::pair<double, int>(overlap[i], i));
    }//first add minimum number of elements
    for (int i = DesiredLocalActivity; i < Num_Columns; ++i) {
        if(overlap[i]>winner.top().first){
              winner.pop();
              winner.push(std::pair<double, int>(overlap[i], i));
        }
    }//winner now contains #DesiredLocalActivity highest overlapping columns
    //now create activity vector
    std::vector<bool> activity(Num_Columns,0);
    for (int i = 0; i < DesiredLocalActivity; ++i) {
        int index = winner.top().second;
        ActColumns[i]=index;
        activity[index]=true;
        winner.pop();
    }
    //change activity log
    return activity;
}

double column::feed_input(std::vector<bool> input){//computes overlap with input
    int overlap=0;
    for(auto &syn : ConnedSynapses){
            overlap+= input[*syn];
    }
    if(overlap<MinOverlap){
        overlap=0;
    }
    return overlap*boosting;
}



