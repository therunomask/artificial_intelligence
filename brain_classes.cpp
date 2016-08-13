#include "brain_header.h"
#include <vector>
#include <deque>

void segment::UpdateCon(std::vector<const cell*> winners ){//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
}


std::vector<bool> layer::current_activation(){
    std::vector<bool>input( p_lower_level->current_prediction());
        //trigger prediction evaluation of lower levels while optaining input
        //
        //missing: propagate expectation to lower level!!
        //
    std::vector<bool> activation;//reserve space for length of ColumnList
    for(auto &pillar : ColumnList){
        //compute activation of columns in the layer
        activation.push_back(pillar.feed_input(input));
    }
}
