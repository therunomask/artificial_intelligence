#include "brain_header.h"
#include <vector>
#include <deque>
#include <queue>


void segment::UpdateCon(std::vector<const cell*> winners ){//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
}


std::vector<bool> layer::activation_learning(void){
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

    std::priority_queue<std::pair<double, size_t>,std::vector<std::pair<double,size_t>>,std::greater<std::pair<double,size_t>>> winner;
    //priority queue orderes its elements automatically
    //long definition to make priority queue order its elements in increasing order
    for (size_t i = 0; i < DesiredLocalActivity; ++i) {
        winner.push(std::pair<double, size_t>(overlap[i], i));
    }//first add minimum number of elements
    for (size_t i = DesiredLocalActivity; i < Num_Columns; ++i) {
        if(overlap[i]>winner.top().first){
              winner.pop();
              winner.push(std::pair<double, size_t>(overlap[i], i));
        }
    }//winner now contains #DesiredLocalActivity highest overlapping columns
    //now create activity vector
    std::vector<bool> activity(Num_Columns,false);
    for (size_t i = 0; i < DesiredLocalActivity; ++i) {
        size_t index = winner.top().second;
        ActColumns[i]=index;
        activity[index]=true;
        winner.pop();
        //learning
        for(size_t l=0;l<input.size();++l){
            if(input[l]){
                ColumnList[index].connectedness[l]+=CondsInc;
                ColumnList[index].connectedness[l]=(ColumnList[index].connectedness[l]>1) ? 1 : ColumnList[index].connectedness[l];
            }else{
                ColumnList[index].connectedness[l]-=CondsInc;
                ColumnList[index].connectedness[l]=(ColumnList[index].connectedness[l]<0) ? 0 : ColumnList[index].connectedness[l];
            }
        }
        //end learning
    }
    //change activity log
    std::transform(ActivityLog.begin(), ActivityLog.end(), ActivityLog.begin(),
                   std::bind1st(std::multiplies<double>(),AverageExp));
    //multiply running average by AverageExp

    std::transform(ActivityLog.begin(), ActivityLog.end(), activity.begin(),
                        ActivityLog.begin(), std::plus<double>());
    //add new activity to running average.
    //last to lines correspond to "ActivityLog= AverageExp*ActivityLog+activity"

    //more learning; activity based
    std::vector<double>::iterator pActivity_max;
    pActivity_max= std::max_element(ActivityLog.begin(), ActivityLog.end());
    int Max_overlap=0;

    double MaxActivity = (*pActivity_max);
    for(size_t i=0;i<Num_Columns;++i){
        //arbitrary boost function!!
        ColumnList[i].boosting=3.0-2.0*(ActivityLog[i]/MaxActivity);
        //optimize w.r.t. this function!!

        if(ColumnList[i].tell_overlap_average()>Max_overlap){
            Max_overlap=ColumnList[i].tell_overlap_average();
        }
    }
    Max_overlap*=AverageOverlapMin;//otimize magic number!

    for(auto& pillar : ColumnList){
        //strenthen pillars that never overlap uniformly
        if(pillar.tell_overlap_average() < Max_overlap){
            for(auto& con: pillar.connectedness){
                con+=SpecialOverlapIncrement;
            }
        }

    }//end of learning

    return activity;
}

double column::feed_input(const std::vector<bool>& input){
    //computes overlap and
    //boosted overlap; also updates running averages
    int overlap=0;
    for(auto &syn : ConnectedSynapses){
            overlap+= input[syn];
    }
    if(overlap<MinOverlap){
        overlap=0;
    }
    Overlap_Average= Average_Exp*Overlap_Average+ overlap;
    return overlap*boosting;
}

std::vector<bool> layer::current_prediction( void ){
        //triggers activation of same level

    //check if a cell predicted activation of a column
    //if so, change its synapses.
    //if not choose cell to predict same activation in the future
    for(auto& active_pillar: ActColumns){
        bool predicted=false;//dummy variable checks of predicting cell is found
        bool is_chosen=false;//also dummy
        for(auto& cell: ColumnList[active_pillar].CellList){
            if(cell.expect==true){
                segment* s=NULL;
                //find segment of cell that signified the end of a sequence
                for(auto& active_segment: cell.ActiveSegments[1]){
                    if(active_segment->EndOfSeq==true){
                        s=active_segment;
                        break;
                    }
                }//if no such segment exists; continue with next cell
                if( s==NULL){
                    continue;
                }
                predicted=true;
                cell.active=true;
                //choose current cell to be the learning cell if it
                //is connected to a cell with learn state on
                for(auto& connected_cells:s->CellAddr){
                    if(connected_cells->learn[1]==true){
                        cell.learn[0]=true;
                        is_chosen=true;
                        break;
                    }
                }
            }
        }
    }

    std::vector<bool> activation_prediction;
    //set for all cells in all columns active, expect, learn=false
    //for next timestep!
    return activation_prediction;
}


