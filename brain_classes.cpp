#include "brain_header.h"
#include <vector>
#include <deque>
#include <queue>
#include <algorithm>


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
        //reset activity to false
        pillar.active=false;
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
    std::transform(ColumnActivityLog.begin(), ColumnActivityLog.end(), ColumnActivityLog.begin(),
                   std::bind1st(std::multiplies<double>(),AverageExp));
    //multiply running average by AverageExp

    std::transform(ColumnActivityLog.begin(), ColumnActivityLog.end(), activity.begin(),
                        ColumnActivityLog.begin(), std::plus<double>());
    //add new activity to running average.
    //last to lines correspond to "ColumnActivityLog= AverageExp*ColumnActivityLog+activity"

    //more learning; activity based
    std::vector<double>::iterator pActivity_max;
    pActivity_max= std::max_element(ColumnActivityLog.begin(), ColumnActivityLog.end());
    int Max_overlap=0;

    double MaxActivity = (*pActivity_max);
    for(size_t i=0;i<Num_Columns;++i){
        //arbitrary boost function!!
        ColumnList[i].boosting=3.0-2.0*(ColumnActivityLog[i]/MaxActivity);
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
        for(auto& activePillarCell: ColumnList[active_pillar].CellList){
            if(activePillarCell.expect[1]==true){
                segment* s=NULL;
                //find segment of cell that signified the end of a sequence
                for(auto& active_segment: activePillarCell.ActiveSegments[1]){
                    if(active_segment->EndOfSeq==true){
                        s=active_segment;
                        break;
                    }
                }//if no such segment exists; continue with next cell
                if( s==NULL){
                    continue;
                }
                predicted=true;
                activePillarCell.active[0]=true;
                //choose current cell to be the learning cell if it
                //is connected to a cell with learn state on
                for(auto& connected_cells : s->Synapse){
                    if(connected_cells.first->learn[1]==true){
                        activePillarCell.learn[0]=true;
                        is_chosen=true;
                        break;
                    }
                }
            }
        }
        if(predicted==false){
            for(auto& PillarCell:ColumnList[active_pillar].CellList){
                PillarCell.active[0]=true;
            }
        }
        if(is_chosen==false){
            //get best matching cell in last timestep
            segment* BestSegment= ColumnList[active_pillar].BestMatchingSegmentInColumn();
            BestSegment->MotherCell->active[0]=true;
            BestSegment->BlindSynapseAdding(this,1);//1= most recent
            BestSegment->EndOfSeq=true;
            SegmentUpdateList[1].push_back(BestSegment);
        }

    }

    for(auto& pillars:ColumnList){
        for(auto& dummycell:pillars.CellList){
            dummycell.expect.erase(dummycell.expect.begin());
            dummycell.expect.push_back(false);
            //check whether cell currently has an active segment
            if(dummycell.ActiveSegments[1].size()!=0){
                //cell predicts now, because of active segment
                dummycell.expect[1]=true;
                //first element of ActiveSegments[1] is most active segment
                //this segment is supposed to learn
                SegmentUpdateList[1].push_back(dummycell.ActiveSegments[1][0]);

                if(dummycell.expect[0]==true){
                //if the cell predicted prediction, the most active segment
                    //of the previous timestep learns
                    SegmentUpdateList[1].push_back(dummycell.ActiveSegments[0][0]);
                }
                else{
                    //find the Segment that matches the activity of the previous
                    //timestep best. Do blind synapse adding for this segment,
                    //we choose this kind of learning, because the current
                    //prediction was not predicted.
                    segment* poBestSegment=dummycell.BestSegment(0);
                    poBestSegment->BlindSynapseAdding(this,0);
                    SegmentUpdateList[1].push_back(poBestSegment);
                }

            }

        }
    }

    //learning, i.e. implementing changes queued up in SegmentUpdateList
    for(auto& dummySegment: SegmentUpdateList[0]){
        //if the cell is active, the synapse should be reinforced.
        //If the cell is inactive, the synapse should only be reinforced if
        //the segment is not at the end of a sequence, otherwise the segment
        //should be weakened
        if(dummySegment->MotherCell->active[0]==true){
            dummySegment->AdaptingSynapses(dummySegment->PositiveLearning);
        }
        else if(dummySegment->EndOfSeq==false&&dummySegment->MotherCell->expect[0]){
            dummySegment->AdaptingSynapses(dummySegment->PositiveLearning);
        }
        else{
            dummySegment->AdaptingSynapses(!(dummySegment->PositiveLearning));
        }
    }

    std::vector<bool> activation_prediction;
    //set for all cells in all columns active, expect, learn=false
    //also delete last element of cell.active and cell.learn and add
    //new element= false
    //for next timestep!
    return activation_prediction;
}

segment* column::BestMatchingSegmentInColumn(void){
    //find the best matching segment of all the cells in the column
    //return that segment
    size_t max_count=0;
    segment* bestSegment= &CellList[0].SegList[0];
    for(auto& dummy_cell: CellList){
        size_t count=0;
        for(auto& dummy_segment:dummy_cell.SegList){
            for(auto& remote_cell: dummy_segment.Synapse){
                if(remote_cell.first->active[1]==true){++count;}
            }
            if(count> max_count){
                max_count=count;
                bestSegment=&dummy_segment;
            }
        }
    }

    return bestSegment;
}


void segment::BlindSynapseAdding(layer* level,size_t t){
    //add all active synapses to a segment that
    //did not connect to a cell in learnstate

    for(auto& remote_cell:level->Three_CellActivityList[t]){
        bool present=false;
        for(auto& dummysynapse: Synapse){
            if(remote_cell==dummysynapse.first){
                present=true;
                break;
            }

        }
        if(present==false){
            AddCell( remote_cell);
        }
    }

}

void cell::UpdateActiveSegments(void){
    //erase first element and insert list of new active segments at the end
    //determine which segments are active by weighted sum of active cells
    //to wich the segments point

    //when updating brain:!!!!!!!!!!!!!!!!!
    //first update all segments of all cells of all columns of all layers,
    //then update update all cells of all columns of all layers
    //then all columns of all layers
    ActiveSegments.erase(ActiveSegments.begin());
    std::vector<segment*> tempSegments;
    size_t max=0;
    double max_sum=0;
    for(auto& dummysegment:SegList){
        double sum=0;
        for(auto& dummysynapse:dummysegment.Synapse){
            if(dummysynapse.first->active[0]==true){
                sum+=dummysynapse.second;
            }
        }
        if(sum>=dummysegment.MinSynapseWeightActivity){
            tempSegments.push_back(&dummysegment);
            //remember largest element to put in front of the vector
            if(max_sum<sum){
                max=tempSegments.size();
                max_sum=sum;
            }

        }
    }

    std::iter_swap(tempSegments.begin(),tempSegments.begin()+max-1);

    ActiveSegments.push_back(tempSegments);
}

segment* cell::BestSegment(size_t t){
    //finds the Segment which matches best at time t
    //and returns a pointer to it
    double max=0;
    segment* pBestSegment;
    for(auto& dummysegment: SegList){
        double sum=0;
        for(auto& dummysynapse: dummysegment.Synapse){
            if(dummysynapse.first->active[t]==true){
                //add activity of synapse only if remote cell is active
                sum+=dummysynapse.second;
            }
        }
        if(sum>=max){
            max=sum;
            pBestSegment=&dummysegment;
        }
    }
    return pBestSegment;
}
void segment::AdaptingSynapses(bool positive){

    //for positive learning reinforce all connections to active cells
    //punish all connections to inactive cells
    if(positive==true){
        //no for loop, because we delete synapses that drop below 0 connectedness
        auto dummySynapse=Synapse.begin();
        while( dummySynapse!=Synapse.end()){
            if(dummySynapse->first->active[0]==true){
                //if remote cell is active-> positive learning is apropriate
                dummySynapse->second=std::min(1.0,dummySynapse->second+Learnincrement);
                ++dummySynapse;
            }
            else{
                //if remote cell is inactive, the synapse was useless, negative learning is apropriate
                dummySynapse->second=dummySynapse->second-Learnincrement;
                if(dummySynapse->second<=0){
                    Synapse.erase(dummySynapse);
                }
                else{
                    ++dummySynapse;
                }
            }
        }
    }
    //for negative learning punish all connections to active cells,
    //since they predicted mistakenly
    else{
        auto dummySynapse=Synapse.begin();
        while( dummySynapse!=Synapse.end()){
            if(dummySynapse->first->active[0]==true){
                //if remote cell is active-> decrement connectedness
                dummySynapse->second=dummySynapse->second-Learnincrement;
                if(dummySynapse->second<=0){
                    Synapse.erase(dummySynapse);
                }
                else{
                    ++dummySynapse;
                }
            }
            else{
                ++dummySynapse;
            }

        }
    }
}

