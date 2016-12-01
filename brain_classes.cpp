#include "brain_header.h"
#include <vector>
#include <deque>
#include <queue>
#include <algorithm>
#include <time.h>


std::vector<layer> BrainConstructionHelper(brain Init_brain,size_t Number_of_Levels, layer& LowestLayer, layer& Top, size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column){

    //maybe not yet finished, check this!

    std::vector<layer> ListOfLayers;
    ListOfLayers.push_back(LowestLayer);

    for(size_t i=1;i<Number_of_Levels-2;++i){
        ListOfLayers.push_back(*(new layer(Number_of_Column_per_Layer, Number_of_Cells_per_Column)));

    }
    ListOfLayers.push_back(Top);

    for(size_t i=1;i<ListOfLayers.size()-1;++i){
// initialize lowest and highest level extra
        ListOfLayers[i].p_lower_level=&(ListOfLayers[i-1]);
        ListOfLayers[i].p_upper_level=&(ListOfLayers[i+1]);

    }

    //now that layers exist, we connect columns to other columns
    //we choose at random which columns to connect to and
    //how strong the connection is
    srand(time(NULL));
    for(layer& DummyLayer: ListOfLayers){
        if(DummyLayer.p_lower_level==NULL){
            continue;
        }
        for(column& DummyPillar: DummyLayer.ColumnList){
            for(size_t SynapseIndex=0;SynapseIndex<DummyLayer.Num_Columns;SynapseIndex++){
                //                                                                                                                                                                                          uniformly distributed between MinoverLap(0.3) and 2 Minoverlap(0.6)
                DummyPillar.ConnectedSynapses.push_back(std::pair<column*,double>(&(DummyLayer.p_lower_level->ColumnList[rand()%DummyLayer.p_lower_level->Num_Columns]),(static_cast<double>(rand())/RAND_MAX+1)*minimal_overlap_under_consideration/(pillars_per_layer*active_pillers_per_pillar)));
            }

        }
    }

    //now that layers exist, we connect segments to cells in their own +-1 level
    //we choose at random which cell to connect to and
    //how strong the connection is
    for(layer& DummyLayer: ListOfLayers){
        if(DummyLayer.p_lower_level==NULL){
            continue;
        }
        for(column& DummyPillar: DummyLayer.ColumnList){
            for(cell& DummyCell: DummyPillar.CellList){
                for(segment& DummySegment: DummyCell.SegList){
                    for(size_t SynapseIndex=0;SynapseIndex<synapses_per_segment;SynapseIndex++){
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer.ColumnList[rand()%DummyLayer.Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer.p_lower_level->ColumnList[rand()%DummyLayer.Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer.p_upper_level->ColumnList[rand()%DummyLayer.Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                    }

                }
            }

        }
    }

    return ListOfLayers;
}

brain::brain(size_t Number_of_Levels, layer& LowestLayer, layer& Top, size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column)
    :NumLevels(Number_of_Levels),
      ListLevels(BrainConstructionHelper(*this, NumLevels ,LowestLayer, Top,Number_of_Column_per_Layer,Number_of_Cells_per_Column))

{


}



layer::layer(size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column)
    :

      DesiredLocalActivity((static_cast<int>(Number_of_Column_per_Layer*FractionOfActiveColumns))),
      ActColumns(std::vector<column*>(DesiredLocalActivity,NULL)),
      TempActColumns(std::vector<column*>(DesiredLocalActivity,NULL)),
      Num_Columns(Number_of_Column_per_Layer),
      p_lower_level(NULL),//still to be initialized for bottommost and highes layer
      p_upper_level(NULL),//still to be initialized for bottommost and highes layer
      ColumnList(std::vector<column>(Number_of_Column_per_Layer,column(this,Number_of_Cells_per_Column))),
      SegmentUpdateList(std::vector<std::vector<segment*>>(2,std::vector<segment*>())),
      CellActivityList(std::vector<cell*>()),
      Three_CellActivityList(std::vector<std::vector<cell*>>(2,std::vector<cell*>()))
{

}
column::column(layer* layer_to_belong_to, size_t Number_of_Cells_per_Column)
    :
    Overlap_Average(Initial_Overlap_Average),
    ConnectedSynapses(std::vector<std::pair<column*,double>>()),//still to be initialized for lowest and highest layer
    MotherLayer(layer_to_belong_to),
    CellList(std::vector<cell>(Number_of_Cells_per_Column,cell(this))),
    active(false),
    ActivityLog(Initial_Activity_log),
    boosting(1.0)
{

}


cell::cell(column* Column_to_belong_to)
    :
      MotherColumn(Column_to_belong_to),
      SegList(std::vector<segment>(synapses_per_segment, segment(this))),
      ActiveSegments( std::vector<std::vector<segment*>>(2,*(new std::vector<segment*>))),
      active(2,false),
      expect(2,false),
      learn(2,false)
{

}

segment::segment(cell* Cell_to_belong_to)
    :
      active(false),
      MotherCell(Cell_to_belong_to),
      Synapse(std::vector< std::pair <cell*,double>>() ),//still to be initialized afterward
      EndOfSeq(false)
{

}


void segment::UpdateCon(std::vector<const cell*> winners ){//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
}

void layer::FindBestColumns(void){
    //finding #DesiredLocalActivity highest overlapping columns

    std::priority_queue<std::pair<double, column*>,std::vector<std::pair<double,column*>>,std::greater<std::pair<double,column*>>> winner;
    //priority queue orderes its elements automatically
    //long definition to make priority queue order its elements in increasing order
    for (size_t i = 0; i < DesiredLocalActivity; ++i) {
        winner.push(std::pair<double, column*>(ColumnList[i].feed_input(),&ColumnList[i] ));
    }//first add minimum number of elements
    for (size_t i = DesiredLocalActivity; i < Num_Columns; ++i) {
        double overlap=ColumnList[i].feed_input();
        if(overlap >winner.top().first){
              winner.pop();
              winner.push(std::pair<double, column*>(overlap, &ColumnList[i]));
        }
    }//winner now contains #DesiredLocalActivity highest overlapping columns


    for (size_t i = 0; i < DesiredLocalActivity; ++i) {
        TempActColumns[i]=winner.top().second;
    }

}

void layer::ConnectedSynapsesUpdate(void){
    //strengthen connections to columns that were active
    //weaken connections to columns that were inactive
    for(auto& pdummyColumn: TempActColumns){

        for(auto& dummy_connected_synapse : pdummyColumn->ConnectedSynapses){
            if(dummy_connected_synapse.first->active==true){
                dummy_connected_synapse.second=std::max(dummy_connected_synapse.second+CondsInc,1.0);
            }else{
                dummy_connected_synapse.second=std::min(dummy_connected_synapse.second-CondsDec,0.0);

            }
        }
    }

}
double layer::ActivityLogUpdateFindMaxActivity(void){
    double MaxActivity = 0;

    for(auto& dummy_pillar:ColumnList){
        //change activity log
        dummy_pillar.ActivityLog =dummy_pillar.ActivityLog* dummy_pillar.Average_Exp+dummy_pillar.active;

        //more learning; activity based
        if(dummy_pillar.ActivityLog>MaxActivity){
            MaxActivity=dummy_pillar.ActivityLog;
        }
    }

    return MaxActivity;

}


void layer::BoostingUpdate_StrenthenWeak(double MaxActivity){

    int Max_overlap=0;

    for(auto& dummy_column:ColumnList){
        //arbitrary boost function!!
        dummy_column.boosting=maximum_boosting-(maximum_boosting-1.0)*(dummy_column.ActivityLog/MaxActivity);
        //optimize w.r.t. this function!!
        if(dummy_column.tell_overlap_average()>Max_overlap){
            Max_overlap=dummy_column.tell_overlap_average();
        }
    }
    Max_overlap *=AverageOverlapMin;//otimize magic number!

    for(auto& pillar : ColumnList){
        //strenthen pillars that never overlap uniformly
        if(pillar.tell_overlap_average() < Max_overlap){
            for(auto& con: pillar.ConnectedSynapses){
                con.second+=SpecialOverlapIncrement;
            }
        }

    }//end of learning

}



double column::feed_input(void){
    //computes overlap and
    //boosted overlap; also updates running averages
    int overlap=0;
    for(auto &syn : ConnectedSynapses){
            overlap+=static_cast<int>( syn.first->active);
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
        for(auto& activePillarCell: active_pillar->CellList){
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
            for(auto& PillarCell:active_pillar->CellList){
                PillarCell.active[0]=true;
            }
        }
        if(is_chosen==false){
            //get best matching cell in last timestep
            segment* BestSegment= active_pillar->BestMatchingSegmentInColumn();
            BestSegment->MotherCell->learn[0]=true;
            BestSegment->BlindSynapseAdding(this,1);//1= most recent
            BestSegment->EndOfSeq=true;
            SegmentUpdateList[1].push_back(BestSegment);
        }

    }
        //check 0==now! <-> change if not
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
                    segment* pBestSegment=dummycell.BestSegment(0);
                    pBestSegment->BlindSynapseAdding(this,0);
                    SegmentUpdateList[1].push_back(pBestSegment);
                }

            }

        }
    }

    //learning, i.e. implementing changes queued up in SegmentUpdateList
    for(auto& dummySegment: SegmentUpdateList[0]){

        //change to waiting for actual activation

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
                dummySynapse->second=std::min(1.0,dummySynapse->second+LearnIncrement);
                ++dummySynapse;
            }
            else{
                //if remote cell is inactive, the synapse was useless, negative learning is apropriate
                dummySynapse->second=dummySynapse->second-LearnIncrement;
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
                dummySynapse->second=dummySynapse->second-LearnIncrement;
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

void brain::update(){
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
 *
 * Do all layers in parallel, i.e. do the layers one after another
 * but don't implement updates until after the last layer.

*/

    for(layer Dummylayer:ListLevels){
        Dummylayer.FindBestColumns();
        Dummylayer.ConnectedSynapsesUpdate();
        double MaxActivity=Dummylayer.ActivityLogUpdateFindMaxActivity();
        Dummylayer.BoostingUpdate_StrenthenWeak( MaxActivity);

    }

}

