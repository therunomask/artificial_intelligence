#include "brain.h"
#include <iostream>
#include <vector>
#include <deque>
#include <queue>
#include <algorithm>
#include <time.h>
#include <unordered_set>
#include <typeinfo>
#include <stdexcept>
#include <thread>


void BrainConstructionHelper(brain& Init_brain,size_t Number_of_Levels, size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column/*,std::vector<bool>(*sensoryinput)(size_t time)*/){

    //maybe not yet finished, check this!

    //std::vector<layer> LayerList;
    //LayerList.push_back(new bottom_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column,Init_brain,sensoryinput));

    //for(size_t i=0;i<Number_of_Levels-2;++i){
    //    LayerList.emplace_back(Number_of_Column_per_Layer, Number_of_Cells_per_Column, Init_brain);
    //}
    //LayerList.push_back(new top_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column, Init_brain));

    Init_brain.LowestLayer.p_upper_level=&(Init_brain.ListOfLevels[0]);
    Init_brain.ListOfLevels[0].p_lower_level=&(Init_brain.LowestLayer);
    Init_brain.ListOfLevels[0].p_upper_level=&(Init_brain.ListOfLevels[1]);
    for(size_t i=1;i<Init_brain.ListOfLevels.size()-1;++i){
    // initialize lowest and highest level extra
        Init_brain.ListOfLevels[i].p_lower_level=&(Init_brain.ListOfLevels[i-1]);
        Init_brain.ListOfLevels[i].p_upper_level=&(Init_brain.ListOfLevels[i+1]);

    }
    Init_brain.ListOfLevels[Init_brain.ListOfLevels.size()-1].p_lower_level=&(Init_brain.ListOfLevels[Init_brain.ListOfLevels.size()-2]);
    Init_brain.ListOfLevels[Init_brain.ListOfLevels.size()-1].p_upper_level=&(Init_brain.HighestLayer);
    Init_brain.HighestLayer.p_lower_level=&(Init_brain.ListOfLevels[Init_brain.ListOfLevels.size()-1]);

    //now that layers exist, we connect columns to other columns
    //we choose at random which columns to connect to and
    //how strong the connection is

    srand(282346120965);//time(NULL));
    //srand(time(NULL));
    for(layer*& DummyLayer: Init_brain.AllLevels){
        if(DummyLayer->p_lower_level==NULL){
            continue;
        }
        for(column& DummyPillar: DummyLayer->ColumnList){
            for(size_t SynapseIndex=0;SynapseIndex<DummyLayer->Num_Columns;SynapseIndex++){
                //                                                                                                                                                                                          uniformly distributed between MinoverLap/2(0.15) and 1.5 Minoverlap(0.45)
                //DummyPillar.ConnectedSynapses.push_back(std::pair<column*,double>(&(DummyLayer.p_lower_level->ColumnList[rand()%DummyLayer.p_lower_level->Num_Columns]),(static_cast<double>(rand())/RAND_MAX+1)*minimal_overlap_under_consideration/(pillars_per_layer*active_pillers_per_pillar)));
                DummyPillar.ConnectedSynapses.push_back(std::pair<column*,double>(&(DummyLayer->p_lower_level->ColumnList[SynapseIndex]),(static_cast<double>(rand())/RAND_MAX+0.5)*minimal_overlap_under_consideration/(pillars_per_layer*active_pillers_per_pillar)));
            }

        }
    }

    //now that layers exist, we connect segments to cells in their own +-1 level
    //we choose at random which cell to connect to and
    //how strong the connection is

    for(layer*& DummyLayer: Init_brain.AllLevels){
        if(DummyLayer->p_lower_level==NULL){
            for(column& DummyPillar: DummyLayer->ColumnList){
                for(cell& DummyCell: DummyPillar.CellList){
                    for(segment& DummySegment: DummyCell.SegList){
                        for(size_t SynapseIndex=0;SynapseIndex<synapses_per_segment;SynapseIndex++){
                            DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                            DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->p_upper_level->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        }

                    }
                }

            }
            continue;
        }
        if(DummyLayer->p_upper_level==NULL){
            for(column& DummyPillar: DummyLayer->ColumnList){
                for(cell& DummyCell: DummyPillar.CellList){
                    for(segment& DummySegment: DummyCell.SegList){
                        for(size_t SynapseIndex=0;SynapseIndex<synapses_per_segment;SynapseIndex++){
                            DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                            DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->p_lower_level->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        }

                    }
                }
            }
            continue;
        }
        for(column& DummyPillar: DummyLayer->ColumnList){
            for(cell& DummyCell: DummyPillar.CellList){
                for(segment& DummySegment: DummyCell.SegList){
                    for(size_t SynapseIndex=0;SynapseIndex<synapses_per_segment;SynapseIndex++){
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->p_lower_level->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                        DummySegment.Synapse.push_back(std::pair<cell*,double>(&(DummyLayer->p_upper_level->ColumnList[rand()%DummyLayer->Num_Columns].CellList[rand()%Number_of_Cells_per_Column]),(static_cast<double>(rand())/RAND_MAX+1)*Minimal_sum_of_synapseweights_for_activity/(cells_per_column*active_pillers_per_pillar*pillars_per_layer)));
                    }

                }
            }

        }
    }


}

brain::brain(size_t Number_of_Levels, size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column,std::vector<bool>(*sensoryinput)(size_t time))
    :NumLevels(Number_of_Levels),
      ListOfLevels(std::vector<layer>()),
      LowestLayer(bottom_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column,*this,sensoryinput)),
      HighestLayer(top_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column,*this)),
      AllLevels(std::vector<layer*>()),
      time(0),
      Martin_Luther(debughelper()),
      max_activation_counter(5),
      max_activation_counter_change(false),
      max_activation_counter_mutex()
{
    ListOfLevels.reserve(Number_of_Levels);
    for(size_t i=0;i<Number_of_Levels-2;++i){
        ListOfLevels.emplace_back(Number_of_Column_per_Layer, Number_of_Cells_per_Column, *this);
    }
    AllLevels.push_back(&LowestLayer);
    for(layer& dummylayer:ListOfLevels){
        AllLevels.push_back(&dummylayer);
    }
    AllLevels.push_back(&HighestLayer);
    BrainConstructionHelper(*this, NumLevels ,Number_of_Column_per_Layer,Number_of_Cells_per_Column/*,sensoryinput*/);


}



layer::layer(size_t Number_of_Column_per_Layer, size_t Number_of_Cells_per_Column, brain& pBrain)
    :
      DesiredLocalActivity((static_cast<int>(Number_of_Column_per_Layer*FractionOfActiveColumns))),
      ActColumns(std::vector<column*>(DesiredLocalActivity,NULL)),
      TempActColumns(std::vector<column*>(DesiredLocalActivity,NULL)),
      MotherBrain(pBrain),
      Num_Columns(Number_of_Column_per_Layer),
      p_lower_level(NULL),//still to be initialized for bottommost and highes layer
      p_upper_level(NULL),//still to be initialized for bottommost and highes layer
      ColumnList(std::vector<column>()),
      CellActivityList(std::vector<cell*>()),
      Three_CellActivityList(std::vector<std::vector<cell*>>(2,std::vector<cell*>())),
      CellUpdateList(std::vector<cell*>()),
      PendingActivity(std::vector<cell*>()),
      PendingExpectation(std::vector<cell*>())
{
    ColumnList.reserve(Number_of_Column_per_Layer);
    for(size_t index=0;index<Number_of_Column_per_Layer;++index){
        ColumnList.emplace_back(*this,Number_of_Cells_per_Column);
    }
    srand(time(NULL));
    for (column*& dummycolumn: ActColumns){
        dummycolumn=&ColumnList[rand()%ColumnList.size()];
    }

}

layer::layer(const layer &dummylayer):
    MotherBrain(dummylayer.MotherBrain),
    Num_Columns(dummylayer.Num_Columns)

{
    std::cout<<"Freeze! You are copying layers! \n";
    throw std::invalid_argument("Don't copy layers \n");

   // return bla;

}

top_layer::top_layer(size_t Number_of_Column_per_layer, size_t Number_of_Cells_per_Column, brain &pBrain):
    layer( Number_of_Column_per_layer, Number_of_Cells_per_Column, pBrain)
{
}

bottom_layer::bottom_layer(size_t Number_of_Column_per_layer, size_t Number_of_Cells_per_Column, brain &pBrain, std::vector<bool>(*sensoryinput)(size_t time))
    :
    layer( Number_of_Column_per_layer, Number_of_Cells_per_Column, pBrain)
{
    external_input=sensoryinput;
}




column::column(layer& layer_to_belong_to, size_t Number_of_Cells_per_Column)
    :
    Overlap_Average(Initial_Overlap_Average),
    ConnectedSynapses(std::vector<std::pair<column*,double>>()),//still to be initialized for lowest and highest layer
    MotherLayer(layer_to_belong_to),
    CellList(std::vector<cell>()),
    active(false),
    expect(false),
    ActivityLog(Initial_Activity_log),
    boosting(1.0)
{
    CellList.reserve(Number_of_Cells_per_Column);
    for(size_t index=0; index<Number_of_Cells_per_Column;++index){
        CellList.emplace_back(*this);
    }

}

column::column(const column &dummycolumn):
    MotherLayer(dummycolumn.MotherLayer)
{
    throw std::invalid_argument("Don't copy columns! \n");
}


cell::cell(column& Column_to_belong_to)
    :
      MotherColumn(Column_to_belong_to),
      SegList(std::deque<segment>()),
      ActiveSegments( std::vector<std::vector<segment*>>(2,*(new std::vector<segment*>))),
      active(2,false),
      expect(2,false),
      SegmentUpdateList(std::vector<SegmentUpdate>())
{
    //SegList.reserve(synapses_per_segment);
    for(size_t index=0; index<synapses_per_segment;++index){
        SegList.emplace_back(*this,1);
    }

}

cell::cell(const cell& dummycell):
    MotherColumn(dummycell.MotherColumn)
{
    throw std::invalid_argument("Don't copy cells! \n");
}

SegmentUpdate::SegmentUpdate(segment*SegmentAddress, std::vector<std::pair<cell *, double> *> active_cells):
    SegmentAddress(SegmentAddress),
    active_cells(active_cells),
    timer(SegmentAddress->ActivationCountdown)
{

}

segment::segment(cell& Cell_to_belong_to, size_t TempActivationCountdown)
    :
      MotherCell(Cell_to_belong_to),
      Synapse(std::deque< std::pair <cell*,double>>() ),//still to be initialized afterward
      ActivationCountdown(TempActivationCountdown)
{

}
segment::segment(const segment &dummysegment):
    MotherCell(dummysegment.MotherCell),
    ActivationCountdown(dummysegment.ActivationCountdown)
{
    throw std::invalid_argument("Don't copy segments! \n");

}

segment segment::operator =(const segment& dummysegment){
    throw std::invalid_argument("Don't operator= segments! \n");
}

std::vector<std::pair<cell*,double>*>  segment::GetActiveCells(){
    std::vector<std::pair<cell*,double>*> tempactive_cells;
    for(std::pair <cell*,double>& dummySynapse: Synapse){
        if(dummySynapse.first->active[1]==true){
            tempactive_cells.emplace_back(&dummySynapse);
        }
    }
    return tempactive_cells;
}


void layer::FindBestColumns(void){
    static size_t counter=0;
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
        winner.pop();
    }


}

void layer::ActiveColumnUpdater(void){
    //implement saved changes

    for(column& dummycolumn:ColumnList){
        dummycolumn.active=false;
    }
    for(size_t i=0;i<TempActColumns.size();++i){
        ActColumns[i]=TempActColumns[i];
        ActColumns[i]->active=true;


    }
    //debughelper statistics
    MotherBrain.Martin_Luther.activation_column[finding_oneself()].push_back(TempActColumns.size());


}


void layer::ConnectedSynapsesUpdate(void){
    //strengthen connections to columns that were active
    //weaken connections to columns that were inactive
    size_t success=0;
    size_t failure=0;

    for(auto& pdummyColumn: ActColumns){

        for(auto& dummy_connected_synapse : pdummyColumn->ConnectedSynapses){
            if(dummy_connected_synapse.first->active==true){
                dummy_connected_synapse.second=std::max(dummy_connected_synapse.second+CondsInc/(pillars_per_layer*active_pillers_per_pillar),1.0);
                //
                ++success;
            }else{
                dummy_connected_synapse.second=std::min(dummy_connected_synapse.second-CondsDec/(pillars_per_layer-pillars_per_layer*active_pillers_per_pillar),0.0);
                //
                ++failure;
            }
        }
    }
    MotherBrain.Martin_Luther.success_column[finding_oneself()].push_back(static_cast<double>(success)/static_cast<double>(success+failure));

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
    //1. too little activity
    for(auto& dummy_column:ColumnList){
        //arbitrary boost function!!
        dummy_column.boosting=maximum_boosting-(maximum_boosting-1.0)*(dummy_column.ActivityLog/MaxActivity);
        //optimize w.r.t. this function!!
        if(dummy_column.tell_overlap_average()>Max_overlap){
            Max_overlap=dummy_column.tell_overlap_average();
        }
    }
    Max_overlap *=AverageOverlapMin;//otimize magic number!
    //2. too little overlap
    for(auto& pillar : ColumnList){
        //strenthen pillars that never overlap uniformly
        if(pillar.tell_overlap_average() < Max_overlap){
            for(auto& con: pillar.ConnectedSynapses){
                con.second+=SpecialOverlapIncrement;
            }
        }

    }//end of learning

}


void layer::Do_SegmentUpdate(){
    //update timers
    //if timer drops to zero or below implement pending update
    //also simultaneously delete used SegmentUpdate objects from the vector
    //and delete element of CellUpdateList if there are no more pending updates for this cell
    for(std::vector<cell*>::iterator itdummycell=CellUpdateList.begin();itdummycell!=CellUpdateList.end();){
        for(std::vector<SegmentUpdate>::iterator itDummyUpdate=(*itdummycell)->SegmentUpdateList.begin();itDummyUpdate!=(*itdummycell)->SegmentUpdateList.end();){
            std::cout<<"pointer position is "<<std::distance((*itdummycell)->SegmentUpdateList.begin(), itDummyUpdate)<<", length of vector is "<<(*itdummycell)->SegmentUpdateList.size()<<std::endl;
            if(itDummyUpdate->timer>0){
                --itDummyUpdate->timer;
                ++itDummyUpdate;
            }else{
                itDummyUpdate->AdaptingSynapses((*itdummycell)->active[0]);
                (*itdummycell)->SegmentUpdateList.erase(itDummyUpdate);
std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

            }
        }
        if((*itdummycell)->SegmentUpdateList.size()==0){
            std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

            CellUpdateList.erase(itdummycell);
        }else{
            ++itdummycell;
        }
    }
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

}

void SegmentUpdate::AdaptingSynapses(bool success){
    //depending on success either we punish all the connections to (previously) active cells
    //or strengthen them.
    //if connection strength drops to zero or below we remove the corresponding synapse from the
    //segment.
    for(std::pair<cell*,double>*& dummySynapse: active_cells){
        if(success==true){
            dummySynapse->second= std::max(static_cast<double>(1),dummySynapse->second+SegmentAddress->LearnIncrement);
            if(SegmentAddress->Synapse.size()<synapses_per_segment){
                SegmentAddress->BlindSynapseAdding(SegmentAddress->ActivationCountdown);
            }
        }else{
            dummySynapse->second= dummySynapse->second- SegmentAddress->LearnIncrement;
            if(dummySynapse->second<=0){
                for(size_t remoteSynapse=0;remoteSynapse<SegmentAddress->Synapse.size();++remoteSynapse ){
                    if(dummySynapse->first==SegmentAddress->Synapse[remoteSynapse].first){
                        SegmentAddress->Synapse.erase(SegmentAddress->Synapse.begin()+remoteSynapse);
                        break;
                    }
                }
            }
        }
    }
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

}





void layer::CellExpectInitiator( void ){


    for(column& pillars:ColumnList){
        for(cell& dummycell:pillars.CellList){
            //check whether cell currently has an active segment

            if(dummycell.ActiveSegments[0].size()!=0){
                //cell predicts now, because of active segment
                PendingExpectation.push_back(&dummycell);
                //first element of ActiveSegments[1] is most active segment
                //this segment is supposed to learn
                CellUpdateList.push_back(&dummycell);
                std::vector<std::pair<cell*,double>*> tempactive_cells=dummycell.ActiveSegments[0][0]->GetActiveCells();
                dummycell.SegmentUpdateList.emplace_back(dummycell.ActiveSegments[0][0],tempactive_cells);

                if(dummycell.expect[1]==false){
                    //if prediction is unexpected,
                    //create a new Segment that matches the activity of the previous
                    //timestep. Do blind synapse adding for this segment,
                    //we choose this kind of learning, because the current
                    //prediction was not predicted.
                    if(MotherBrain.max_activation_counter==dummycell.ActiveSegments[0][0]->ActivationCountdown){
                        MotherBrain.max_activation_counter_change=true;
                    }
                    dummycell.SegList.emplace_back(dummycell,dummycell.ActiveSegments[0][0]->ActivationCountdown+1);
                    dummycell.SegList[dummycell.SegList.size()-1].BlindSynapseAdding(1);//1 = last timestep
                    dummycell.SegmentUpdateList.emplace_back(&dummycell.SegList[dummycell.SegList.size()-1],tempactive_cells);

                }
            }

        }
    }

}



void layer::CellUpdater(void){
    //debughelper statistics
    MotherBrain.Martin_Luther.activation_cell[finding_oneself()].push_back(static_cast<double>(PendingActivity.size())/(cells_per_column*pillars_per_layer));


    CellActivityList=std::vector<cell*>();

    for(column& dummyColumn: ColumnList){
        dummyColumn.expect=false;
        for(cell& dummyCell: dummyColumn.CellList){
        //shift time by one step
            if(dummyCell.active.size()<MotherBrain.max_activation_counter){
                dummyCell.active.push_back(false);
                dummyCell.expect.push_back(false);
            }
            for(size_t index=dummyCell.active.size()-1;index>0;--index){
                dummyCell.active[index]=dummyCell.active[index-1];
                dummyCell.expect[index]=dummyCell.expect[index-1];
            }
            dummyCell.active[0]=false;
            dummyCell.expect[0]=false;
        }
    }
  //set those cells to active,expect,learn that were identified to do so previously
    for(cell*& pdummyCell:PendingActivity){
        CellActivityList.push_back(pdummyCell);
        pdummyCell->active[0]=true;


    }
    PendingActivity=std::vector<cell*>();
    for(cell*& pdummyCell:PendingExpectation){
        pdummyCell->expect[0]=true;
        pdummyCell->MotherColumn.expect=true;
    }
    PendingExpectation=std::vector<cell*>();



}

void layer::CellLearnInitiator(void){
    //check if the activity of a column was predicted;
    //if so, activate the responsible cell and reward
    //the responsible segment
    //if not, activate all the cells in the column
    //and add a segment pointing to the active
    //cells of last round
    for(column*& DummyColumn:ActColumns){
        if(DummyColumn->expect==true){
            //this function finds the segment that fits best the
            //activity of the last timestep with SequenceCounter==1
            //and also leads all the expecting Cells to activation
            segment* BestSegment=DummyColumn->BestMatchingSegmentInColumnActivateCells();
            BestSegment->MotherCell.SegmentUpdateList.emplace_back(BestSegment,BestSegment->GetActiveCells());
        }
        else{
            //activity came unexpected, so we
            //activate all the cells in the column.
            //find the cell with fewest segments
            //and add a segment to this cell, which
            //is connected to active cells of the last round
            size_t SegmentCounter=0;
            size_t SegmentMinCounter=-1;//maximal value of size_t
            cell* PoorestCell=&DummyColumn->CellList[0];
            for(cell& DummyCell:DummyColumn->CellList){
                PendingActivity.push_back(&DummyCell);
                SegmentCounter=DummyCell.SegList.size();
                if(SegmentCounter<SegmentMinCounter){
                    PoorestCell=&DummyCell;
                    SegmentMinCounter=SegmentCounter;
                }
            }
            PoorestCell->SegList.emplace_back(*PoorestCell,1);
            PoorestCell->SegList[PoorestCell->SegList.size()-1].BlindSynapseAdding(0);

        }
    }



}


void layer::Three_CellListUpdater(void){
    if(Three_CellActivityList.size()>MotherBrain.max_activation_counter){
        Three_CellActivityList.pop_back();
    }
    Three_CellActivityList.insert(Three_CellActivityList.begin(),CellActivityList);
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_lower_level->CellActivityList.begin(),p_lower_level->CellActivityList.end());
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_upper_level->CellActivityList.begin(),p_upper_level->CellActivityList.end());
}
void top_layer::Three_CellListUpdater(void){
    if(Three_CellActivityList.size()>MotherBrain.max_activation_counter){
        Three_CellActivityList.pop_back();
    }
    Three_CellActivityList.insert(Three_CellActivityList.begin(),CellActivityList);
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_lower_level->CellActivityList.begin(),p_lower_level->CellActivityList.end());
}

void layer::forgetting(){
    //we forget by lowering the strength of all synapses uniformly in each round
    //synapses that are not rewarded often enough will be deleted.
    //segments that still predict sometimes will be reimbursed with blindsynapseadding.

    //segments that hardly ever predict correctly will drop below synapses_per_segment/2
    //synapses and deleted completely
    /*
    for(column& DummyColumn: ColumnList){
        for(cell& DummyCell: DummyColumn.CellList){
            for(size_t DummySegmentIndex=0;DummySegmentIndex<DummyCell.SegList.size();){
                for(size_t SynIndex=0;SynIndex<DummyCell.SegList[DummySegmentIndex].Synapse.size();){
                    DummyCell.SegList[DummySegmentIndex].Synapse[SynIndex].second-=Forgetfulness;
                    if(DummyCell.SegList[DummySegmentIndex].Synapse[SynIndex].second<=0){
                      //->  DummyCell.SegList[DummySegmentIndex].Synapse.erase(DummyCell.SegList[DummySegmentIndex].Synapse.begin()+SynIndex);
                    }else{
                        ++SynIndex;
                    }
                }// #####################################################
                if(DummyCell.SegList[DummySegmentIndex].Synapse.size()<synapses_per_segment/2 ){
                    std::vector<segment> dummydeque;
                    dummydeque.emplace_back(DummyCell,3);
                    //->dummydeque.erase(dummydeque.begin());
                    DummyCell.SegList.begin();
                    //DummyCell.SegList.erase(DummyCell.SegList.begin()+DummySegmentIndex);
                }else{
                    ++DummySegmentIndex;
                }
            }
        }
    }*/
}


void bottom_layer:: FindBestColumns(){

    std::vector<bool> sensory_input=external_input(MotherBrain.time);

    TempActColumns=std::vector<column*>();

    for(size_t index=0;index<Num_Columns;++index){

        if(sensory_input[index]==true) TempActColumns.push_back(&ColumnList[index]);

    }
}

void bottom_layer::Three_CellListUpdater(){
    if(Three_CellActivityList.size()>MotherBrain.max_activation_counter){
        Three_CellActivityList.pop_back();
    }
    Three_CellActivityList.insert(Three_CellActivityList.begin(),CellActivityList);
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_upper_level->CellActivityList.begin(),p_upper_level->CellActivityList.end());
}

double column::feed_input(void){
    //computes overlap and
    //boosted overlap; also updates running averages
    int overlap=0;
    for(auto &syn : ConnectedSynapses){//v this silly number is 0.3
        if(syn.second>minimal_overlap_under_consideration/(pillars_per_layer*active_pillers_per_pillar)){
            overlap+=static_cast<int>( syn.first->active||syn.first->expect);
        }
    }
    if(overlap<MinOverlap){
        overlap=0;
    }
    Overlap_Average= Average_Exp*Overlap_Average+ overlap;
    return overlap*boosting;
}

segment* column::BestMatchingSegmentInColumnActivateCells(void){
    //find the best matching segment of all the cells in the column,
    //which led the column to expect its activation (i.e. cell.expect==true)
    //return that segment
    double max_count=0;
    segment* pbestSegment= &CellList[0].SegList[0];
    for(cell& dummy_cell: CellList){
        if(dummy_cell.expect[0]==true){
            MotherLayer.PendingActivity.push_back(&dummy_cell);
            for(segment& dummy_segment:dummy_cell.SegList){
                double count=0;
                for(auto& remote_cell: dummy_segment.Synapse){
                    if(remote_cell.first->active[1]==true&&dummy_segment.ActivationCountdown==1)
                    {count+=remote_cell.second;}
                }
                if(count> max_count){
                    max_count=count;
                    pbestSegment=&dummy_segment;
                }
            }
        }
    }

    return pbestSegment;
}


inline void segment::BlindSynapseAdding(size_t t){
    //add all active synapses to a segment

    for(auto& remote_cell:MotherCell.MotherColumn.MotherLayer.Three_CellActivityList[t]){
        bool alreadythere=false;
        for(auto dummySynapse: Synapse){
            if(dummySynapse.first==remote_cell){
                alreadythere=true;
                break;
            }
        }
        if(!alreadythere){
            AddCell( remote_cell);
        }
    }

}

void cell::UpdateActiveSegments(void){
    //erase last element and insert list of new active segments at the front
    //determine which segments are active by weighted sum of active cells
    //to wich the segments point

    if(MotherColumn.MotherLayer.MotherBrain.max_activation_counter<ActiveSegments.size()){
        ActiveSegments.pop_back();
    }

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

    if(max>0){//something will be changed
        std::iter_swap(tempSegments.begin(),tempSegments.begin()+max-1);
    }

    ActiveSegments.insert(ActiveSegments.begin(),tempSegments);
}

segment* cell::BestSegmentInCell(size_t t){
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




void UpdateInitialiser(layer* DummyLayer){
   //prepares DummyLayer for the upcomming updates

    DummyLayer->FindBestColumns();

    double MaxActivity=DummyLayer->ActivityLogUpdateFindMaxActivity();
    if(DummyLayer!=&DummyLayer->MotherBrain.LowestLayer){
        DummyLayer->BoostingUpdate_StrenthenWeak( MaxActivity);
    }
    DummyLayer->CellExpectInitiator();


    DummyLayer->CellLearnInitiator();
}

void ThreadUpdater(layer* DummyLayer){
    //execute pending updates for DummyLayer

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    DummyLayer->ActiveColumnUpdater();

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    if(DummyLayer!=&DummyLayer->MotherBrain.LowestLayer){
        DummyLayer->ConnectedSynapsesUpdate();
    }

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    DummyLayer->CellUpdater();

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    DummyLayer->Do_SegmentUpdate();
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

}



/* #######################################################################
 *
 *
 *
 *
 *
 * debug functions
 *
 *
 *
 *
 * ########################################################################

*/


void brain::inventory(){

    bool dbool=false;
    for(size_t index1=0;index1<AllLevels.size();++index1){
        if(&AllLevels[index1]->MotherBrain!=this&&dbool==false){
            std::cout<<index1<<"-th layer does not find its mother\n";
            dbool=true;
        }
        for(size_t index2=0;index2<AllLevels[index1]->ColumnList.size();++index2){
            if(&AllLevels[index1]->ColumnList[index2].MotherLayer!=AllLevels[index1]&&dbool==false){
                std::cout<<index2<<"-th column of "<<index1<<"-th layer does not find its mother\n";
                dbool=true;
            }
            for(size_t index3=0;index3<AllLevels[index1]->ColumnList[index2].CellList.size();++index3){
                if(&AllLevels[index1]->ColumnList[index2].CellList[index3].MotherColumn!=&AllLevels[index1]->ColumnList[index2]&&dbool==false){
                    std::cout<<index3<<"-th cell of "<<index2<<"-th column of "<<index1<<"-th layer does not find its mother\n";
                    dbool=true;
                }
                for(size_t index4=0;index4<AllLevels[index1]->ColumnList[index2].CellList[index3].SegList.size();++index4){
                    if(&AllLevels[index1]->ColumnList[index2].CellList[index3].SegList[index4].MotherCell!=&AllLevels[index1]->ColumnList[index2].CellList[index3]&&dbool==false){
                        std::cout<<index4<<"-th segment of "<<index3<<"-th cell of "<<index2<<"-th column of "<<index1<<"-th layer does not find its mother\n";
                        dbool=true;
                    }

                }
            }
        }

    }
    if(dbool==false){
        std::cout<<"everything works! \n";
    }
    dbool=false;
    for(size_t i=0;i<AllLevels.size()-1;++i){
        if(AllLevels[i]->p_upper_level!=AllLevels[i+1]){
            dbool=true;
            break;
        }
        if(AllLevels[i+1]->p_lower_level!=AllLevels[i]){
            dbool=true;
            break;
        }
    }
    if(dbool==false){
        std::cout<<"everything works2! \n";
    }
}

size_t layer::finding_oneself(){
    for(size_t i=0;i<MotherBrain.AllLevels.size();++i){
        if(MotherBrain.AllLevels[i]==this){
            return i;
        }
    }
    std::cout<<"layer got lost finding itself \n";
    std::abort();

}
void layer::who_am_I(void){
    std::cout<<"I am layer "<<finding_oneself()<<std::endl;
}

size_t column::finding_oneself(){

    for(size_t i=0;i<MotherLayer.ColumnList.size();++i){
        if(&MotherLayer.ColumnList[i]==this){
            return i;
        }
    }
    throw std::invalid_argument("column got lost finding itself \n");

}
void column::who_am_I(void){

    std::cout<<"I am column "<<finding_oneself();
    std::cout<<" of layer "<<MotherLayer.finding_oneself()<<std::endl;
}
size_t cell::finding_oneself(){
    for(size_t i=0;i<MotherColumn.CellList.size();++i){
        if(&MotherColumn.CellList[i]==this){
            return i;
        }
    }
    throw std::invalid_argument("cell got lost finding itself \n");

}
void cell::who_am_I(void){
    std::cout<<"I am cell "<<finding_oneself();
    std::cout<<" of column "<<MotherColumn.finding_oneself();
    std::cout<<" of layer "<<MotherColumn.MotherLayer.finding_oneself()<<std::endl;
}
size_t segment::finding_oneself(){
    for(size_t i=0;i<MotherCell.SegList.size();++i){
        if(&MotherCell.SegList[i]==this){
            return i;
        }
    }
    throw std::invalid_argument("segment got lost finding itself \n");

}
void segment::who_am_I(void){
    std::cout<<"I am segment "<<finding_oneself()<<" of cell "<<MotherCell.finding_oneself()<<" of column "<<MotherCell.MotherColumn.finding_oneself()<<" of layer "<<MotherCell.MotherColumn.MotherLayer.finding_oneself()<<std::endl;
}

debughelper::debughelper(void):
    success_column(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    success_cell(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    activation_column(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    activation_cell(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    avg_synapses_per_segment(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>()))
{
    log.open("BrainLog.txt");
}

void debughelper::tell(std::vector<std::vector<double> >* dummyvec){
    //prints the name of one of the types of statistics debughelper is taking during runtime
    //and the averages over the data itself

    std::string name;

    if(dummyvec==&success_column){ name="success_column";
    }else if(dummyvec==&success_cell){
        name="success_cell";
    }else if (dummyvec==&activation_column){
        name="activation_column";

    } else if(dummyvec==&activation_cell){
        name="activation_cell";

    }else if(dummyvec==&avg_synapses_per_segment){
        name="avg_synapses_per_segment";

    }

    std::cout<<"This is the debughelper speaking, taking the average of "<<name<<".\n";

    for(size_t l=0;l<dummyvec->size();++l){
        std::cout<<"In layer "<<l<<" the statistics is. \n";
        for(size_t t=0;t<(*dummyvec)[l].size()/100;++t){
            std::cout<<" average number at time "<<t<<" is ";
            double average=0;
            for(size_t u=0;u<100;++u){
                average+=(*dummyvec)[l][100*t+u];
            }
            std::cout<<average/100<<std::endl;
        }
    }
}


//end of debugging functions




void brain::update(){
/*updates:
 * layer::actcolumns ←
 * layer::SegmentUpdateList ←
 * layer::CellActivityList ←
 * layer::Three_CellActivityList ←
 * column::Overlap_Average ←
 * column::ConnectedSynapses ←
 * column::active ←
 * column::ActivityLog ←
 * column::boosting ←
 * cell::ActiveSegments ←
 * cell::active ←
 * cell::expect ←
 * cell::learn ←
 * segment::Synapse ←
 * segment::EndOfSeq!!
 *
 * Do all layers in parallel, i.e. do the layers one after another
 * but don't implement updates until after the last layer.

*/    
    {//create Threads only locally in here
        if(multithreadding==true){
            std::vector<std::thread> Threads;

            Threads.reserve(layers_per_brain);
            for(layer*& DummyLayer:AllLevels){
                Threads.emplace_back(UpdateInitialiser, DummyLayer);
            }
            for(size_t i=0;i<AllLevels.size();++i){
                Threads[i].join();
            }
        }else{
            for(layer*& DummyLayer:AllLevels){
                UpdateInitialiser(DummyLayer);
            }
        }

    }//end of the first generation of threads
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;


    if(max_activation_counter_change==true){
        ++max_activation_counter;
        max_activation_counter_change=false;
    }

    {//begin of 2nd generation of threads
        //updates!
         if(multithreadding==true){
            std::vector<std::thread> Threads;
            Threads.reserve(layers_per_brain);
            for(layer*& DummyLayer:AllLevels){
                Threads.emplace_back(ThreadUpdater, DummyLayer);
            }
            for(size_t i=0;i<AllLevels.size();++i){
                Threads[i].join();
            }
         }else{
             for(layer*& DummyLayer:AllLevels){
                ThreadUpdater(DummyLayer);
             }
         }


    }
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    //needs an extra loop due to interference
    for(layer*& DummyLayer:AllLevels){

        for(column& DummyColumn : DummyLayer->ColumnList){
            for(cell& DummyCell: DummyColumn.CellList){

                DummyCell.UpdateActiveSegments();
            }
        }

        DummyLayer->Three_CellListUpdater();
    }

    //update inner clock
    ++time;
}


