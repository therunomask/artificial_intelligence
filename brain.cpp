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

    //srand(282346120965);//time(NULL));
    srand(time(NULL));
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
      time(0),
      Martin_Luther(debughelper()),
      ListOfLevels(std::vector<layer>()),
      LowestLayer(bottom_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column,*this,sensoryinput)),
      HighestLayer(top_layer(Number_of_Column_per_Layer,Number_of_Cells_per_Column,*this)),
      AllLevels(std::vector<layer*>())
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
      PendingExpectation(std::vector<cell*>()),
      PendingLearning(std::vector<cell*>())
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
      learn(2,false),
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
    MotherCell(dummysegment.MotherCell)
{
    throw std::invalid_argument("Don't copy segments! \n");

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
            if(itDummyUpdate->timer>0){
                --itDummyUpdate->timer;
                ++itDummyUpdate;
            }else{
                itDummyUpdate->AdaptingSynapses((*itdummycell)->active[0]);
                (*itdummycell)->SegmentUpdateList.erase(itDummyUpdate);
            }
        }
        if((*itdummycell)->SegmentUpdateList.size()==0){
            CellUpdateList.erase(itdummycell);
        }else{
            ++itdummycell;
        }
    }

}
void SegmentUpdate::AdaptingSynapses(bool success){
    //depending on success either we punish all the connections to (previously) active cells
    //or strengthen them.
    //if connection strength drops to zero or below we remove the corresponding synapse from the
    //segment.
    for(std::pair<cell*,double>*& dummySynapse: active_cells){
        if(success==true){
            dummySynapse->second= std::max(static_cast<double>(1),dummySynapse->second+SegmentAddress->LearnIncrement);
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
}

//void layer::SegmentUpdater(void){

//    //debughelper; statistics
//    static int success=0;
//    static int failure=0;
//    static bool loop=false;
//    bool outer_most=false;

//    bool erased=false;


//    for(cell* dummycell:CellUpdateList){
//        if(dummycell->active[0]==true){

//            if(dummycell->SegmentUpdateList[dummycell->SegmentUpdateList.size()-1]->EndOfSeq==true){
//                for(segment*& dummySegment: dummycell->SegmentUpdateList){
//             //reward segments if cell is active and activity was predicted(EndOfSeq is true)
//                    dummySegment->AdaptingSynapses(dummySegment->PositiveLearning);
//                    //debughelper; statistics
//                    ++success;
//                }
//                dummycell->SegmentUpdateList=std::vector<segment*>();
//                CellUpdateList.erase(dummycell);
//                erased=true;
//                break;
//            }
//        }
//        else if(dummycell->expect[0]==false){

//            for(segment*& dummySegment: dummycell->SegmentUpdateList){
//            //punish segments if cell is not active and not predicting; but was previously predicting
//                dummySegment->AdaptingSynapses(!(dummySegment->PositiveLearning));
//                //debughelper; statistics
//                ++failure;
//            }
//            dummycell->SegmentUpdateList=std::vector<segment*>();
//            CellUpdateList.erase(dummycell);
//            erased=true;
//            break;
//        }

//    }
//    //debughelper statistics
//    if(loop==false){
//        outer_most=true;
//    }
//    if(erased==true){

//        //debughelper statistics
//        loop=true;


//        this->SegmentUpdater();
//    }
//    //debughelper statistics
//    if(success+failure!=0&&outer_most==true){
//        MotherBrain.Martin_Luther.success_cell[finding_oneself()].push_back(static_cast<double>(success)/static_cast<double>(success+failure));
//    } else if (outer_most==true){
//        MotherBrain.Martin_Luther.success_cell[finding_oneself()].push_back(0);

//    }
//    //debughelper statistics

//    if(outer_most==true){
//        loop=false;
//        success=0;
//        failure=0;
//    }
//}




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
            dummyCell.active[1]=dummyCell.active[0];
            dummyCell.active[0]=false;
            dummyCell.expect[1]=dummyCell.expect[0];
            dummyCell.expect[0]=false;
            dummyCell.learn[1]=dummyCell.learn[0];
            dummyCell.learn[0]=false;
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
    for(cell*& pdummyCell:PendingLearning){
        pdummyCell->learn[0]=true;
    }
    PendingLearning=std::vector<cell*>();


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

//void layer::CellLearnInitiator(void){
//    //check if a cell predicted activation of a column
//    //if so, change its synapses.
//    //if not choose cell to predict same activation in the future

//    for(column*& active_pillar: ActColumns){
//        bool predicted=false;//dummy variable checks of predicting cell is found
//        bool is_chosen=false;//also dummy

//        for(cell& activePillarCell: active_pillar->CellList){

//            if(activePillarCell.expect[1]==true){

//                segment* s=NULL;
//                //find segment of cell that signified the end of a sequence
//                for(segment*& active_segment: activePillarCell.ActiveSegments[1]){
//                    if(active_segment->EndOfSeq==true){
//                        s=active_segment;
//                        break;
//                    }
//                }//if no such segment exists; continue with next cell
//                if( s==NULL){
//                    continue;
//                }

//                predicted=true;
//                PendingActivity.push_back(&activePillarCell);

//                //choose current cell to be the learning cell if it
//                //is connected to a cell with learn state on
//                for(auto& connected_cells : s->Synapse){
//                    if(connected_cells.first->learn[1]==true){
//                        PendingLearning.push_back(&activePillarCell);
//                        is_chosen=true;
//                        break;
//                    }

//                }
//            }
//        }

//        if(predicted==false){

//            for(auto& PillarCell:active_pillar->CellList){
//                PendingActivity.push_back(&PillarCell);
//            }
//        }

//        if(is_chosen==false){

//            //get best matching cell in last timestep
//            segment* pBestSegment= active_pillar->BestMatchingSegmentInColumn();
//            PendingLearning.push_back(&(pBestSegment->MotherCell));
//            pBestSegment->BlindSynapseAdding(this,1);//1= last timestep
//            pBestSegment->EndOfSeq=true;


//            pBestSegment->MotherCell.SegmentUpdateList.push_back( pBestSegment);
//            CellUpdateList.insert(& pBestSegment->MotherCell);

//        }

//    }
//}

void layer::Three_CellListUpdater(void){
    Three_CellActivityList.pop_back();
    Three_CellActivityList.insert(Three_CellActivityList.begin(),CellActivityList);
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_lower_level->CellActivityList.begin(),p_lower_level->CellActivityList.end());
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_upper_level->CellActivityList.begin(),p_upper_level->CellActivityList.end());
}
void top_layer::Three_CellListUpdater(void){
    Three_CellActivityList.pop_back();
    Three_CellActivityList.insert(Three_CellActivityList.begin(),CellActivityList);
    Three_CellActivityList[0].insert(Three_CellActivityList[0].begin(),p_lower_level->CellActivityList.begin(),p_lower_level->CellActivityList.end());
}


void bottom_layer:: FindBestColumns(){

    std::vector<bool> sensory_input=external_input(MotherBrain.time);

    TempActColumns=std::vector<column*>();

    for(size_t index=0;index<Num_Columns;++index){

        if(sensory_input[index]==true) TempActColumns.push_back(&ColumnList[index]);

    }
}

void bottom_layer::Three_CellListUpdater(){
    Three_CellActivityList.pop_back();
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
        AddCell( remote_cell);
    }

}

void cell::UpdateActiveSegments(void){
    //erase last element and insert list of new active segments at the front
    //determine which segments are active by weighted sum of active cells
    //to wich the segments point


    ActiveSegments.pop_back();

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



//void segment::AdaptingSynapses(bool positive){

//    //for positive learning reinforce all connections to active cells
//    //punish all connections to inactive cells
//    if(positive==true){
//        //no for loop, because we delete synapses that drop below 0 connectedness
//        auto dummySynapse=Synapse.begin();
//        while( dummySynapse!=Synapse.end()){
//            if(dummySynapse->first->active[0]==true){
//                //if remote cell is active-> positive learning is apropriate
//                dummySynapse->second=std::min(1.0,dummySynapse->second+LearnIncrement);
//                ++dummySynapse;
//            }
//            else{
//                //if remote cell is inactive, the synapse was useless, negative learning is apropriate
//                dummySynapse->second=dummySynapse->second-LearnIncrement;
//                if(dummySynapse->second<=0){
//                    Synapse.erase(dummySynapse);
//                }
//                else{
//                    ++dummySynapse;
//                }
//            }
//        }
//    }
//    //for negative learning punish all connections to active cells,
//    //since they predicted mistakenly
//    else{
//        auto dummySynapse=Synapse.begin();
//        while( dummySynapse!=Synapse.end()){
//            if(dummySynapse->first->active[0]==true){
//                //if remote cell is active-> decrement connectedness
//                dummySynapse->second=dummySynapse->second-LearnIncrement;
//                if(dummySynapse->second<=0){
//                    Synapse.erase(dummySynapse);
//                }
//                else{
//                    ++dummySynapse;
//                }
//            }
//            else{
//                ++dummySynapse;
//            }

//        }
//    }
//}


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

    DummyLayer->ActiveColumnUpdater();

    if(DummyLayer!=&DummyLayer->MotherBrain.LowestLayer){
        DummyLayer->ConnectedSynapsesUpdate();
    }

    DummyLayer->CellUpdater();

    DummyLayer->Do_SegmentUpdate();

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
        std::vector<std::thread> Threads;
        Threads.reserve(layers_per_brain);
        for(layer*& DummyLayer:AllLevels){
            Threads.emplace_back(UpdateInitialiser, DummyLayer);
        }
        for(size_t i=0;i<AllLevels.size();++i){
            Threads[i].join();
        }
//    for(layer*& DummyLayer:AllLevels){
//        DummyLayer->FindBestColumns();

//        double MaxActivity=DummyLayer->ActivityLogUpdateFindMaxActivity();
//        if(DummyLayer!=&DummyLayer->MotherBrain.LowestLayer){
//            DummyLayer->BoostingUpdate_StrenthenWeak( MaxActivity);
//        }
//        DummyLayer->CellExpectInitiator();


//        DummyLayer->CellLearnInitiator();
//    }

    }//end of the first generation of threads

    {//begin of 2nd generation of threads
        //updates!
        std::vector<std::thread> Threads;
        Threads.reserve(layers_per_brain);
        for(layer*& DummyLayer:AllLevels){
            Threads.emplace_back(ThreadUpdater, DummyLayer);
        }
        for(size_t i=0;i<AllLevels.size();++i){
            Threads[i].join();
        }
//        for(layer*& DummyLayer:AllLevels){

//            DummyLayer->ActiveColumnUpdater();

//            if(DummyLayer!=&LowestLayer){
//                DummyLayer->ConnectedSynapsesUpdate();
//            }

//            DummyLayer->CellUpdater();

//            DummyLayer->SegmentUpdater();
//        }

    }
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


