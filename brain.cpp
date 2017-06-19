#include"brain.h"
#include"bottomlayer.h"
#include"cell.h"
#include"column.h"
#include"toplayer.h"
#include"debughelper.h"
#include"layer.h"
#include"segment.h"
#include"segmentupdate.h"
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
#include <assert.h>

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
            std::vector<size_t> Indices;
            for(size_t i=0;i<DummyLayer->Num_Columns;++i){
                Indices.push_back(i);
            }
            for(size_t SynapseIndex=0;SynapseIndex<Maximum_Connectedness;SynapseIndex++){

                size_t chosenIndex=rand()%Indices.size();
    //                                                                                                                                                                                                          uniformly distributed between MinoverLap/2(0.15) and 1.5 Minoverlap(0.45)
                DummyPillar.ConnectedSynapses.push_back(std::pair<column*,double>(&(DummyLayer->p_lower_level->ColumnList[Indices[chosenIndex]]),0.3*(static_cast<double>(rand())/RAND_MAX+0.5)*minimal_overlap_under_consideration/(pillars_per_layer*active_pillers_per_pillar)));
                Indices.erase(Indices.begin()+ chosenIndex);
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
      Martin_Luther(*this),
      max_activation_counter(5),
      max_activation_counter_change(false)
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
    srand(282346120965);
    //srand(time(NULL));
    for (column*& dummycolumn: ActColumns){
        dummycolumn=&ColumnList[rand()%ColumnList.size()];
    }

}

layer::layer(const layer &dummylayer):
    MotherBrain(dummylayer.MotherBrain),
    Num_Columns(dummylayer.Num_Columns)

{

    throw std::invalid_argument("Don't copy layers \n");

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
    active(std::vector<bool>(2,false)),
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
      active(2,false),
      expect(2,false)
{
    //SegList.reserve(synapses_per_segment);
    for(size_t index=0; index<segments_per_cell;++index){
        SegList.emplace_back(*this,1);
    }

}

cell::cell(const cell& dummycell):
    MotherColumn(dummycell.MotherColumn)
{
    throw std::invalid_argument("Don't copy cells! \n");
}

SegmentUpdate::SegmentUpdate(std::vector<cell*> active_cells, size_t countdown):
    active_cells(active_cells),
    timer(countdown)
{

}

segment::segment(cell& Cell_to_belong_to, size_t TempActivationCountdown)
    :
      MotherCell(Cell_to_belong_to),
      Synapse(std::vector< std::pair <cell*,double>>() ),//still to be initialized afterward
      ActivationCountdown(TempActivationCountdown),
      SegmentUpdateList(std::vector<SegmentUpdate>())
{

}
segment::segment(const segment &dummysegment):
    MotherCell(dummysegment.MotherCell),
    Synapse(dummysegment.Synapse),
    ActivationCountdown(dummysegment.ActivationCountdown),
    SegmentUpdateList(dummysegment.SegmentUpdateList)
{

}

void segment::operator=(const segment& rhsSegment){
    //we only want this to be used in forgetting(). So only in erasing
    //segments from the same cell.
    if(&MotherCell!=&rhsSegment.MotherCell){
        std::cout<<"segments belong to different cells! Abort!\n";
    }
    assert(&MotherCell==&rhsSegment.MotherCell);
    this->Synapse=rhsSegment.Synapse;
    this->ActivationCountdown=rhsSegment.ActivationCountdown;

}

//segment::~segment(void){
//}

std::vector<cell*>  segment::GetActiveCells(size_t time){
    std::vector<cell*> tempactive_cells;

    for(std::pair <cell*,double>& dummySynapse: Synapse){
        if(dummySynapse.first->active[time]==true){
            tempactive_cells.emplace_back(dummySynapse.first);
        }
    }
    return tempactive_cells;
}

inline void segment::AddCell( cell*  const newcell){
    //adds new cell to segment with
    //standard connectedness close to disconnection
    double con=InitCon;
    Synapse.push_back(std::pair <cell*,double>(newcell,con));
    //CellCon.push_back(InitCon);
}

void layer::FindBestColumns(void){
    //finding #DesiredLocalActivity highest overlapping columns

    std::priority_queue<std::pair<double, column*>,std::vector<std::pair<double,column*>>,std::greater<std::pair<double,column*>>> winner;
    //priority queue orderes its elements automatically
    //long definition to make priority queue order its elements in increasing order

    for (size_t i = 0; i < Num_Columns; ++i) {
        double overlap=ColumnList[i].feed_input();
        if(overlap>0){
            if(winner.size()<DesiredLocalActivity){
                winner.push(std::pair<double, column*>(overlap, &ColumnList[i]));
            }
            else if(overlap >winner.top().first){
                  winner.pop();
                  winner.push(std::pair<double, column*>(overlap, &ColumnList[i]));
            }
        }
    }//winner now contains #DesiredLocalActivity highest overlapping columns

    size_t NumberOfWinners=winner.size();
    for (size_t i = 0; i<NumberOfWinners; ++i) {
        TempActColumns[i]=winner.top().second;
        winner.pop();
    }
    if(NumberOfWinners<DesiredLocalActivity){
        //if there are not enough winners, we create future winners
        std::vector<column*> loser(ColumnSynapseAdding(DesiredLocalActivity-NumberOfWinners));

        for(size_t i=NumberOfWinners;i<DesiredLocalActivity;++i){
            TempActColumns[i]=loser[i-NumberOfWinners];
        }

    }



//    for(size_t i=0;i<TempActColumns.size();++i){
//        MotherBrain.Martin_Luther<<"at time "<<MotherBrain.time<<" ";
//        TempActColumns[i]->who_am_I();
//    }

}

std::vector<column * > layer::ColumnSynapseAdding(size_t ColumnsToFind){
    //We would like to asign to columns, which were inactive so far, the unencountered input.

std::priority_queue<std::pair<double, column *>, std::vector<std::pair<double, column *> > > LazyColumns;
    for(column& DummyColumn:ColumnList){
        //sort by inactivity
        double AverageActivity = DummyColumn.ActivityLog;
        if(LazyColumns.size()<ColumnsToFind){
            LazyColumns.push(std::pair<double, column*>(AverageActivity, &DummyColumn));
        }
        else if(AverageActivity <LazyColumns.top().first){
              LazyColumns.pop();
              LazyColumns.push(std::pair<double, column*>(AverageActivity, &DummyColumn));
        }
    }

    //copy queue into vector, for iteration
    std::vector<column*> ReturnColumns;
    while (LazyColumns.size()!=0){
        ReturnColumns.push_back(LazyColumns.top().second);
        LazyColumns.pop();
    }

    //add synapses to these columns
    for(auto& DummyColumn:ReturnColumns){
        for(column*& RemoteColumn: p_lower_level->ActColumns){
            bool AlreadyThere=false;
            for(auto& DummySynapse: DummyColumn->ConnectedSynapses){
                if(RemoteColumn==DummySynapse.first){
                    AlreadyThere=true;
                    break;
                }
            }
            if(AlreadyThere==false){
                DummyColumn->ConnectedSynapses.push_back(std::pair<column*,double>(RemoteColumn,ColumnInitialConnectedness));
            }
        }
    }


    return ReturnColumns;
}

void layer::ActiveColumnUpdater(void){
    //implement saved changes

    for(column& dummycolumn:ColumnList){
        dummycolumn.active[1]=dummycolumn.active[0];
        dummycolumn.active[0]=false;
    }


    for(size_t i=0;i<TempActColumns.size();++i){
        ActColumns[i]=TempActColumns[i];
        ActColumns[i]->active[0]=true;
    }


    //debughelper statistics
    //MotherBrain.Martin_Luther.activation_column[finding_oneself()].push_back(TempActColumns.size());


}


void layer::ConnectedSynapsesUpdate(void){
    //strengthen connections to columns that were active
    //weaken connections homogenously only if there are too many strong connections

    for(auto& pdummyColumn: ActColumns){

        double TotalConecctedness=0;

        for(auto& dummy_connected_synapse : pdummyColumn->ConnectedSynapses){
            if(dummy_connected_synapse.first->active[1]==true){
                dummy_connected_synapse.second=std::min(dummy_connected_synapse.second+CondsInc/(pillars_per_layer*active_pillers_per_pillar),1.0);
            }
            TotalConecctedness+=dummy_connected_synapse.second;
        }
        if(TotalConecctedness>MaximumConnectedness){
            double decrement=(TotalConecctedness-MaximumConnectedness)/pdummyColumn->ConnectedSynapses.size();
            for(size_t i=0 ; i<pdummyColumn->ConnectedSynapses.size();){
                pdummyColumn->ConnectedSynapses[i].second-=decrement;
                if(pdummyColumn->ConnectedSynapses[i].second<1/100000000.0){
                    pdummyColumn->ConnectedSynapses.erase(pdummyColumn->ConnectedSynapses.begin()+i);
                }else{
                    ++i;
                }
            }
        }

    }

}
double layer::ActivityLogUpdateFindMaxActivity(void){
    double MaxActivity = 0;

    for(auto& dummy_pillar:ColumnList){
        //change activity log
        dummy_pillar.ActivityLog =dummy_pillar.ActivityLog* dummy_pillar.Average_Exp+dummy_pillar.active[0];

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
    for(column& pillar : ColumnList){
        //strenthen pillars that never overlap uniformly
        if(pillar.tell_overlap_average() < Max_overlap){
            for(auto& con: pillar.ConnectedSynapses){
                con.second+=SpecialOverlapIncrement;
                con.second=std::min(con.second,1.0);
            }
        }

    }//end of learning

}


void layer::Do_SegmentUpdate(){
    //update timers
    //if timer drops to zero or below implement pending update
    //also simultaneously delete used SegmentUpdate objects from the vector
    //and delete element of CellUpdateList if there are no more pending updates for this cell
    size_t success=0;
    size_t failure=0;

    for(std::vector<cell*>::iterator itdummycell=CellUpdateList.begin();itdummycell!=CellUpdateList.end();){
        bool StayInCellUpdateList=false;
        for(segment& dummysegment:(*itdummycell)->SegList){
            for(std::vector<SegmentUpdate>::iterator itDummyUpdate=dummysegment.SegmentUpdateList.begin();itDummyUpdate!=dummysegment.SegmentUpdateList.end();){

                if(itDummyUpdate->timer>1){

                    --itDummyUpdate->timer;
                    ++itDummyUpdate;
                    PendingExpectation.emplace_back(*itdummycell);

                }else{

                    dummysegment.AdaptingSynapses((*itdummycell)->active[0],*itDummyUpdate);
                    if((*itdummycell)->active[0]==true){
                        ++success;
                    }else{
                        ++failure;
                    }
                    dummysegment.SegmentUpdateList.erase(itDummyUpdate);

                }
            }

            if(dummysegment.SegmentUpdateList.size()!=0){
               StayInCellUpdateList=true;
            }

        }
        if(StayInCellUpdateList==false){
            CellUpdateList.erase(itdummycell);

        }else{

            ++itdummycell;

        }
    }


    MotherBrain.Martin_Luther.success_column[finding_oneself()].push_back(static_cast<double>(success)/static_cast<double>(success+failure));

}

void segment::AdaptingSynapses(bool success, SegmentUpdate& ToBeUpdated){
    //depending on success either we punish all the connections to (previously) active cells
    //or strengthen them.
    //if connection strength drops to zero or below we remove the corresponding synapse from the
    //segment.

    for(cell*& remote_cell: ToBeUpdated.active_cells){
        for(auto dummySynapse:Synapse){

            if(remote_cell==dummySynapse.first){
                if(success==true){
                    dummySynapse.second= std::min(static_cast<double>(1),dummySynapse.second+LearnIncrement);
                    if(Synapse.size()<synapses_per_segment){
                        BlindSynapseAdding(ActivationCountdown-1);
                    }
                }else{
                    dummySynapse.second= dummySynapse.second- LearnIncrement;
                    if(dummySynapse.second<=0){
                        for(size_t remoteSynapse=0;remoteSynapse<Synapse.size();++remoteSynapse ){
                            if(dummySynapse.first==Synapse[remoteSynapse].first){
                                Synapse.erase(Synapse.begin()+remoteSynapse);
                                break;
                            }
                        }
                    }

                }
                break;
            }
        }
    }


}





void layer::CellExpectInitiator( void ){
//because we would like to predict one more timestep
//, after marking an active column with an active cell with an active segment for
   //updates, we create another segment. the new segment should predict the prediction
    //of the old one

    for(column& pillars:ColumnList){
        for(cell& dummycell:pillars.CellList){

            //check whether cell currently has an active segment           
            segment* BestSegment= dummycell.BestSegmentInCell(0);
            if(BestSegment!=NULL){
                //cell predicts now, because of active segment
                PendingExpectation.push_back(&dummycell);

                //this segment is supposed to learn
                CellUpdateList.push_back(&dummycell);
                BestSegment->SegmentUpdateList.emplace_back(BestSegment->GetActiveCells(0),BestSegment->ActivationCountdown);


                if(dummycell.expect[0]==false){
                    //if prediction is unexpected,
                    //create a new Segment that matches the activity of the previous
                    //timestep. Do blind synapse adding for this segment,
                    if(MotherBrain.max_activation_counter==BestSegment->ActivationCountdown+1){
                        MotherBrain.max_activation_counter_change=true;
                    }

                    dummycell.SegList.emplace_back(dummycell,BestSegment->ActivationCountdown+1);
                    dummycell.SegList[dummycell.SegList.size()-1].BlindSynapseAdding(1);

                }
            }
        }
    }
}



void layer::CellUpdater(void){


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




    for(column*& DummyColumn:TempActColumns){
        if(DummyColumn->expect==true){
            //this function finds the segment that fits best the
            //activity of the last timestep with SequenceCounter==1
            //and also leads all the expecting Cells to activation
            for(size_t timer=0; timer<Three_CellActivityList.size();++timer){
                segment* BestSegment=DummyColumn->BestMatchingSegmentInColumnActivateCells(timer);
                BestSegment->SegmentUpdateList.emplace_back(BestSegment->GetActiveCells(timer),BestSegment->ActivationCountdown);
            }
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

    for(column& DummyColumn: ColumnList){
        for(cell& DummyCell: DummyColumn.CellList){
            for(size_t DummySegmentIndex=0;DummySegmentIndex<DummyCell.SegList.size();){
                for(size_t SynIndex=0;SynIndex<DummyCell.SegList[DummySegmentIndex].Synapse.size();){
                    if(DummyCell.SegList[DummySegmentIndex].Synapse[SynIndex].second<1&&
                            DummyCell.SegList[DummySegmentIndex].SegmentUpdateList.size()==0){

                        DummyCell.SegList[DummySegmentIndex].Synapse[SynIndex].second-=Forgetfulness;
                    }
                    if(DummyCell.SegList[DummySegmentIndex].Synapse[SynIndex].second<=0){

                        DummyCell.SegList[DummySegmentIndex].Synapse.erase(DummyCell.SegList[DummySegmentIndex].Synapse.begin()+SynIndex);

                    }else{
                        ++SynIndex;
                    }
                }

                if(DummyCell.SegList[DummySegmentIndex].Synapse.size()<synapses_per_segment/2 ){

                    DummyCell.SegList.erase(DummyCell.SegList.begin()+DummySegmentIndex);

                }else{
                    ++DummySegmentIndex;
                }
            }
        }
    }
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
            overlap+=static_cast<int>( syn.first->active[0]);//||syn.first->expect);
        }
    }
    if(overlap<MinOverlap){
        overlap=0;
    }
    Overlap_Average= Average_Exp*Overlap_Average+ overlap;
    return overlap*boosting;
}

segment* column::BestMatchingSegmentInColumnActivateCells(size_t time){
    //find the best matching segment of all the cells in the column,
    //which led the column to expect its activation (i.e. cell.expect==true)
    //return that segment.
    //also activate all the cells that expected activation of the column
    double max_count=0;
    segment* pbestSegment=NULL;
    for(cell& dummyCell:CellList){
        if(dummyCell.SegList.size()>0){
            pbestSegment= &dummyCell.SegList[0];
            break;
        }
    }

    if(pbestSegment==NULL){
        throw std::invalid_argument("expecting column has no segment \n");
    }
    for(cell& dummy_cell: CellList){
        if(dummy_cell.expect[0]==true){
            bool DontDuplicateCells=false;
            for(cell*& extraCell: MotherLayer.PendingActivity){
                if(extraCell==&dummy_cell){
                    DontDuplicateCells=true;
                }
            }
            if(DontDuplicateCells==false){
                MotherLayer.PendingActivity.push_back(&dummy_cell);
            }
            for(segment& dummy_segment:dummy_cell.SegList){
                double count=0;
                for(auto& remote_cell: dummy_segment.Synapse){
                    if(remote_cell.first->active[time]==true&&dummy_segment.ActivationCountdown==time+1)
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


inline void segment::BlindSynapseAdding(size_t time){
    //add all active synapses to a segment

    for(auto& remote_cell:MotherCell.MotherColumn.MotherLayer.Three_CellActivityList[time]){
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


segment* cell::BestSegmentInCell(size_t t){
    //finds the Segment which matches best at time t
    //and returns a pointer to it
    double max=0;
    segment* pBestSegment=NULL;
    for(segment& dummysegment: SegList){
        double sum=0;
        for(auto& dummysynapse: dummysegment.Synapse){

            if(dummysynapse.first->active[t]==true){
                //add activity of synapse only if remote cell is active
                sum+=dummysynapse.second;
            }
        }
        if(sum>=max&&sum>=dummysegment.MinSynapseWeightActivity){
            max=sum;
            pBestSegment=&dummysegment;
        }
    }
    return pBestSegment;
}




void UpdateInitialiser(layer* DummyLayer){

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
    if(DummyLayer!=&DummyLayer->MotherBrain.LowestLayer){
        DummyLayer->ConnectedSynapsesUpdate();
    }

    DummyLayer->ActiveColumnUpdater();

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
    for(layer*& dummylevel:AllLevels){
        dummylevel->who_am_I();
        for(column& dummycolumn:dummylevel->ColumnList){
            dummycolumn.who_am_I();
            for(auto& dummysynapse : dummycolumn.ConnectedSynapses){
                Martin_Luther<<dummysynapse.second<<" and pointing to ";
                dummysynapse.first->who_am_I();
            }
            for(cell& dummycell:dummycolumn.CellList){
                dummycell.who_am_I();
                for(segment& dummysegment:dummycell.SegList){
                    dummysegment.who_am_I();
                    for(auto& dummysynapse : dummysegment.Synapse){
                        Martin_Luther<<dummysynapse.second<<" and pointing to ";
                        dummysynapse.first->who_am_I();
                    }
                }
            }
        }
    }
}

void brain::pointerCheck(){

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
    //std::cout<<"layer got lost finding itself \n";
    std::abort();

}
void layer::who_am_I(void){
    MotherBrain.Martin_Luther<<"I am layer "<<finding_oneself()<<std::endl;
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

    MotherLayer.MotherBrain.Martin_Luther<<"I am column "<<finding_oneself();
    MotherLayer.MotherBrain.Martin_Luther<<" of layer "<<MotherLayer.finding_oneself()<<std::endl;
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
    MotherColumn.MotherLayer.MotherBrain.Martin_Luther<<"I am cell "<<finding_oneself();
    MotherColumn.MotherLayer.MotherBrain.Martin_Luther<<" of column "<<MotherColumn.finding_oneself();
    MotherColumn.MotherLayer.MotherBrain.Martin_Luther<<" of layer "<<MotherColumn.MotherLayer.finding_oneself()<<std::endl;
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
    MotherCell.MotherColumn.MotherLayer.MotherBrain.Martin_Luther<<"I am segment with activationcountdown "<<ActivationCountdown<<", "<<finding_oneself()<<" of cell "<<MotherCell.finding_oneself()<<" of column "<<MotherCell.MotherColumn.finding_oneself()<<" of layer "<<MotherCell.MotherColumn.MotherLayer.finding_oneself()<<std::endl;
}

debughelper::debughelper(brain& BrainToBelongTo):
    success_column(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    success_cell(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    activation_column(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    activation_cell(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    avg_synapses_per_segment(std::vector<std::vector<double>>(layers_per_brain,std::vector<double>())),
    Log("BrainLog.txt"),
    Motherbrain(BrainToBelongTo)
{
    Log<<"oeffne Logfile\n";
    Log.flush();
}



template<typename T>
std::ofstream& debughelper::operator<<(T message){
    Log<<message;
    Log.flush();
    return Log;
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
    Log<<"This is the debughelper speaking, taking the average of "<<name<<".\n";


    for(size_t l=0;l<dummyvec->size();++l){
        for(size_t t=0;t<(*dummyvec)[l].size()/100;++t){
            double average=0;
            size_t nans=0;
            for(size_t u=0;u<100;++u){
                if((*dummyvec)[l][100*t+u]!=(*dummyvec)[l][100*t+u]){
                    //std::cout<<"NaN at "<<100*t+u<<std::endl;
                    ++nans;
                }else{
                    average+=(*dummyvec)[l][100*t+u];
                }

            }
            if(nans!=100){
                Log<<"layer "<<l<<"time "<<t*100<<"avareage "<<average/(100-nans)<<std::endl;
            }
        }
    }
}

size_t debughelper::count_All_Segments(void){
    size_t allSegments=0;
    for(layer*& dummyLayer:Motherbrain.AllLevels){
        for(column& dummyColumn:dummyLayer->ColumnList){
            for(cell& dummyCell: dummyColumn.CellList){
                allSegments+=dummyCell.SegList.size();
            }
        }
    }
    return allSegments;
}

size_t debughelper::count_All_Synapses(void){
    size_t allSynapses=0;
    for(layer*& dummyLayer:Motherbrain.AllLevels){
        for(column& dummyColumn:dummyLayer->ColumnList){
            for(cell& dummyCell: dummyColumn.CellList){
                for(segment& dummySegment: dummyCell.SegList){
                    allSynapses+=dummySegment.Synapse.size();
                }
            }
        }
    }
    return allSynapses;
}

void debughelper::checkConnectivity(void ){
    for(layer*& dummyLayer: Motherbrain.AllLevels){
        for(column& dummycolumn: dummyLayer->ColumnList){
            for(auto& connection:dummycolumn.ConnectedSynapses){
                if(connection.second<0||connection.second>1){
                    Log<<"this column has a wrong connection of "<<connection.second<<std::endl;
                    dummycolumn.who_am_I();

                }
                for(cell& dummycell:dummycolumn.CellList){
                    for(segment& dummysegment:dummycell.SegList){
                        for(auto& dummysynapse: dummysegment.Synapse){
                            if(dummysynapse.second<0||dummysynapse.second>1){
                                Log<<"this segment has a wrong connection of "<<dummysynapse.second<<std::endl;
                                dummysegment.who_am_I();

                            }
                        }
                    }
                }
            }
        }
    }
}



void debughelper::ThreeCellActivityTester(std::vector<std::vector<cell*>> ThreeCellActivityList,size_t LayerNumber){
    static std::vector<std::vector<std::vector<cell*>>> ThreeCellActivityListCopy(Motherbrain.AllLevels.size(),std::vector<std::vector<cell*>>(Motherbrain.max_activation_counter,std::vector<cell*>()));
    if(Motherbrain.time>Motherbrain.max_activation_counter){
        for(size_t iTime=0;iTime<ThreeCellActivityList.size()-1;++iTime){
            if(ThreeCellActivityList[iTime+1].size()!=ThreeCellActivityListCopy[LayerNumber][iTime].size()){
                Log<<"It is now "<<Motherbrain.time<<" the ThreeCellActivityList changed wrongly!\n";
                Log<<"It has the wrong length at position "<<iTime+1<<std::endl;
                throw std::invalid_argument("ThreeCellActivityList Changed wrongly.\n");

            }

            for(size_t iCell=0;iCell<ThreeCellActivityList[0].size();++iCell){
                if(ThreeCellActivityList[0][iCell]->active[0]==false){
                    Log<<"Cell:\n";
                    ThreeCellActivityList[0][iCell]->who_am_I();
                    Log<< " in ThreeCellActivityList is not active!!!\n ";
                }
            }
            for(size_t iCell=0;iCell<ThreeCellActivityList[iTime+1].size();++iCell){


                if(ThreeCellActivityList[iTime+1][iCell]!=ThreeCellActivityListCopy[LayerNumber][iTime][iCell]){
                    Log<<"It is now "<<Motherbrain.time<<" the ThreeCellActivityList changed wrongly!\n";
                    Log<<"It has the wrong length at position "<<iTime+1<<", "<<iCell<<std::endl;
                    throw std::invalid_argument("ThreeCellActivityList Changed wrongly.\n");
                }
            }
        }
    }
    ThreeCellActivityListCopy[LayerNumber]=ThreeCellActivityList;
}

void debughelper::totalColumnConnection(){
    double totalConnection=0;
    for(layer& DummyLayer: Motherbrain.ListOfLevels){
        for(column& dummycolumn: DummyLayer.ColumnList){
            for(auto& DummySynapse: dummycolumn.ConnectedSynapses){
                totalConnection+= DummySynapse.second;
            }
        }
    }
    Log<<"It is now "<<Motherbrain.time<<" the current total column connection is "<<totalConnection<<std::endl;
}


//end of debugging functions


void debughelper::HowStatic(size_t period){
    //finish comparison between activity of cells and columns now and period timesteps ago
    static size_t RealPeriod= period;
    static std::vector<std::vector<std::vector<bool>>> ColumnAct(Motherbrain.AllLevels.size(),std::vector<std::vector<bool>>(Motherbrain.AllLevels[0]->ColumnList.size(),std::vector<bool>(RealPeriod,false)));
    static std::vector<std::vector<std::vector<std::vector< bool>>>> CellAct(Motherbrain.AllLevels.size(), std::vector<std::vector<std::vector<bool>>>(Motherbrain.AllLevels[0]->ColumnList.size(),std::vector<std::vector<bool>>(Motherbrain.AllLevels[0]->ColumnList[0].CellList.size(),std::vector<bool>(RealPeriod,false))));

    assert(period==RealPeriod);


    size_t ColumnCount=0;
    size_t CellCount=0;

    for(size_t LayerIndex=0;LayerIndex<Motherbrain.AllLevels.size();++LayerIndex){
        for(size_t ColumnIndex=0;ColumnIndex<Motherbrain.AllLevels[LayerIndex]->ColumnList.size();++ColumnIndex){
            ColumnAct[LayerIndex][ColumnIndex].insert(ColumnAct[LayerIndex][ColumnIndex].begin(), Motherbrain.AllLevels[LayerIndex]->ColumnList[ColumnIndex].active[0]);
            if(ColumnAct[LayerIndex][ColumnIndex][0]==ColumnAct[LayerIndex][ColumnIndex].back()){
                ++ColumnCount;
            }
            ColumnAct[LayerIndex][ColumnIndex].pop_back();
            for(size_t CellIndex=0;CellIndex<Motherbrain.AllLevels[LayerIndex]->ColumnList[ColumnIndex].CellList.size();++CellIndex){
                CellAct[LayerIndex][ColumnIndex][CellIndex].insert(CellAct[LayerIndex][ColumnIndex][CellIndex].begin(),Motherbrain.AllLevels[LayerIndex]->ColumnList[ColumnIndex].CellList[CellIndex].active[0]);
                if(CellAct[LayerIndex][ColumnIndex][CellIndex][0]==CellAct[LayerIndex][ColumnIndex][CellIndex].back()){
                    ++CellCount;
                }
                else{
                    Log<<"Cell number: "<<LayerIndex<<","<<ColumnIndex<<","<<CellIndex<<" changed \n";
                }
                CellAct[LayerIndex][ColumnIndex][CellIndex].pop_back();
            }
        }
    }
    double ColumnAverage=static_cast<double>( ColumnCount)/(ColumnAct[0].size()*ColumnAct.size());
    double CellAverage= static_cast<double>(CellCount)/(CellAct[0][0].size()*CellAct[0].size()*CellAct.size());
    Log<<"At time "<<Motherbrain.time<<", the columns are as static as "<<ColumnAverage<<", while the cells are as static as "<<CellAverage<<std::endl;
}

void brain::update(){
    //Martin_Luther.checkConnectivity();


    if(time%10==3&&time>400){
        std::cout<<"bla \n";
    }
    Martin_Luther.HowStatic(10);
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
    //needs an extra loop due to interference

    //needs an extra loop due to interference
    for(layer*& DummyLayer:AllLevels){
        DummyLayer->forgetting();
        DummyLayer->Three_CellListUpdater();
    }

    for(size_t LayerIndex=0;LayerIndex<AllLevels.size();++LayerIndex ){
        Martin_Luther.ThreeCellActivityTester(AllLevels[LayerIndex]->Three_CellActivityList,LayerIndex);
    }

    //update inner clock
    ++time;
}


