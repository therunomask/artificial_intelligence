#ifndef SEGMENT_H
#define SEGMENT_H
#include"cell.h"
#include"segmentupdate.h"
#include"magicnumbers.h"


#include<vector>
#include<deque>



class segment{
private:
    static constexpr double InitCon=initial_connectedness;
            //cells in the segment

public:
    segment(cell& Cell_to_belong_to, size_t TempActivationCountdown);
    segment(const segment& dummysegment);
    void operator=(const segment& dummysegment);


    cell& MotherCell;                            //3*activeCollumns per layer
    static constexpr double MinSynapseWeightActivity=Minimal_sum_of_synapseweights_for_activity;
    std::vector< std::pair <cell*,double>> Synapse;//pointer to adresses of cells in the segment
    //and connectedness values of the synapses in the segment
    size_t ActivationCountdown;//1 if prediction's due to this segment
            //should activate the column in question in the next
            //timestep rather than predicting prediction of activation
    std::vector<SegmentUpdate> SegmentUpdateList;

    static constexpr bool PositiveLearning= true;
    static constexpr double LearnIncrement= Learning_Increment;//possibly change to harder
                     // punishment for larger segments

    void AddCell( cell*  const newcell);


    void BlindSynapseAdding(size_t t);
    std::vector<cell *> GetActiveCells(size_t time);


    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);
    void AdaptingSynapses(bool positive, SegmentUpdate& ToBeUpdated);//increase connectedness
    //of winners, decrease connectedness
    //of all other cells in the segment, disconnect them if
    //connection too weak.
};

#endif // SEGMENT_H
