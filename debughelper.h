#ifndef DEBUGHELPER_H
#define DEBUGHELPER_H

#include"brain.h"



class debughelper{
private:
public:
    debughelper(brain &BrainToBelongTo);
    std::vector<std::vector<double>> success_column;
    std::vector<std::vector<double>> success_cell;
    std::vector<std::vector<double>> activation_column;
    std::vector<std::vector<double>> activation_cell;
    std::vector<std::vector<double>> avg_synapses_per_segment;
    std::ofstream Log;
    brain& Motherbrain;

    template<typename T>
    std::ofstream& operator<< ( T message);
    //std::ofstream& operator<< ( double bla);
    void tell(std::vector<std::vector<double> > *dummyvec);
    size_t count_All_Segments(void);
    size_t count_All_Synapses(void);
    void totalColumnConnection(void);

    void HowStatic(size_t period);

    void checkConnectivity(void );
    void ThreeCellActivityTester(std::vector<std::vector<cell *> > ThreeCellActivityList,size_t LayerNumber);
};


#endif // DEBUGHELPER_H
