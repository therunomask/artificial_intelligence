#ifndef MAGICNUMBERS_H
#define MAGICNUMBERS_H


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//          optimize with respect to these numbers
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//general dimensions of the system
#define layers_per_brain                            6//8 with multithreadding
#define active_pillers_per_pillar                   0.1 //0.02
#define cells_per_column                            3
#define pillars_per_layer                           2 / active_pillers_per_pillar
#define segments_per_cell                           cells_per_column/active_pillers_per_pillar/7
#define synapses_per_segment                        3*active_pillers_per_pillar*pillars_per_layer
//end of general dimensions

//synapse
#define initial_connectedness                       0.5
#define Minimal_sum_of_synapseweights_for_activity  0.3*cells_per_column*active_pillers_per_pillar*pillars_per_layer
#define Learning_Increment                          0.07

//column
#define column_geometric_factor                     99/100.0
#define minimal_overlap_under_consideration         0.25*pillars_per_layer*active_pillers_per_pillar
#define Initial_Activity_log                        2.5
//(1/(1-column_geometric_factor) ^(1/active_pillers_per_pillar)))
#define Initial_Overlap_Average                     Initial_Activity_log* active_pillers_per_pillar*pillars_per_layer
#define ColumnInitialConnectedness                  0.5


//layer
#define Learning_Increment_spatial                  0.1
#define Maximum_Connectedness                       2//1/active_pillers_per_pillar
#define Average_Overlap_lower_boundary              0.01
#define Homogenous_Overlap_Increment                0.0
#define Forgetfulness                               Learning_Increment/60 //highly dependent on current model 3*200
                                                                                                               // ^ frequency of input repitition
//magic boosting function
#define maximum_boosting                            1.0//3.0




#endif // MAGICNUMBERS_H
