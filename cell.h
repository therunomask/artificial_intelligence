#ifndef CELL_H
#define CELL_H
#include"column.h"
#include"segment.h"
#include"vector"
#include"deque"
#include"magicnumbers.h"
class column;


class cell{

public:
    cell(column& Column_to_belong_to);
    cell(const cell& dummycell);

    column& MotherColumn;
    std::deque<segment> SegList;//list of segments (net of horizontal
            // connections) of this cell
    std::vector<bool> active;//saves activity of last few timesteps
    std::vector<bool> expect;//predicts input due to past experience and dendrite information

    segment* BestSegmentInCell(size_t t);
    //debugging after this mark
    void who_am_I(void);
    size_t finding_oneself(void);
};





#endif // CELL_H
