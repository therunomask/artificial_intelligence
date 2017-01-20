#include <iostream>
#include <vector>
#include <deque>

#include "brain.cpp"

std::vector<bool> senses(size_t time){
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    std::vector<bool> output(pillars_per_layer,false);
    //for(size_t index; index<pillars_per_layer;++index){
    //    output.push_back(false);
    //}
    output[time%pillars_per_layer]=true;
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;


    return output;
}



int main(int argc, char *argv[])
{

    size_t End=100;
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    brain joseph(layers_per_brain,pillars_per_layer,cells_per_column, senses);

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    for(size_t t=0;t<End;++t){
        std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;
        joseph.update();
        std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;
        std::cout<<"it is now "<<joseph.time<<" o'clock"<<std::endl;
    }

    return 0;
}
