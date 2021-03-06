#include <iostream>
#include <vector>
#include <deque>

#include "brain.h"


std::vector<bool> senses(size_t time){
    std::vector<bool> output(pillars_per_layer,false);

    //time =1;
    output[2*time%20]=true;
    output[(2*time+1)%20]=true;
    //output[(time+2)%100]=true;
    //output[(time+3)%100]=true;
    return output;
}



int main(int argc, char *argv[])
{

    clock_t timer;
    size_t End=2000;

    brain joseph(layers_per_brain,pillars_per_layer,cells_per_column, senses);


    timer=clock();
    for(size_t t=0;t<End;++t){
        joseph.update();
        std::cout<<"it is now "<<joseph.time<<" o'clock"<<std::endl;
    }
    timer= clock() - timer;
    std::cout<<"we achieve "<<CLOCKS_PER_SEC*End/static_cast<float>(timer)<<"step steps per second"<<std::endl;
    //joseph.Martin_Luther.tell(&joseph.Martin_Luther.activation_column);
//    joseph.Martin_Luther.tell(&joseph.Martin_Luther.success_cell);
    joseph.Martin_Luther.tell(&joseph.Martin_Luther.success_column);
    std::cout<<"Maximal length of chain of segments is "<<joseph.max_activation_counter<<std::endl;
    return 0;
}
