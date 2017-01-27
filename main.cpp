#include <iostream>
#include <vector>
#include <deque>

#include "brain.cpp"


std::vector<bool> senses(size_t time){
    std::vector<bool> output(pillars_per_layer,false);

    time =1;
    output[time%20]=true;
    output[(time+1)%20]=true;
    output[(time+2)%20]=true;
    output[(time+3)%20]=true;


    return output;
}



int main(int argc, char *argv[])
{

    clock_t timer;
    size_t End=1000;
    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    brain joseph(layers_per_brain,pillars_per_layer,cells_per_column, senses);

    std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;

    timer=clock();
    for(size_t t=0;t<End;++t){
        //std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;
        joseph.update();
        //std::cout<<"still working at line "<<__LINE__<<" in function "<<__FUNCTION__<<std::endl;
        //std::cout<<"it is now "<<joseph.time<<" o'clock"<<std::endl;
    }
    timer= clock() - timer;
    std::cout<<"we achieve "<<CLOCKS_PER_SEC*End/static_cast<float>(timer)<<"step steps per second"<<std::endl;
    joseph.Martin_Luther.tell(joseph.Martin_Luther.activation_column,"activation of columns");
    return 0;
}
