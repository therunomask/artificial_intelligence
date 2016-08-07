#include <iostream>
#include <vector>
#include <deque>

#include "brain_classes.cpp"





int main(int argc, char *argv[])
{
    std::vector<int> l{1,2};
    for(auto &k : l){
        std::cout<<k<<"vector entry"<<std::endl;
    }
    std::cout << "Hello World!" << std::endl;
    return 0;
}
