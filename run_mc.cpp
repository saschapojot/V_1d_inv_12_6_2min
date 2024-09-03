#include "./mc_subroutine/mc_read_load_compute.hpp"
#include "./potentialFunction/potentialFunctionPrototype.hpp"

int main(int argc, char *argv[])
{

    if (argc != 2) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }

    std::string coefsStr="25,80,20,30,0.824655597149,1.07465559715,0.07,0.05,50,0.9,5";

    auto mcObj=mc_computation(std::string(argv[1]));
    mcObj.init_and_run();
}