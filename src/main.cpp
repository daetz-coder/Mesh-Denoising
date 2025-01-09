// src/main.cpp

#include "denoise_obj.hpp"
#include <iostream>

int main(int argc, char* argv[]) {
    if(argc != 3) {
        std::cerr << "Usage: denoise_obj <input.obj> <output.obj>\n";
        return 1;
    }

    std::string input_obj = argv[1];
    std::string output_obj = argv[2];

    denoise_obj(input_obj, output_obj);

    return 0;
}
