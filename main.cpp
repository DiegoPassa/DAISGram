#include <iostream>

#include "tensor.h"
#include "DAISGram.h"
#include "libbmp.h"

int main(int, char**) {
    Tensor t{3, 2, 3}, t1{3, 2, 3};
    t.init_random(1, 4);
    t1.init_random(2, 6);

    std::cout << t << t1;
    std::cout << t.concat(t1, 2);
    t = t.concat(t1, 1);
    std::cout << t;

    std::cout << t;
    std::cout << t.padding(2, 2);

    DAISGram andy;
    andy.load_image("../images/flower_hires.bmp");

    andy = andy.warhol();
    andy.save_image("../andy.bmp");

    return 0;
}

   

    
