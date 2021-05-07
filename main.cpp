#include <iostream>

#include "tensor.h"
#include "DAISGram.h"
#include "libbmp.h"

int main(int, char**) {

    //TEST POSITIVO
    /* DAISGram andy;
    andy.load_image("../images/flower_hires.bmp");

    andy = andy.warhol();
    andy.save_image("../andy.bmp"); */

    //TEST POSITIVO
    /* t.read_file("../tensors/t_4_30_2_progressive.txt");

    std::cout << t;

    t.write_file("../tensoree.txt"); */

    //TEST POSITIVO
    /* DAISGram grayscale;
    grayscale.load_image("../images/dais.bmp");
    grayscale = grayscale.grayscale();
    grayscale.save_image("../grayscale.bmp"); */

    //TEST POSITIVO
    /* DAISGram blend, im_2;
    blend.load_image("../images/blend/blend_a.bmp");
    im_2.load_image("../images/blend/blend_b.bmp");
    blend = blend.blend(im_2, 0.75);
    blend.save_image("../blend.bmp"); */

    return 0;
}

   

    
