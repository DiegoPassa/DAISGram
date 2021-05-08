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

    //TEST POSITIVO
    /* DAISGram bright;
    bright.load_image("../images/dais.bmp");
    bright = bright.brighten(20);
    bright.save_image("../bright.bmp"); */

    /*DAISGram random;
    random.generate_random(100, 100, 3);
    random.save_image("random.bmp"); */

    //TEST POSITIVO
    /*DAISGram greenscreen;
    DAISGram bkg;
    greenscreen.load_image("./images/greenscreen/gs_4.bmp");
    bkg.load_image("./images/greenscreen/gs_4_bkg.bmp");
    int green[] {230, 226, 221};
    float threshold[] {70, 70, 55};
    greenscreen = greenscreen.greenscreen(bkg, green, threshold);
    greenscreen.save_image("greenscreen.bmp");*/

    //TEST POSITIVO 
    /*DAISGram equalized;
    equalized.load_image("./images/equalize/hill.bmp");
    equalized = equalized.equalize();
    equalized.save_image("./equalized.bmp");*/
    
    return 0;
}

   

    
