#include <iostream>

#include "DAISGram.h"
#include "libbmp.h"
#include "tensor.h"

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
    /*DAISGram blend, im_2, toCompare;
    blend.load_image("../images/blend/blend_a.bmp");
    im_2.load_image("../images/blend/blend_b.bmp");
    toCompare.load_image("../results/blend/blend_0.50.bmp");
    blend = blend.blend(im_2, 0.50);
    blend.save_image("../blend.bmp"); 

    if (blend == toCompare)
        cout << "uguali";*/

    //TEST POSITIVO
    /* DAISGram bright;
    bright.load_image("../images/dais.bmp");
    bright = bright.brighten(20);
    bright.save_image("../bright.bmp"); */

    /*DAISGram random;
    random.generate_random(100, 100, 3);
    random.save_image("random.bmp"); */

    //TEST POSITIVO
    /* DAISGram greenscreen;
    DAISGram bkg;
    greenscreen.load_image("../images/greenscreen/gs_2.bmp");
    bkg.load_image("../images/greenscreen/gs_2_bkg.bmp");
    int green[] {144,208,49};
    float threshold[] {100,100,50};
    greenscreen = greenscreen.greenscreen(bkg, green, threshold);
    greenscreen.save_image("../greenscreen.bmp"); */

    //TEST POSITIVO
    /* DAISGram sharp;
    sharp.load_image("../images/flower_hires.bmp");
    sharp = sharp.sharpen();
    sharp.save_image("../sharp.bmp"); */

    //TEST POSITIVO
     DAISGram edge;
    edge.load_image("../images/dais.bmp");
    edge = edge.edge();
    edge.save_tensor_to_file("../edge_provafinale.txt");
    edge.save_image("../edge_g.bmp"); 

    //TEST POSITIVO
    /* DAISGram emboss;
    emboss.load_image("../images/flower_hires.bmp");
    emboss = emboss.emboss();
    emboss.save_image("../emboss.bmp"); */

    //TEST POSITIVO
    /* DAISGram smooth;
    smooth.load_image("../images/dais.bmp");
    smooth = smooth.smooth(7);
    smooth.save_image("../smooth.bmp"); */

    //TEST POSITIVO
    /*DAISGram equalized, toCompare;
    equalized.load_image("../images/dais.bmp");
    toCompare.load_image("../results/dais_equalize.bmp");
    //equalized = equalized.grayscale();
    equalized = equalized.equalize();
    equalized.save_image("../equalized.bmp"); 

    if (equalized == toCompare)
        cout << "test positivo\n"; */

    //TEST POSITIVO
    /*DAISGram sobel_h, sobel_v, f_sobel;
    sobel_h.load_image("./images/dais.bmp");
    sobel_h = sobel_h.sobel();
    sobel_h.save_image("./sobel_h.bmp");
    
    sobel_v.load_image("./images/dais.bmp");
    sobel_v = sobel_v.sobel(false);
    sobel_v.save_image("./sobel_v.bmp");

    f_sobel.load_image("./images/seba.bmp");
    f_sobel = f_sobel.full_sobel();
    f_sobel.save_image("./f_sobel.bmp");*/

    //TEST POSITIVO
    /* DAISGram flipped;
    flipped.load_image("../images/dais.bmp");
    flipped = flipped.flip(false);
    flipped.save_image("../flip.bmp");
    flipped.load_image("../images/dais.bmp");
    flipped = flipped.flip(true);
    flipped.save_image("../flip_v.bmp"); */

    //TEST POSITIVO
    /* DAISGram inverted;
    inverted.load_image("../images/dais.bmp");
    inverted = inverted.invert_colours();
    inverted.save_image("../inverted.bmp"); */

    //TEST OSITIVO
    /* DAISGram color_equalization;
    color_equalization.load_image("../images/seba.bmp");
    color_equalization = color_equalization.equalize();
    color_equalization.save_image("../c_eq.bmp"); */

    return 0;
}
