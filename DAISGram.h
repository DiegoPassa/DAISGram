#ifndef DAISGRAM_H
#define DAISGRAM_H

#include <iostream>
#include <string>

#include "dais_exc.h"
#include "libbmp.h"
#include "tensor.h"

using namespace std;

class DAISGram {
   private:
    Tensor data;

   public:
    //DAISGram();

    //~DAISGram();

    /**
     * Load a bitmap from file
     *
     * @param filename String containing the path of the file
     */
    void load_image(string filename);

    /**
     * Save a DAISGram object to a bitmap file.
     * 
     * Data is clamped to 0,255 before saving it.
     *
     * @param filename String containing the path where to store the image.
     */
    void save_image(string filename);

    /**
     * Get rows
     *
     * @return returns the number of rows in the image
     */
    int getRows() const;

    /**
     * Get columns
     *
     * @return returns the number of columns in the image
     */
    int getCols() const;

    /**
     * Get depth
     *
     * @return returns the number of channels in the image
     */
    int getDepth() const;

    /**
     * Brighten the image
     * 
     * It sums the bright variable to all the values in the image.
     * 
     * Before returning the image, the corresponding tensor should be clamped in [0,255]
     * 
     * @param bright the amount of bright to add (if negative the image gets darker)
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram brighten(float bright) const;

    /**
     * Create a grayscale version of the object
     * 
     * A grayscale image is produced by substituting each pixel with its average on all the channel
     *  
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram grayscale() const;

    /**
     * Create a Warhol effect on the image
     * 
     * This function returns a composition of 4 different images in which the:
     * - top left is the original image
     * - top right is the original image in which the Red and Green channel are swapped
     * - bottom left is the original image in which the Blue and Green channel are swapped
     * - bottom right is the original image in which the Red and Blue channel are swapped
     *  
     * The output image is twice the dimensions of the original one.
     * 
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram warhol() const;

    /**
     * Sharpen the image
     * 
     * This function makes the image sharper by convolving it with a sharp filter
     * 
     * filter[3][3]
     *    0  -1  0
     *    -1  5 -1
     *    0  -1  0
     *  
     * Before returning the image, the corresponding tensor should be clamped in [0,255]
     * 
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram sharpen() const;

    /**
     * Emboss the image
     * 
     * This function makes the image embossed (a light 3D effect) by convolving it with an
     * embossing filter
     * 
     * filter[3][3]
     *    -2 -1  0
     *    -1  1  1
     *     0  1  2
     * 
     * Before returning the image, the corresponding tensor should be clamped in [0,255]
     *  
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram emboss() const;

    /**
     * Smooth the image
     * 
     * This function remove the noise in an image using convolution and an average filter
     * of size h*h:
     * 
     * c = 1/(h*h)
     * 
     * filter[3][3]
     *    c c c
     *    c c c
     *    c c c
     *  
     * @param h the size of the filter
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram smooth(int h = 3) const;

    /**
     * Edges of an image
     * 
     * This function extract the edges of an image by using the convolution 
     * operator and the following filter
     * 
     * 
     * filter[3][3]
     * -1  -1  -1
     * -1   8  -1
     * -1  -1  -1
     * 
     * Remeber to convert the image to grayscale before running the convolution.
     * 
     * Before returning the image, the corresponding tensor should be clamped in [0,255]
     *  
     * @return returns a new DAISGram containing the modified object
     */
    DAISGram edge() const;

    /**
     * Blend with anoter image
     * 
     * This function generate a new DAISGram which is the composition 
     * of the object and another DAISGram object
     * 
     * The composition follows this convex combination:
     * results = alpha*this + (1-alpha)*rhs 
     * 
     * rhs and this obejct MUST have the same dimensions.
     * 
     * @param rhs The second image involved in the blending
     * @param alpha The parameter of the convex combination  
     * @return returns a new DAISGram containing the blending of the two images.
     */
    DAISGram blend(const DAISGram& rhs, float alpha = 0.5) const;

    /**
     * Green Screen
     * 
     * This function substitutes a pixel with the corresponding one in a background image 
     * if its colors are in the surrounding (+- threshold) of a given color (rgb).
     * 
     * (rgb - threshold) <= pixel <= (rgb + threshold)
     * 
     * 
     * @param bkg The second image used as background
     * @param rgb[] The color to substitute (rgb[0] = RED, rgb[1]=GREEN, rgb[2]=BLUE) 
     * @param threshold[] The threshold to add/remove for each color (threshold[0] = RED, threshold[1]=GREEN, threshold[2]=BLUE) 
     * @return returns a new DAISGram containing the result.
     */
    DAISGram greenscreen(DAISGram& bkg, int rgb[], float threshold[]) const;

    /**
     * Equalize
     * 
     * Stretch the distribution of colors of the image in order to use the full range of intesities.
     * 
     * See https://it.wikipedia.org/wiki/Equalizzazione_dell%27istogramma
     * 
     * @return returns a new DAISGram containing the equalized image.
     */
    DAISGram equalize() const;

    /**
     * Generate Random Image
     * 
     * Generate a random image from nois
     * 
     * @param h height of the image
     * @param w width of the image
     * @param d number of channels
     * @return returns a new DAISGram containing the generated image.
     */
    void generate_random(int h, int w, int d);

    //------------------------ FUNZIONI DI UTILITY AGGIUNTIVE -------------------------------

    /**
     * operatore che permette di confrontare due immagini usando l'operator == tra i relativi tensori 
     * 
     * @param rhs 
     * @return true if equals
     * @return false if not equals
     */
    bool operator==(const DAISGram& rhs) const;

    /**
     * metodo usato per fare il round dei dati del tensore dell'immagine
     * necessario in fase di test per verificare l'effettiva correttezza 
     * del risultato rispetto a quello datoci dal professore
     * 
     * @return DAISGram 
     */
    DAISGram round() const;

    /**
     * salva il tensore dell'immagine su file
     * 
     * @param filename 
     */
    void save_tensor_to_file(string filename) const;

    //------------------------- FUNZIONI AGGIUNTIVE PER LA MANIPOLAZIONE DI IMMAGINI --------------------------

    /**
     * implementazione del filtro sobel (verticale o orizzontale) 
     * utile per il riconiscimento dei contorni
     * 
     * @param horizontal 
     * @return DAISGram 
     */
    DAISGram sobel(bool horizontal = true) const;

    /**
     * implementazione completa del filtro sobel data dall'unione del sobel verticale con quello
     * orizzontale (calcolando la norma pixel per pixel dei due)
     * 
     * @return DAISGram 
     */
    DAISGram full_sobel() const;

    /**
     * specchia l'immagine verticalmente o orizzontalmente
     * 
     * @param vertical 
     * @return DAISGram 
     */
    DAISGram flip(bool vertical = true) const;

    /**
     * inverte i colori dell'immagine 
     * 
     * @return DAISGram 
     */
    DAISGram invert_colours() const;

    /**
     * converte l'encoding dell'immagine da RGB a HSV (hue, saturation, value)
     * usata per equalizzare le immagini colorate
     * 
     * @return DAISGram 
     */
    DAISGram convert_rgb_to_hsv() const;

    /**
     * converte l'encoding dell'immagine da HSV a RGB
     * 
     * @return DAISGram 
     */
    DAISGram convert_hsv_to_rgb() const;

    /**
     * equalizzazione per le immagini colorate, il metodo converte l'immagine in formato HSV
     * e esegue l'equalizzazione sul canale value e infine riconverte l'immagine in formato RGB,
     * ciò permette una più corretta equalizzazione dei colori.
     * 
     * @return DAISGram 
     */
    DAISGram color_equalize() const;

    /**
     * riduce la risoluzione dell'immagine rendendola più "pixellosa"
     * 
     * @param pixels 
     * @return DAISGram 
     */
    DAISGram pixelate(int pixels = 8) const;

    /**
     * riproduce l'immagine usando caratteri ascii e la salva su un file .txt
     * 
     * @param filename 
     */
    void asciiArt(string filename) const;
};

#endif
