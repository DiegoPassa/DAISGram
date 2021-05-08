#include "DAISGram.h"

#include <iostream>
#include <string>

#include "dais_exc.h"
#include "libbmp.h"
#include "tensor.h"

using namespace std;

/**
 * Load a bitmap from file
 *
 * @param filename String containing the path of the file
 */
void DAISGram::load_image(string filename) {
    BmpImg img = BmpImg();

    img.read(filename.c_str());

    const int h = img.get_height();
    const int w = img.get_width();

    data = Tensor(h, w, 3, 0.0);

    for (int i = 0; i < img.get_height(); i++) {
        for (int j = 0; j < img.get_width(); j++) {
            data(i, j, 0) = (float)img.red_at(j, i);
            data(i, j, 1) = (float)img.green_at(j, i);
            data(i, j, 2) = (float)img.blue_at(j, i);
        }
    }
}

/**
 * Save a DAISGram object to a bitmap file.
 * 
 * Data is clamped to 0,255 before saving it.
 *
 * @param filename String containing the path where to store the image.
 */
void DAISGram::save_image(string filename) {
    data.clamp(0, 255);

    BmpImg img = BmpImg(getCols(), getRows());

    img.init(getCols(), getRows());

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            img.set_pixel(j, i, (unsigned char)data(i, j, 0), (unsigned char)data(i, j, 1), (unsigned char)data(i, j, 2));
        }
    }

    img.write(filename);
}

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
void DAISGram::generate_random(int h, int w, int d) {
    data = Tensor(h, w, d, 0.0);
    data.init_random(128, 50);
    data.rescale(255);
}

/**
 * Get rows
 *
 * @return returns the number of rows in the image
 */
int DAISGram::getRows() {
    return data.rows();
}

/**
 * Get columns
 *
 * @return returns the number of columns in the image
 */
int DAISGram::getCols() {
    return data.cols();
}

/**
 * Get depth
 *
 * @return returns the number of channels in the image
 */
int DAISGram::getDepth() {
    return data.depth();
}

DAISGram DAISGram::brighten(float bright) {
    DAISGram brightened;
    brightened.data = data;
    for (int i = 0; i < data.rows(); i++)
        for (int j = 0; j < data.cols(); j++)
            for (int k = 0; k < data.depth(); k++)
                brightened.data(i, j, k) += bright;

    brightened.data.clamp(0, 255);

    return brightened;
}

DAISGram DAISGram::grayscale() {
    DAISGram gray;
    gray.data = data;
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            int media = 0;
            for (int k = 0; k < data.depth(); k++) {
                media += data(i, j, k);
            }
            media /= data.depth();
            for (int k = 0; k < data.depth(); k++) {
                gray.data(i, j, k) = media;
            }
        }
    }
    return gray;
}

DAISGram DAISGram::warhol() {
    DAISGram result;

    Tensor red{data.subset(0, data.rows(), 0, data.cols(), 0, 1)};
    Tensor green{data.subset(0, data.rows(), 0, data.cols(), 1, 2)};
    Tensor blue{data.subset(0, data.rows(), 0, data.cols(), 2, 3)};
    Tensor new_t, bottom;

    //parte superiore
    result.data = data;
    new_t = green.concat(red, 2).concat(blue, 2);
    result.data = result.data.concat(new_t, 1);

    //parte inferiore
    new_t = red.concat(blue, 2).concat(green, 2);
    bottom = new_t;
    new_t = blue.concat(green, 2).concat(red, 2);
    bottom = bottom.concat(new_t, 1);

    result.data = result.data.concat(bottom, 0);
    return result;
}

DAISGram DAISGram::sharpen() {
    Tensor filter;
    float f[3 * 3] = {0, -1, 0, -1, 5, -1, 0, -1, 0};

    filter.init_filter(f, 3, 3);
    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage.data = data.convolve(filter);
    
    return newImage;
}

DAISGram DAISGram::blend(const DAISGram& rhs, float alpha) {
    if (alpha < 0 || alpha > 1)
        throw(invalid_parameter());

    DAISGram new_d;
    new_d.data = data * alpha + rhs.data * (1 - alpha);

    return new_d;
}

DAISGram DAISGram::greenscreen(DAISGram& bkg, int rgb[], float threshold[]) {
    if (getRows() != bkg.getRows() || getCols() != bkg.getCols() || getDepth() != bkg.getDepth()) {
        throw(dimension_mismatch());
    }

    DAISGram newImage;
    newImage.data = data;
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            bool flag = true;
            for (int k = 0; k < data.depth(); k++) {
                if (!flag || !(data(i, j, k) >= (rgb[k] - threshold[k]) && data(i, j, k) <= (rgb[k] + threshold[k]))) {
                    flag = false;
                }
            }
            if (flag) {
                for (int k = 0; k < data.depth(); k++) {
                    newImage.data(i, j, k) = bkg.data(i, j, k);
                }
            }
        }
    }
    return newImage;
}

DAISGram DAISGram::equalize() {
    DAISGram equalized;
    equalized.data = data;

    for (int i = 0; i < equalized.data.depth(); i++) {
        int istogram[256] = {};
        for (int j = 0; j < equalized.data.rows(); j++) {
            for (int k = 0; k < equalized.data.cols(); k++) {
                int t = equalized.data(j, k, i);
                ++istogram[t];
            }
        }

        int cdf[256] = {};
        int cdf_min = 0;
        int t = 0;

        for (int j = 0; j < 256; j++) {
            while (istogram[j] == 0)
                ++j;
            if (t == 0)
                cdf_min = istogram[j];
            cdf[j] = istogram[j] + t;
            t = cdf[j];
        }

        for (int j = 0; j < equalized.data.rows(); j++) {
            for (int k = 0; k < equalized.data.cols(); k++) {
                int v = equalized.data(j, k, i);
                equalized.data(j, k, i) = (cdf[v] - cdf_min) * 255 / (equalized.data.rows() * equalized.data.cols() - 1);
            }
        }
    }

    return equalized;
}

DAISGram DAISGram::edge() {
    Tensor filter;
    float f[3 * 3] = {-1, -1, -1, -1, 8, -1, -1, -1, -1};

    filter.init_filter(f, 3, 3);
    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage = grayscale();
    newImage.data = data.convolve(filter);

    return newImage;
}