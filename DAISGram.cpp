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

    for (int i = 0; i < img.get_height(); ++i) {
        for (int j = 0; j < img.get_width(); ++j) {
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
int DAISGram::getRows() const {
    return data.rows();
}

/**
 * Get columns
 *
 * @return returns the number of columns in the image
 */
int DAISGram::getCols() const {
    return data.cols();
}

/**
 * Get depth
 *
 * @return returns the number of channels in the image
 */
int DAISGram::getDepth() const {
    return data.depth();
}

DAISGram DAISGram::brighten(float bright) const {
    DAISGram brightened;
    brightened.data = data;
    for (int i = 0; i < data.rows(); i++)
        for (int j = 0; j < data.cols(); j++)
            for (int k = 0; k < data.depth(); k++)
                brightened.data(i, j, k) += bright;

    brightened.data.clamp(0, 255);

    return brightened;
}

DAISGram DAISGram::grayscale() const {
    DAISGram gray;
    gray.data = data;
    for (int i = 0; i < data.rows(); i++) {
        for (int j = 0; j < data.cols(); j++) {
            float media = 0;
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

DAISGram DAISGram::warhol() const {
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

DAISGram DAISGram::sharpen() const {
    Tensor filter;
    float f[3 * 3] = {0, -1, 0,
                      -1, 5, -1,
                      0, -1, 0};

    filter.init_filter(f, 3, 3);
    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage.data = data.convolve(filter);

    newImage.data.clamp(0, 255);
    //newImage.data.rescale(255);

    return newImage;
}

DAISGram DAISGram::emboss() const {
    Tensor filter;
    float f[3 * 3] = {-2, -1, 0,
                      -1, 1, 1,
                      0, 1, 2};

    filter.init_filter(f, 3, 3);
    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage.data = data.convolve(filter);

    newImage.data.clamp(0, 255);
    //newImage.data.rescale(255);

    return newImage;
}

DAISGram DAISGram::edge() const {
    Tensor filter;
    float f[3 * 3] = {-1, -1, -1,
                      -1, 8, -1,
                      -1, -1, -1};

    filter.init_filter(f, 3, 3);
    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage = grayscale();
    newImage.data = newImage.data.convolve(filter);

    newImage.data.clamp(0, 255);
    //newImage.data.rescale(255);

    return newImage;
}

DAISGram DAISGram::smooth(int h) const {
    float c = 1.f / (h * h);
    Tensor filter{h, h, 3, c};

    DAISGram newImage;
    newImage.data = data.convolve(filter);

    newImage.data.clamp(0, 255);
    //newImage.data.rescale(255);

    return newImage;
}

DAISGram DAISGram::blend(const DAISGram& rhs, float alpha) const {
    if (alpha < 0 || alpha > 1)
        throw(invalid_parameter());

    DAISGram new_d;
    new_d.data = data * alpha + rhs.data * (1 - alpha);

    return new_d;
}

DAISGram DAISGram::greenscreen(DAISGram& bkg, int rgb[], float threshold[]) const {
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

DAISGram DAISGram::equalize() const {
    DAISGram equalized;
    equalized.data = data;

    for (int i = 0; i < equalized.data.depth(); i++) {
        vector<int> istogram(256, 0);
        for (int j = 0; j < equalized.data.rows(); j++) {
            for (int k = 0; k < equalized.data.cols(); k++) {
                int t = equalized.data(j, k, i);
                ++istogram[t];
            }
        }

        vector<int> cdf(256, 0);
        int cdf_min = 0;
        int t = 0;

        for (int j = 0; j < 256; j++) {
            if (istogram[j] != 0) {
                if (t == 0)
                    cdf_min = istogram[j];

                cdf[j] = istogram[j] + t;
                t = cdf[j];
            }
        }

        for (int j = 0; j < equalized.data.rows(); j++) {
            for (int k = 0; k < equalized.data.cols(); k++) {
                int v = equalized.data(j, k, i);
                equalized.data(j, k, i) = (cdf[v] - cdf_min) * 255 / (equalized.data.rows() * equalized.data.cols() - cdf_min);
            }
        }
    }

    return equalized;
}

bool DAISGram::operator==(const DAISGram& rhs) const {
    bool equals = false;

    if (data == rhs.data) {
        equals = true;
    }

    return equals;
}

DAISGram DAISGram::round() const {
    DAISGram rounded{*this};

    for (int k = 0; k < getDepth(); k++) {
        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                rounded.data(i, j, k) = (int)rounded.data(i, j, k);
            }
        }
    }

    return rounded;
}

void DAISGram::save_tensor_to_file(string filename) const {
    data.write_file(filename);
}

DAISGram DAISGram::sobel(bool horizontal) const {
    Tensor filter;

    if (horizontal) {
        float f[3 * 3] = {-1, -2, -1,
                          0, 0, 0,
                          1, 2, 1};
        filter.init_filter(f, 3, 3);
    } else {
        float f[3 * 3] = {-1, 0, 1,
                          -2, 0, 2,
                          -1, 0, 1};
        filter.init_filter(f, 3, 3);
    }

    filter = filter.concat(filter, 2).concat(filter, 2);

    DAISGram newImage;
    newImage = grayscale();
    newImage.data = newImage.data.convolve(filter);

    newImage.data.clamp(0, 255);
    //newImage.data.rescale(255);

    return newImage;
}

DAISGram DAISGram::full_sobel() const {
    DAISGram sh, sv, s;
    s.data.init(getRows(), getCols(), getDepth());
    sh = sobel();
    sv = sobel(false);

    for (int k = 0; k < getDepth(); k++) {
        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                s.data(i, j, k) = sqrtf(sh.data(i, j, k) * sh.data(i, j, k) + sv.data(i, j, k) * sv.data(i, j, k));
            }
        }
    }

    s.data.clamp(0, 255);
    //s.data.rescale(255);

    return s;
}

DAISGram DAISGram::flip(bool vertical) const {
    DAISGram flipped;
    flipped.data.init(getRows(), getCols(), getDepth());
    for (int k = 0; k < getDepth(); k++) {
        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                if (vertical)
                    flipped.data(getRows() - i - 1, j, k) = data(i, j, k);
                else
                    flipped.data(i, getCols() - j - 1, k) = data(i, j, k);
            }
        }
    }
    return flipped;
}

DAISGram DAISGram::invert_colours() const {
    DAISGram inverted;
    inverted.data = data;
    for (int k = 0; k < getDepth(); k++) {
        for (int i = 0; i < getRows(); i++) {
            for (int j = 0; j < getCols(); j++) {
                inverted.data(i, j, k) = 255 - inverted.data(i, j, k);
            }
        }
    }
    return inverted;
}

DAISGram DAISGram::convert_rgb_to_hsv() const {
    DAISGram image_converted;
    image_converted.data.init(getRows(), getCols(), getDepth());

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            float r = data(i, j, 0) / 255.f;
            float g = data(i, j, 1) / 255.f;
            float b = data(i, j, 2) / 255.f;

            float max_rgb = max(r, max(g, b));
            float min_rgb = min(r, min(g, b));
            float diff = max_rgb - min_rgb;

            //hue
            if (max_rgb == min_rgb) {
                image_converted.data(i, j, 0) = 0;
            } else if (max_rgb == r) {
                image_converted.data(i, j, 0) = fmod((60 * ((g - b) / diff) + 360), 360);
            } else if (max_rgb == g) {
                image_converted.data(i, j, 0) = fmod((60 * ((b - r) / diff) + 120), 360);
            } else if (max_rgb == b) {
                image_converted.data(i, j, 0) = fmod((60 * ((r - g) / diff) + 240), 360);
            }

            //saturation
            if (max_rgb == 0) {
                image_converted.data(i, j, 1) = 0;
            } else {
                image_converted.data(i, j, 1) = (diff / max_rgb) * 100;
            }

            //value
            image_converted.data(i, j, 2) = max_rgb * 100;
        }
    }

    return image_converted;
}

DAISGram DAISGram::convert_hsv_to_rgb() const {
    DAISGram image_converted;
    image_converted.data.init(getRows(), getCols(), getDepth());

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++) {
            float H = data(i, j, 0);
            float S = data(i, j, 1);
            float V = data(i, j, 2);

            if (H > 360 || H < 0 || S > 100 || S < 0 || V > 100 || V < 0) {
                throw(index_out_of_bound());
            }

            float s = S / 100;
            float v = V / 100;
            float C = s * v;
            float X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
            float m = v - C;

            float r, g, b;
            if (H >= 0 && H < 60) {
                r = C, g = X, b = 0;
            } else if (H >= 60 && H < 120) {
                r = X, g = C, b = 0;
            } else if (H >= 120 && H < 180) {
                r = 0, g = C, b = X;
            } else if (H >= 180 && H < 240) {
                r = 0, g = X, b = C;
            } else if (H >= 240 && H < 300) {
                r = X, g = 0, b = C;
            } else {
                r = C, g = 0, b = X;
            }

            image_converted.data(i, j, 0) = (r + m) * 255;
            image_converted.data(i, j, 1) = (g + m) * 255;
            image_converted.data(i, j, 2) = (b + m) * 255;
        }
    }

    return image_converted;
}

DAISGram DAISGram::color_equalize() const {
    DAISGram imageEqualized;
    imageEqualized = convert_rgb_to_hsv();

    vector<int> istogram(101, 0);

    for (int j = 0; j < imageEqualized.data.rows(); j++) {
        for (int k = 0; k < imageEqualized.data.cols(); k++) {
            int t = imageEqualized.data(j, k, 2);
            ++istogram[t];
        }
    }

    vector<int> cdf(istogram.size(), 0);
    int cdf_min = 0;
    int t = 0;

    for (size_t j = 0; j < cdf.size(); j++) {
        if (istogram[j] != 0) {
            if (t == 0)
                cdf_min = istogram[j];

            cdf[j] = istogram[j] + t;
            t = cdf[j];
        }
    }

    for (int j = 0; j < imageEqualized.data.rows(); j++) {
        for (int k = 0; k < imageEqualized.data.cols(); k++) {
            int v = imageEqualized.data(j, k, 2);
            imageEqualized.data(j, k, 2) = (cdf[v] - cdf_min) * 100 / (imageEqualized.data.rows() * imageEqualized.data.cols() - cdf_min);
        }
    }

    return imageEqualized.convert_hsv_to_rgb();
}

DAISGram DAISGram::pixelate(int pixels) const {
    DAISGram pixelated;
    pixelated.data = data;

    for (int i = 0; i < getRows(); i += pixels)
        for (int j = 0; j < getCols(); j += pixels)
            for (int k = 0; k < getDepth(); k++)
                for (int pr = 0; pr < pixels && pr + i < getRows(); pr++)
                    for (int pc = 0; pc < pixels && pc + j < getCols(); pc++)
                        pixelated.data(pr + i, pc + j, k) = data(i, j, k);

    return pixelated;
}

void DAISGram::asciiArt(string filename) const {
    ofstream ost{filename};
    char map[11] = " .,:;ox%#@";

    DAISGram gray;
    gray.data = data;
    gray = gray.grayscale();

    for (int i = 0; i < getRows(); i++) {
        for (int j = 0; j < getCols(); j++)
            ost << map[(255 - (int)gray.data(i, j, 0)) * 10 / 256];
        ost << endl;
    }
}