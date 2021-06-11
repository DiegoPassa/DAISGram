#include <math.h>
#include <stdio.h>
#include <string.h>

#include "DAISGram.h"
#include "libbmp.h"
#include "tensor.h"

void show_help() {
    printf("*** DAISGram ***\n");
    printf("\targ 1: input file name (img1) \n");
    printf("\targ 2: input file name (img2) \n");
    printf("\targ 3: operazione da effettuare (gray, brighten, blend, sharp, edge, emboss, smooth, warhol, equalize, chromakey) \n");
    printf("\targ 4: output file name\n");
    printf("\targ 5: file da comparare con il risultato ottenuto (se non si vuole comparare nulla inserire '-'\n");
    printf(
        "\targ 6: Diversi significati in funzione dell'operazione:\n"
        "\t\t- [sobel/flip] parametro  per la direzione (orizzontale o verticale)\n"
        "\t\t- [asciiart] file .txt su cui salvare il risultato"
        "\t\t\n");
    printf(
        "\targ 7: Diversi significati in funzione dell'operazione (default 3):\n"
        "\t\t- [smooth]: kernel size \n"
        "\t\t- [brighten]: valore bright per aumentare la luminosità \n"
        "\t\t- [pixelate]: pixel da raggruppare"
        "\t\t\n");
    printf(
        "\targ 8: Diversi significati in funzione dell'operazione (default 1.0):\n"
        "\t\t- [blend] parametro alpha per il blending di due immagini");
    printf("\n");
}

int main(int argc, char* argv[]) {
    char* fn_in_1; /* file 1 */
    char* fn_in_2; /* file 2 */
    char* fn_in_comp;
    char* operation; /* operazione da eseguire */
    char* fn_out;    /* output file */

    int k_size = 3;    /* kernel size */
    float alpha = 1.;  /* alpha della blend */
    char* dir_or_file; /*direzione per sobel e flip o file per asciiart */

    /* variabili di appoggio per le computazioni */
    DAISGram b, c, img, toCompare;

    if (argc < 5) {
        show_help();
        return 0;
    }

    fn_in_1 = argv[1];   /* file 1 */
    fn_in_2 = argv[2];   /* file 2 */
    operation = argv[3]; /* operazione da eseguire */
    fn_out = argv[4];    /* output file */
    fn_in_comp = argv[5];
    dir_or_file = argv[6];

    if (argc > 7) {
        k_size = atoi(argv[7]);
    }

    if (argc > 8) {
        alpha = atof(argv[8]);
    }

    cout << "Prima immagine processata: " << fn_in_1 << "\n";
    cout << "Operazione: " << operation << "\n";
    cout << "Seconda immagine processata: " << fn_in_2 << "\n";
    cout << "Parametri: " << k_size << ", " << alpha << "\n";
    cout << "Immagine da confrontare: " << fn_in_comp << "\n";
    cout << "direzione/file .txt per asciiart: " << dir_or_file << "\n";
    cout << "File di output: " << fn_out << "\n\n";

    b.load_image(fn_in_1); /* leggi il file di input */

    if (strcmp(operation, "brighten") == 0) {
        img = b.brighten(k_size); /* aumenta la luminosità */
    } else if (strcmp(operation, "blend") == 0) {
        c.load_image(fn_in_2);
        img = b.blend(c, alpha); /* effettua il blending di due immagini */
    } else if (strcmp(operation, "gray") == 0) {
        img = b.grayscale();
    } else if (strcmp(operation, "equalize") == 0) {
        img = b.equalize();
    } else if (strcmp(operation, "chromakey") == 0) {
        c.load_image(fn_in_2);
        int r_, g_, b_;
        float thr, thg, thb;
        cout << "Enter green-screen parameters (int RGB[3]) and (float RGB Threshold[3])" << endl;
        cin >> r_ >> g_ >> b_ >> thr >> thg >> thb;
        int rgb[3] = {r_, g_, b_};
        float th[3] = {thr, thg, thb};
        img = b.greenscreen(c, rgb, th);
    } else if (strcmp(operation, "sharp") == 0) {
        img = b.sharpen();
    } else if (strcmp(operation, "edge") == 0) {
        img = b.edge();
    } else if (strcmp(operation, "emboss") == 0) {
        img = b.emboss();
    } else if (strcmp(operation, "smooth") == 0) {
        img = b.smooth(k_size);
    } else if (strcmp(operation, "warhol") == 0) {
        img = b.warhol();
    } else if (strcmp(operation, "sobel") == 0) {
        bool dir = true;

        if (strcmp(dir_or_file, "verticale") == 0) dir = false;

        img = b.sobel(dir);
    } else if (strcmp(operation, "full_sobel") == 0) {
        img = b.full_sobel();
    } else if (strcmp(operation, "pixelate") == 0) {
        img = b.pixelate(k_size);
    } else if (strcmp(operation, "color_eq") == 0) {
        img = b.color_equalize();
    } else if (strcmp(operation, "flip") == 0) {
        bool dir = true;

        if (strcmp(dir_or_file, "orizzontale") == 0) dir = false;

        img = b.flip(dir);
    } else if (strcmp(operation, "invert_col") == 0) {
        img = b.invert_colours();
    } else if (strcmp(operation, "asciiart") == 0) {
        string file = dir_or_file;
        b.asciiArt(file);
    } else {
        throw(unknown_operation());
    }

    //img.save_tensor_to_file("./img.txt");
    //toCompare.save_tensor_to_file("./comp.txt");

    if (strcmp(operation, "asciiart") != 0) {
        if (strcmp(fn_in_comp, "-") != 0) {
            try {
                toCompare.load_image(fn_in_comp);  //immagine corretta da confrontare
                img = img.round();

                if (img == toCompare)
                    cout << "--------------------------RISULTATO CORRETTO-------------------------------\n\n";
                else
                    cout << "---------------------------RISULTATO ERRATO--------------------------------\n\n";
            } catch (const dimension_mismatch& e) {
                std::cerr << "Impossibile confrontare le immagini a causa della diverse dimensione delle due." << '\n';
            }
        }

        img.save_image(fn_out);
    }

    return 0; /* ciao a tutti!*/
}
