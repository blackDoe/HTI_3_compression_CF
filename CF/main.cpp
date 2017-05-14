#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fichiers.h"
#include "matrix.h"
#include "pred.h"
#include "dct.h"


int SaveIntImage_pgm(char *nom, int **im, int Height, int Width) {

    int i, j;
    unsigned char **ima_uc = alocamuc(Height, Width);
    int maxim = 0;

    for (i = 0; i < Height; i++)
        for (j = 0; j < Width; j++)
            if (abs(im[i][j]) > maxim)
                maxim = abs(im[i][j]);

    //  fprintf(stderr, "\nnom = %s  maxim = %d\n", nom, maxim);

    for (i = 0; i < Height; i++)
        for (j = 0; j < Width; j++)
            ima_uc[i][j] = (unsigned char)(abs(im[i][j]) * 255 / maxim);

    ecriture_pgm(nom, ima_uc, Width, Height);

    dalocuc(ima_uc, Height);
    return EXIT_FAILURE;
}

int SaveIntImage_pgm_tronc(char *nom, int **im, int Height, int Width) {

    int i, j;

    unsigned char **ima_uc = alocamuc(Height, Width);
    for (i = 0; i < Height; i++)
        for (j = 0; j < Width; j++)
        {
            if (abs(im[i][j]) >= 255.0)
                ima_uc[i][j] = 255;
            else
                ima_uc[i][j] = (unsigned char)(abs(im[i][j]));
        }

    ecriture_pgm(nom, ima_uc, Width, Height);

    dalocuc(ima_uc, Height);
    return EXIT_FAILURE;
}


void My_dct2dim(double **xd, double **tdct1, int H, int W) {

    // Valeurs utiles : on ne les calcule qu'une fois car long à écrire
    double sqrth = sqrt((double) (H));
    double sqrtw = sqrt((double) (W));
    double sqrt2 = sqrt(2.);

    int i, j;   // Pour la définition des matrices utiles

    double **matDct1 = alocamd(H, W);
    for (j = 0; j < W; j++)
        matDct1[0][j] = 1 / sqrth;
    for (i = 1; i < H; i++)
        for (j = 0; j < W; j++)
            matDct1[i][j] = (sqrt2 / sqrth) * cos(M_PI * (2*j + 1) * i / (2*H));

    double **matDct2 = alocamd(W, H);
    for (i = 0; i < H; i++) {
        matDct2[i][0] = 1 / sqrtw;
        for (j = 1; j < W; j++)
            matDct2[i][j] = (sqrt2 / sqrtw) * cos(M_PI * (2*i + 1) * j / (2*W));
    }

    // Image transition
    double **tdct0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            tdct0[i][j] = tdct1[i][j] = 0.0;

    int a, b;

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (i = 0; i < H; i++)
                tdct0[a][b] += matDct1[a][i] * xd[i][b];

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (j = 0; j < W; j++)
                tdct1[a][b] += tdct0[a][j] * matDct2[j][b];

}

void My_dct2dim_inv(double **tdct1, double **xrecd1, int H, int W) {

    // Valeurs utiles : on ne les calcule qu'une fois car long à écrire
    double sqrth = sqrt((double)(H));
    double sqrtw = sqrt((double)(W));
    double sqrt2 = sqrt(2.);

    int i, j;   // Pour la définition des matrices utiles

    // Transposée de matDct1
    double **matDct1Inv = alocamd(W, H);
    for (j = 0; j < W; j++)
        matDct1Inv[j][0] = 1 / sqrth;
    for (i = 1 ; i < H; i++)
        for(j = 0; j < W; j++)
            matDct1Inv[j][i] = (sqrt2/sqrth) * cos(M_PI * (2*j + 1) * i / (2*H));

    // Transposée de matDct2
    double **matDct2Inv = alocamd(H, W);
    for (i = 0; i < H; i++) {
        matDct2Inv[0][i] = 1 / sqrtw;
        for (j = 1; j < W; j++)
            matDct2Inv[j][i] = (sqrt2/sqrtw) * cos(M_PI * (2*i + 1) * j / (2*W));
    }

    double **xrecd0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xrecd0[i][j] = xrecd1[i][j] = 0.0;

    int a, b;

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (i = 0; i < H; i++)
                xrecd0[a][b] += matDct1Inv[a][i] * tdct1[i][b];

    for (a = 0; a < H; a++)
        for (b = 0; b <W; b++)
            for (j = 0; j < W; j++)
                xrecd1[a][b] += xrecd0[a][j] * matDct2Inv[j][b];

}

void My_dct2dim_bloc(double **xd, double **tdct1, int H, int W, int Bx, int By, int step) {

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtBx = sqrt((double)(Bx));
    double sqrtBy = sqrt((double)(By));
    double sqrt2 = sqrt(2.);

    int i, j;   // Pour la définition des matrices utiles

    double **matDct1= alocamd(Bx, By);
    for (j = 0; j < By; j++)
        matDct1[0][j] = 1 / sqrtBx;
    for (i = 1; i < Bx; i++)
        for(j = 0; j < By; j++)
            matDct1[i][j] = (sqrt2/sqrtBx) * cos(M_PI * (2*j + 1) * i / (2*Bx));

    double **matDct2 = alocamd(By, Bx);
    for (i = 0; i < Bx; i++) {
        matDct2[i][0] = 1/sqrtBy;
        for ( j = 1; j < By; j++)
            matDct2[i][j] = (sqrt2/sqrtBy) * cos(M_PI * (2*i + 1) * j / (2*By));
    }

    // Image transition test 1
    double **tdct0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            tdct0[i][j] = tdct1[i][j] = 0.0;

    int nbx, nby;   // Indice du bloc courant
    int a, b;       // Indice du pixel courant dans le bloc

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (i = 0; i < Bx; i++)
                        tdct0[a+(nbx*Bx)][b+(nby*By)] += matDct1[a][i] * xd[i+(nbx*Bx)][b+(nby*By)];

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (j = 0; j < Bx; j++)
                        tdct1[a+(nbx*Bx)][b+(nby*By)] += tdct0[a+(nbx*Bx)][j+(nby*By)] * matDct2[j][b];

}

void My_dct2dim_bloc_inv(double **tdct1, double **xrecd1, int H, int W, int Bx, int By) {

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtBx = sqrt((double)(Bx));
    double sqrtBy = sqrt((double)(By));
    double sqrt2 = sqrt(2.);

    int i, j;    // Pour la définition des matrices utiles

    // Trasnposée de matDct1 par blocs
    double **matDct1Inv = alocamd(Bx, By);
    for (j = 0; j < By; j++)
        matDct1Inv[j][0] = 1 / sqrtBx;
    for (i = 1; i < Bx; i++)
        for(j = 0; j < By; j++)
            matDct1Inv[j][i] = (sqrt2/sqrtBx) * cos(M_PI * (2*j + 1) * i / (2*Bx));

    // Transposée de matDct2 par bloc
    double **matDct2Inv = alocamd(By, Bx);
    for (i = 0; i < Bx; i++) {
        matDct2Inv[0][i] = 1/sqrtBy;
        for (j = 1; j < By; j++)
            matDct2Inv[j][i] = (sqrt2/sqrtBy) * cos(M_PI * (2*i + 1) * j / (2*By));
    }

    // Image transitoire
    double **xrecd0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xrecd0[i][j] = xrecd1[i][j] = 0.0;

    int nbx, nby;       // Indice du bloc courant
    int a, b;           // Indice du pixel courant dans le bloc

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (i = 0; i < Bx; i++)
                        xrecd0[a+(nbx*Bx)][b+(nby*By)] += matDct1Inv[a][i] * tdct1[i+(nbx*Bx)][b+(nby*By)];

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (j = 0; j < Bx; j++)
                        xrecd1[a+(nbx*Bx)][b+(nby*By)] += xrecd0[a+(nbx*Bx)][j+(nby*By)] * matDct2Inv[j][b];

}


void My_dst2dim(double **xd, double **tdst1, int H, int W) {

    // Valeurs utiles : on ne les calcule qu'une fois car long à écrire
    double sqrtfh = sqrt(2./(double)(H+1));
    double sqrtfw = sqrt(2./(double)(W+1));

    int i, j;   // Pour la définition des matrices utiles

    double **matDst1 = alocamd(H, W);
    for (i = 0 ; i < H; i++)
        for (j = 0; j < W; j++)
            matDst1[i][j] = sqrtfh * sin(M_PI * (j+1) * (i+1) / (H+1));

    double **matDst2 = alocamd(W, H);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            matDst2[i][j] = sqrtfw * sin(M_PI * (i+1) * (j+1) / (W+1));

    // Image transition
    double **tdst0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            tdst0[i][j] = tdst1[i][j] = 0.0;

    int a, b;

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (i = 0; i < H; i++)
                tdst0[a][b] += matDst1[a][i] * xd[i][b];

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (j = 0; j < W; j++)
                tdst1[a][b] += tdst0[a][j] * matDst2[j][b];

}

void My_dst2dim_inv(double **tdst1, double **xrecd1, int H, int W) {

    // Valeurs utiles : on ne les calcule qu'une fois car long à écrire
    double sqrtfh = sqrt(2./(double)(H+1));
    double sqrtfw = sqrt(2./(double)(W+1));

    int i, j;   // Pour la définition des matrices utiles

    // Transposée de matDst1
    double **matDst1Inv = alocamd(W, H);
    for (i = 0 ; i < H; i++)
        for(j = 0; j < W; j++)
            matDst1Inv[j][i] = sqrtfh * sin(M_PI * (j+1) * (i+1) / (H+1));

    // Transposée de matDst2
    double **matDst2Inv = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            matDst2Inv[j][i] = sqrtfw * sin(M_PI * (i+1) * (j+1) / (W+1));

    double **xrecd0 = alocamd(H, W);
    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++)
            xrecd0[i][j] = xrecd1[i][j] = 0.0;

    int a, b;

    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (i = 0; i < H; i++)
                xrecd0[a][b] += matDst1Inv[a][i] * tdst1[i][b];

    for (a = 0; a < H; a++)
        for (b = 0; b <W; b++)
            for (j = 0; j < W; j++)
                xrecd1[a][b] += xrecd0[a][j] * matDst2Inv[j][b];

}

void My_dst2dim_bloc(double **xd, double **tdst1, int H, int W, int Bx, int By, int step) {

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtfx = sqrt(2./(double)(Bx+1));
    double sqrtfy = sqrt(2./(double)(By+1));

    int i, j;   // Pour la définition des matrices utiles

    double **matDst1= alocamd(Bx, By);
    for (i = 0; i < Bx; i++)
        for(j = 0; j < By; j++)
            matDst1[i][j] = sqrtfx * sin(M_PI * (j+1) * (i+1) / (Bx+1));

    double **matDst2 = alocamd(By, Bx);
    for (i = 0; i < Bx; i++)
        for (j = 0; j < By; j++)
            matDst2[i][j] = sqrtfy * sin(M_PI * (i+1) * (j+1) / (By+1));

    // Image transition test 1
    double **tdst0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            tdst0[i][j] = tdst1[i][j] = 0.0;

    int nbx, nby;   // Indice du bloc courant
    int a, b;       // Indice du pixel courant dans le bloc

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (i = 0; i < Bx; i++)
                        tdst0[a+(nbx*Bx)][b+(nby*By)] += matDst1[a][i] * xd[i+(nbx*Bx)][b+(nby*By)];

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (j = 0; j < Bx; j++)
                        tdst1[a+(nbx*Bx)][b+(nby*By)] += tdst0[a+(nbx*Bx)][j+(nby*By)] * matDst2[j][b];

}

void My_dst2dim_bloc_inv(double **tdst1, double **xrecd1, int H, int W, int Bx, int By) {

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtfx = sqrt(2./(double)(Bx+1));
    double sqrtfy = sqrt(2./(double)(By+1));

    int i, j;    // Pour la définition des matrices utiles

    // Transposée de matDst1 par blocs
    double **matDst1Inv = alocamd(Bx, By);
    for (i = 0; i < Bx; i++)
        for(j = 0; j < By; j++)
            matDst1Inv[j][i] = sqrtfx * sin(M_PI * (j+1) * (i+1) / (Bx+1));

    // Transposée de matDst2 par bloc
    double **matDst2Inv = alocamd(By, Bx);
    for (i = 0; i < Bx; i++)
        for (j = 0; j < By; j++)
            matDst2Inv[j][i] = sqrtfy * sin(M_PI * (i+1) * (j+1) / (By+1));

    // Image transitoire
    double **xrecd0 = alocamd(H, W);
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xrecd0[i][j] = xrecd1[i][j] = 0.0;

    int nbx, nby;       // Indice du bloc courant
    int a, b;           // Indice du pixel courant dans le bloc

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (i = 0; i < Bx; i++)
                        xrecd0[a+(nbx*Bx)][b+(nby*By)] += matDst1Inv[a][i] * tdst1[i+(nbx*Bx)][b+(nby*By)];

    for (nbx = 0; nbx < (H/Bx); nbx++)
        for (nby = 0; nby < (W/By); nby++)
            for (a = 0; a < Bx; a++)
                for (b = 0; b <By; b++)
                    for (j = 0; j < Bx; j++)
                        xrecd1[a+(nbx*Bx)][b+(nby*By)] += xrecd0[a+(nbx*Bx)][j+(nby*By)] * matDst2Inv[j][b];

}


void quantification(double **tdt, int **err, int H, int W, int step) {

    int i, j;

    for (i = 1; i < H; i++)
        for (j = 1; j < W; j++) {
            err[i][j] = quantiz(tdt[i][j], step);
            tdt[i][j] = (double)(err[i][j]);
        }

    for (j = 0; j < W; j++) {
        err[0][j] = quantiz(tdt[0][j], 1);
        tdt[0][j] = (double)(err[0][j]);
    }

    for (i = 0; i < H; i++) {
        err[i][0] = quantiz(tdt[i][0], 1);
        tdt[i][0] = (double)(err[i][0]);
    }

}

void finale(double **xrecd, unsigned char **xrec, int H, int W) {

    int i, j;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            if (xrecd[i][j] < 0.0)
                xrec[i][j] = 0;
            else if (xrecd[i][j] > 255.0)
                xrec[i][j] = 255;
            else
                xrec[i][j] = (unsigned char)(xrecd[i][j]);
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////        /////////////////////////////////////////////////////////////
////////////////////////////////////////////////     MAIN     //////////////////////////////////////////////////////////
///////////////////////////////////////////////////        /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main (int argc, char *argv[]) {

    char nom[200], nom_out[200], nom_err[300];
    char nom_out_My[200], nom_err_My[300];
    int W, H;   // Dimensions de l'image: H = Height, W = Width
    int i, j;

    if (argc != 3) {
        fprintf(stderr, "Utilisation: %s <nom_image_pgm> <pas_quantification> \n", argv[0]);
        exit(0);
    }

    strcpy(nom, argv[1]);
    strcpy(nom_out, argv[1]); strcat(nom_out, ".out");
    strcpy(nom_err, argv[1]); strcat(nom_err, ".err");
    strcpy(nom_out_My, argv[1]); strcat(nom_out_My, "_My.out");
    strcpy(nom_err_My, argv[1]); strcat(nom_err_My, "_My.err");

    int step = atoi(argv[2]);

    lecture_dim(nom, &W, &H);

    fprintf(stderr, "Width = %d  Height = %d\n", W, H);

    unsigned char **x = alocamuc(H, W);
    unsigned char **xrec = alocamuc(H, W);  // Image finale en utilisant la dct donnée
    unsigned char **xrec1 = alocamuc(H, W); // Notre image finale

    double **xrecd = alocamd(H, W);     // Image reconstruite
    double **xrecd1 = alocamd(H, W);    // Notre image reconstruite

    lecture_pgm(nom, x);    // x contient l'image initiale

    int **err = alocami(H, W);      // Allocation de la matrice des erreurs de prediction
    int **err1 = alocami(H, W);     // Allocation de notre matrice des erreurs de prediction

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err[i][j] = (int)(x[i][j]);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err1[i][j] = (int)(x[i][j]);

    fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));

    /////////////////////////////////////////   CODE DCT  //////////////////////////////////////////////////////////////
    /*
    double **tdct = alocamd(H, W);  // image transformee
    double **tdct1 = alocamd(H, W);  // Mon image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xd[i][j] = (double)(x[i][j]);


    // Application de la compression dct
    dct2dim(xd, tdct, H, W);
    fprintf(stderr, "OK DCT directe\n");

    My_dct2dim(xd, tdct1, H, W);
    fprintf(stderr, "\nMa DCT directe effectuée\n");


    // Quantification des coefficients
    quantification(tdct, err, H, W, step);
    quantification(tdct1, err1, H, W, step);

    // Calcul de l'entropie de la compression
    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie en utilisant la dct donnée = %g [bits/pixel]\n", entro);
    double Myentro = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", Myentro);


    // Application de la dct inverse (décompression)
    dct2dim_inv(tdct, xrecd, H, W);
    My_dct2dim_inv(tdct1, xrecd1, H, W);
    fprintf(stderr, "\nMa DCT inverse effectuée\n");


    // écriture de la matrice image finale
    finale(xrecd, xrec, H, W);
    SaveIntImage_pgm_tronc(nom_err, err, H, W);
    ecriture_pgm(nom_out, xrec, W, H);

    finale(xrecd1, xrec1, H, W);
    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);


    // libération de mémoire
    dalocd(xrecd, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    dalocd(xrecd1, H);
    dalocd(tdct1, H);
    */
    ////////////////////////////////////////    FIN CODE DCT     ///////////////////////////////////////////////////////

    //////////////////////////////////////    DCT PAR BLOCS    /////////////////////////////////////////////////////////
    /*
    double **tdct = alocamd(H, W);  // image transformee
    double **tdct1 = alocamd(H, W); // Mon image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
      for (j = 0; j < W; j++)
        xd[i][j] = (double)(x[i][j]);



    int Bx = 32, By = 32;

    // Application de la dct par bloc (compression)
    dct2dim_bloc(xd, tdct, H, W, Bx, By, step);
    My_dct2dim_bloc(xd, tdct1, H, W, Bx, By, step);
    fprintf(stderr, "\nMa DCT par bloc directe effectuée\n");


    // définition des matrices d'erreurs (et quantification)
    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++)
            err[i][j] = (int)tdct[i][j];


    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            err1[i][j] = quantiz(tdct1[i][j], step);
            tdct1[i][j] = (double) (err[i][j]);
        }


    // Calcul de l'entropie de la compression
    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie de la dct par bloc donnée = %g [bits/pixel]\n", entro);

    double entro1 = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", entro1);


    // Application de la dct par bloc inverse (décompression)
    dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By);
    My_dct2dim_bloc_inv(tdct1, xrecd1, H, W, Bx, By);
    fprintf(stderr, "\nMa DCT par block inverse effectuée\n");


    // écriture de la matrice image finale
    finale(xrecd, xrec, H, W);
    SaveIntImage_pgm_tronc(nom_err, err, H, W);
    ecriture_pgm(nom_out, xrec, W, H);

    finale(xrecd1, xrec1, H, W);
    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);


    // libération de mémoire
    dalocd(xrecd, H);
    dalocd(xrecd1, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    dalocd(tdct1, H);
    */
    ////////////////////////////////////////     FIN DCT PAR BLOCS    //////////////////////////////////////////////////

    /////////////////////////////////////////   CODE TSIN  /////////////////////////////////////////////////////////////
    /*
    double **tdst1 = alocamd(H, W);  // Mon image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xd[i][j] = (double)(x[i][j]);

    // Application de la compression dst
    My_dst2dim(xd, tdst1, H, W);
    fprintf(stderr, "\nMa DST directe effectuée\n");

    // Quantification des coefficients
    quantification(tdst1, err1, H, W, step);

    // Calcul de l'entropie de la compression
    double Myentro = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", Myentro);

    // Application de la dst inverse (décompression)
    My_dst2dim_inv(tdst1, xrecd1, H, W);
    fprintf(stderr, "\nMa DST inverse effectuée\n");

    // écriture de la matrice image finale
    finale(xrecd1, xrec1, H, W);
    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);

    // libération de mémoire
    dalocd(xd, H);
    dalocd(xrecd1, H);
    dalocd(tdst1, H);
    */
    ////////////////////////////////////////    FIN CODE TSIN     //////////////////////////////////////////////////////

    //////////////////////////////////////    TSIN PAR BLOCS    ////////////////////////////////////////////////////////
    /*
    double **tdst1 = alocamd(H, W); // Mon image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
      for (j = 0; j < W; j++)
        xd[i][j] = (double)(x[i][j]);



    int Bx = 8, By = 8;

    // Application de la dct par bloc (compression)
    My_dst2dim_bloc(xd, tdst1, H, W, Bx, By, step);
    fprintf(stderr, "\nMa DST par bloc directe effectuée\n");


    // définition des matrices d'erreurs (et quantification)
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            err1[i][j] = quantiz(tdst1[i][j], step);
            tdst1[i][j] = (double) (err1[i][j]);
        }


    // Calcul de l'entropie de la compression
    double entro1 = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", entro1);


    // Application de la dct par bloc inverse (décompression)
    My_dst2dim_bloc_inv(tdst1, xrecd1, H, W, Bx, By);
    fprintf(stderr, "\nMa DST par block inverse effectuée\n");


    // écriture de la matrice image finale
    finale(xrecd1, xrec1, H, W);
    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);


    // libération de mémoire
    dalocd(xrecd1, H);
    dalocd(xd, H);
    dalocd(tdst1, H);
    */
    ///////////////////////////////////////     FIN TSIN  PAR BLOCS    /////////////////////////////////////////////////

    dalocuc(x, H);
    dalocuc(xrec, H);
    daloci(err, H);
    return EXIT_FAILURE;
}

