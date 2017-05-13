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

    return 1;   //EXIT_SUCCESS or EXIT_FAILURE is better
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

    return 1;   //EXIT_SUCCESS or EXIT_FAILURE is better
}


void My_dct2dim(double **xd, double **tdct1, int H, int W){

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrth = sqrt((double)H);
    double sqrtw = sqrt((double)W);
    double sqrt2 = sqrt(2.);

    // Matrices utiles
    int i, j;
    double **matDct1= alocamd(H, W);

    for ( j = 0; j < W; j++)
        matDct1[0][j] = 1 / sqrth;
    for( i = 1 ; i < H; i++){
        for(j = 0; j < W; j++){
            matDct1[i][j]= (sqrt2 / sqrth) * cos(M_PI * (2 * j + 1) * i / (2*H));
        }
    }


    double **matDct2 = alocamd(W, H);

    for ( i = 0; i < H; i++){
        matDct2[i][0] = 1 / sqrtw;
        for ( j = 1; j < W; j++){
            matDct2[i][j]= (sqrt2 / sqrtw) * cos(M_PI * (2 * i + 1) * j / (2*W));
        }
    }

    double **tdct0 = alocamd(H, W); // image transition

    for (i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            tdct0[i][j] = 0.0;
            tdct1[i][j] = 0.0;
        }

    int a, b;

    for (a = 0; a < H; a++)
        for (b = 0; b <W; b++)
            for (i = 0; i < H; i++)
                tdct0[a][b] += matDct1[a][i] * xd[i][b];


    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            for (j = 0; j < W; j++)
                tdct1[a][b] += tdct0[a][j] * matDct2[j][b];

}


void My_dct2dim_inv(double **tdct1, double **xrecd1, int H, int W){

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrth = sqrt((double)H);
    double sqrtw = sqrt((double)W);
    double sqrt2 = sqrt(2.);

    // Matrices utiles
    int i, j;
    double **matDct1Inv = alocamd(W, H);        // transposée de matDct1

    for ( j = 0; j < W; j++)
        matDct1Inv[j][0] = 1 / sqrth;
    for( i = 1 ; i < H; i++){
        for(j = 0; j < W; j++){
            matDct1Inv[j][i]= (sqrt2 / sqrth) * cos(M_PI * (2 * j + 1) * i / (2*H));
        }
    }


    double **matDct2Inv = alocamd(H, W);        // transposée de matDct2

    for ( i = 0; i < H; i++){
        matDct2Inv[0][i] = 1 / sqrtw;
        for ( j = 1; j < W; j++){
            matDct2Inv[j][i]= (sqrt2 / sqrtw) * cos(M_PI * (2 * i + 1) * j / (2*W));
        }
    }


    double **xrecd0 = alocamd(H, W);

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            xrecd0[i][j] = 0.0;
            xrecd1[i][j] = 0.0;
        }


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


void quantification(double **tdct, int **err, int H, int W, int step){

    int i, j;
    for (i = 1; i < H; i++)
        for (j = 1; j < W; j++) {
            err[i][j] = quantiz(tdct[i][j], step);
            tdct[i][j] = (double)(err[i][j]);
        }
    for (j = 0; j < W; j++) {
        err[0][j] = quantiz(tdct[0][j], 1);
        tdct[0][j] = (double)(err[0][j]);
    }
    for (i = 0; i < H; i++) {
        err[i][0] = quantiz(tdct[i][0], 1);
        tdct[i][0] = (double)(err[i][0]);
    }

}


void finale(double **xrecd, unsigned char **xrec, int H, int W){

    int i, j;
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            if (xrecd[i][j] < 0.0)
                xrec[i][j] = 0;
            else if (xrecd[i][j] > 255.0)
                xrec[i][j] = 255;
            else
                xrec[i][j] = (unsigned char) (xrecd[i][j]);
        }
}


void My_dct2dim_bloc(double **xd, double **tdct1, int H, int W, int Bx, int By, int step){

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtBx = sqrt((double)Bx);
    double sqrtBy = sqrt((double)By);
    double sqrt2 = sqrt(2.);

    // Matrices utiles
    int i, j;
    double **matDct1= alocamd(Bx, By);

    for ( j = 0; j < By; j++)
        matDct1[0][j] = 1 / sqrtBx;
    for( i = 1 ; i < Bx; i++){
        for(j = 0; j < By; j++){
            matDct1[i][j]= (sqrt2 / sqrtBx) * cos(M_PI * (2 * j + 1) * i / (2*Bx));
        }
    }

    double **matDct2 = alocamd(By, Bx);

    for ( i = 0; i < Bx; i++){
        matDct2[i][0] = 1/sqrtBy;
        for ( j = 1; j < By; j++){
            matDct2[i][j]= (sqrt2 / sqrtBy) * cos(M_PI * (2 * i + 1) * j / (2*By));
        }
    }


    double **tdct0 = alocamd(H, W); // image transition test 1

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++){
            tdct0[i][j] = 0.0;
            tdct1[i][j] = 0.0;
        }

    int nbx, nby;       // indice du bloc courant
    int a, b;           // indice du pixel courant dans le bloc

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


void My_dct2dim_bloc_inv(double **tdct1, double **xrecd1, int H, int W, int Bx, int By){

    // Valeurs utiles : on ne les calcule qu'une fois car long
    double sqrtBx = sqrt((double)Bx);
    double sqrtBy = sqrt((double)By);
    double sqrt2 = sqrt(2.);

    // Matrices utiles
    int i, j;
    double **matDct1Inv = alocamd(Bx, By);      // transposée de matDct1 par bloc

    for ( j = 0; j < By; j++)
        matDct1Inv[j][0] = 1 / sqrtBx;
    for( i = 1 ; i < Bx; i++){
        for(j = 0; j < By; j++){
            matDct1Inv[j][i]= (sqrt2 / sqrtBx) * cos(M_PI * (2 * j + 1) * i / (2*Bx));
        }
    }

    double **matDct2Inv = alocamd(By, Bx);      // transposée de matDct2 par bloc

    for ( i = 0; i < Bx; i++){
        matDct2Inv[0][i] = 1/sqrtBy;
        for ( j = 1; j < By; j++){
            matDct2Inv[j][i]= (sqrt2 / sqrtBy) * cos(M_PI * (2 * i + 1) * j / (2*By));
        }
    }


    double **xrecd0 = alocamd(H, W); // image transisitoire

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            xrecd0[i][j] = 0.0;
            xrecd1[i][j] = 0.0;
        }

    int nbx, nby;       // indice du bloc courant
    int a, b;           // indice du pixel courant dans le bloc

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



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////        /////////////////////////////////////////////////////////////
////////////////////////////////////////     MAIN     //////////////////////////////////////////////////////////
///////////////////////////////////////////        /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////


int main (int argc, char *argv[]){
    char nom[200], nom_out[200], nom_err[300];
    char nom_out_My[200], nom_err_My[300];
    int W, H; // les dimensions de l'image: H = Height, W = Width
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
    unsigned char **xrec = alocamuc(H, W);      // image finale en utilisant la dct donnée
    unsigned char **xrec1 = alocamuc(H, W);     // Mon image finale

    double **xrecd = alocamd(H, W);     // image reconstruite
    double **xrecd1 = alocamd(H, W);    // Mon image reconstruite

    lecture_pgm(nom, x);    // x contient l'image initiale

    int **err = alocami(H, W);      // allocation de la matrice des erreurs de prediction
    int **err1 = alocami(H, W);     // allocation de notre matrice des erreurs de prediction

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err[i][j] = (int)(x[i][j]);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err1[i][j] = (int)(x[i][j]);

    fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));


    /////////////////////////////////////////   CODE DCT  ///////////////////////////////////////////////////////////
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

    ////////////////////////////////////////    FIN CODE DCT     /////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////////////////////////////////////    DCT PAR BLOCS    ///////////////////////////////////////////////////////

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


    ////////////////////////////////////////     FIN DCT PAR BLOCS    /////////////////////////////////////////////////


    dalocuc(x, H);
    dalocuc(xrec, H);
    daloci(err, H);

    return 1;
}

