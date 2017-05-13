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
    double **matDct1Inv = alocamd(W, H);


    for ( j = 0; j < W; j++)
        matDct1Inv[j][0] = 1 / sqrth;
    for( i = 1 ; i < H; i++){
        for(j = 0; j < W; j++){
            matDct1Inv[j][i]= (sqrt2 / sqrth) * cos(M_PI * (2 * j + 1) * i / (2*H));
        }
    }


    double **matDct2Inv = alocamd(H, W);

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



////////////////////////////////////////     MAIN     //////////////////////////////////////////////////////////


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
    unsigned char **xrec = alocamuc(H, W);
    unsigned char **xrec1 = alocamuc(H, W);

    double **xrecd = alocamd(H, W); // image reconstruite
    double **xrecd1 = alocamd(H, W); // Mon image reconstruite

    lecture_pgm(nom, x); // x contient l'image initiale

    int **err = alocami(H, W); // allocation de la matrice des erreurs de prediction
    int **err1 = alocami(H, W);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err[i][j] = (int)(x[i][j]);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err1[i][j] = (int)(x[i][j]);

    fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));


    /////////////////////////////////////////   CODE DCT  ///////////////////////////////////////////////
/*


    double **tdct = alocamd(H, W);  // image transformee
    double **tdct1 = alocamd(H, W);  // Mon image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xd[i][j] = (double)(x[i][j]);

    dct2dim(xd, tdct, H, W);
    fprintf(stderr, "OK DCT directe\n");


    My_dct2dim(xd, tdct1, H, W);
    fprintf(stderr, "Ma DCT directe effectuée\n");

    ////////////////////////////////////    TEST 1 bis     ///////////////////////////////////////////////////


    //  Exemple de produit matriciel classique C = A*B
    /*  avec C de taille (H, W) , A de taille (H, K) et B de taille (K, W)
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            int z;
            for (z = 0; z < K; z++)
                C[i][j] += A[i][z] * B[z][j];
        }
    */


    ///////////////////////////////////////    Comparaison avec la "vraie" dct  ////////////////////////////////////////

    /*
    fprintf(stderr, "\ntdct[3] - tdct1[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %f ", tdct[3][j] - tdct1[3][j]);
    fprintf(stderr, "]\n");
    */
    /*
    fprintf(stderr, "\ntdct1[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", tdct1[3][j]);
    fprintf(stderr, "]\n");
    */
    /*
    fprintf(stderr, "\ntdct1bis[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %f ", tdct1bis[3][j]);
    fprintf(stderr, "]\n");
    */
    /*
    fprintf(stderr, "\ntdct2[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", tdct2[3][j]);
    fprintf(stderr, "]\n");
    */

    //////////////////////////////////////////////////////////////////////////////////////////////////////

/*

    //  Quantification des coefficients de la DCT ... avec tdct

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

    //  Quantification des coefficients de la DCT ... avec tdct1

    for (i = 1; i < H; i++)
        for(j = 1; j < W; j++) {
            err1[i][j] = quantiz(tdct1[i][j], step);
            tdct1[i][j] = (double)(err1[i][j]);
        }
    for (j = 0; j < W; j++) {
        err1[0][j] =  quantiz(tdct1[0][j], 1);
        tdct1[0][j] = (double)(err1[0][j]);
        }
    for (i = 0; i < H; i++) {
        err1[i][0] =  quantiz(tdct1[i][0], 1);
        tdct1[i][0] = (double)(err1[i][0]);
        }

    //  FIN Quantification des coefficients de la DCT

    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

    double Myentro = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", Myentro);

    dct2dim_inv(tdct, xrecd, H, W);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            if (xrecd[i][j] < 0.0)
                xrec[i][j] = 0;
            else if (xrecd[i][j] > 255.0)
                xrec[i][j] = 255;
            else
                xrec[i][j] = (unsigned char) (xrecd[i][j]);
        }

    SaveIntImage_pgm_tronc(nom_err, err, H, W);

    ecriture_pgm(nom_out, xrec, W, H);

    ////////////////////////////////////////  dct2dim_inv  TEST1bis    ////////////////////////////////////////////////


    My_dct2dim_inv(tdct1, xrecd1, H, W);
    fprintf(stderr, "My DCT inverse effectuée\n");

    /*
    fprintf(stderr, "\nxrecd[3] - xrecd1[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", abs(xrecd[3][j] - xrecd1[3][j]));
    fprintf(stderr, "]\n");
     */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*


    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            if (xrecd1[i][j] < 0.0)
                xrec1[i][j] = 0;
            else if (xrecd1[i][j] > 255.0)
                xrec1[i][j] = 255;
            else
                xrec1[i][j] = (unsigned char) (xrecd1[i][j]);
        }

    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);

    /*
    fprintf(stderr, "\nxrecd[3] - xrecd1[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %d ", (int)(xrecd[3][j] - xrecd1[3][j]));
    fprintf(stderr, "]\n");
    */

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    dalocd(xrecd, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    //
    dalocd(xrecd1, H);
    dalocd(tdct1, H);

    //  FIN CODE DCT

*/

    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////


  // DCT PAR BLOCS

    double **tdct = alocamd(H, W);  // image transformee
    double **xd = alocamd(H, W);    // image initiale en double


    for (i = 0; i < H; i++)
      for (j = 0; j < W; j++)
        xd[i][j] = (double)(x[i][j]);


    int Bx = 4, By = 4;

    dct2dim_bloc(xd, tdct, H, W, Bx, By, step);

    ///////////////////////////////    dct2dim_bloc TEST 1     ////////////////////////////////////////////////////

    double **matDct1= alocamd(Bx, By);

    for ( j = 0; j < By; j++)
        matDct1[0][j] = 1 / (sqrt(Bx));
    for( i = 1 ; i < Bx; i++){
        for(j = 0; j < By; j++){
            matDct1[i][j]= (sqrt(2.) / sqrt((double) Bx)) * cos(M_PI * (2 * j + 1) * i / (2*Bx));
        }
    }


    double **matDct2 = alocamd(By, Bx);

    for ( i = 0; i < Bx; i++){
        matDct2[i][0] = 1/(sqrt(By));
        for ( j = 1; j < By; j++){
            matDct2[i][j]= (sqrt(2.) / sqrt((double) By)) * cos(M_PI * (2 * i + 1) * j / (2*By));
        }
    }

    // matDct2 par transposition
    /*
for ( i = 0; i < H; i++)
    for ( j = 0; j < W; j++)
        matDct2[i][j] = matDct1[j][i];
*/


    double **tdct0 = alocamd(H, W); // image transition test 1
    double **tdct1 = alocamd(H, W); // image transformee test 1

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++){
            tdct0[i][j] = 0.0;
            tdct1[i][j] = 0.0;
        }

    int nbx, nby;
    int a, b;

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

    // arrondi de tdct1 en entier
    /*
    for (a = 0; a < H; a++)
        for (b = 0; b < W; b++)
            tdct1[a][b] = (int)(tdct1[a][b] + 0.5);
    */

    //TEST
    /*
    fprintf(stderr, "\ntdct1[0] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", tdct1[0][j]);
    fprintf(stderr, "]\n");

    fprintf(stderr, "\ntdct[0] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", tdct[0][j]);
    fprintf(stderr, "]\n");
    */

    fprintf(stderr, "DCT block test1 directe effectuée\n");




    ////////////////////////////////////////////////////////////////////////////////////////////////////////////


    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++)
            err[i][j] = (int)tdct[i][j];



    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            err1[i][j] = quantiz(tdct1[i][j], step);
            tdct1[i][j] = (double) (err[i][j]);
        }

    // TEST
    /*
    fprintf(stderr, "\nerr1[0] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", err1[0][j]);
    fprintf(stderr, "]\n");

    fprintf(stderr, "\nerr[0] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", err[0][j]);
    fprintf(stderr, "]\n");
     */


    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

    double entro1 = calc_entropie(err1, H, W);
    fprintf(stderr, "\nMon entropie = %g [bits/pixel]\n", entro1);

    dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By);


    /////////////////////////////////////    DCT par block inverse    ////////////////////////////////////////////

    double **matDct1Inv = alocamd(Bx, By);

    for (i = 0; i < By; i++)
        for (j = 0; j < Bx; j++)
            matDct1Inv[i][j] = matDct1[j][i];

    double **matDct2Inv = alocamd(By, Bx);

    for(i = 0; i < By; i++)
        for(j = 0; j < Bx; j++)
            matDct2Inv[i][j] = matDct2[j][i];


    double **xrecd0 = alocamd(H, W); // image transisitoire
  //  double **xrecd1 = alocamd(H, W); // image reconstruite


    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            xrecd0[i][j] = 0.0;
            xrecd1[i][j] = 0.0;
        }

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

    // TEST
    /*
    fprintf(stderr, "\nxrecd0[17] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", xrecd0[17][j]);
    fprintf(stderr, "]\n");


    fprintf(stderr, "\nxrecd1[17] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", xrecd1[17][j]);
    fprintf(stderr, "]\n");


    fprintf(stderr, "\nxrecd[17] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %d ", (int) (xrecd[17][j]));
    fprintf(stderr, "]\n");
     */


    fprintf(stderr, "DCT1 par block inverse effectuée\n");


    //////////////////////////////////////////////////////////////////////////////////////////////


    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            if (xrecd[i][j] < 0.0)
                xrec[i][j] = 0;
            else if (xrecd[i][j] > 255.0)
                xrec[i][j] = 255;
            else
                xrec[i][j] = (unsigned char) (xrecd[i][j]);
        }

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            if (xrecd1[i][j] < 0.0)
                xrec1[i][j] = 0;
            else if (xrecd1[i][j] > 255.0)
                xrec1[i][j] = 255;
            else
                xrec1[i][j] = (unsigned char) (xrecd1[i][j]);
        }

    SaveIntImage_pgm_tronc(nom_err, err, H, W);
    ecriture_pgm(nom_out, xrec, W, H);

    SaveIntImage_pgm_tronc(nom_err_My, err1, H, W);
    ecriture_pgm(nom_out_My, xrec1, W, H);

    dalocd(xrecd, H);
    dalocd(xrecd0,H);
    dalocd(xrecd1, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    dalocd(tdct0, H);
    dalocd(tdct1, H);


    //  FIN DCT PAR BLOCS




    //calcul entropie
    /*
    entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);
    */

    SaveIntImage_pgm(nom_err, err, H, W);

    ecriture_pgm(nom_out, xrec, W, H);

    //FIN PREDICTION

    dalocuc(x, H);
    dalocuc(xrec, H);
    daloci(err, H);

    return 1;
}

