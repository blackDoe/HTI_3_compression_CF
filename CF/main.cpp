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



////////////////////////////////////////     MAIN     //////////////////////////////////////////////////////////


int main (int argc, char *argv[]){
    char nom[200], nom_out[200], nom_err[300];
    char nom_out1bis[200], nom_err1bis[300];
    char nom_out1[200], nom_err1[300];
    int W, H; // les dimensions de l'image: H = Height, W = Width
    int i, j;

    if (argc != 3) {
        fprintf(stderr, "Utilisation: %s <nom_image_pgm> <pas_quantification> \n", argv[0]);
        exit(0);
    }

    strcpy(nom, argv[1]);
    strcpy(nom_out, argv[1]); strcat(nom_out, ".out");
    strcpy(nom_err, argv[1]); strcat(nom_err, ".err");
    strcpy(nom_out1bis, argv[1]); strcat(nom_out1bis, "_1bis.out");
    strcpy(nom_err1bis, argv[1]); strcat(nom_err1bis, "_1bis.err");
    strcpy(nom_out1, argv[1]); strcat(nom_out1, "_1block.out");
    strcpy(nom_err1, argv[1]); strcat(nom_err1, "_1block.err");

    int step = atoi(argv[2]);

    lecture_dim(nom, &W, &H);

    fprintf(stderr, "Width = %d  Height = %d\n", W, H);

    unsigned char **x = alocamuc(H, W);
    unsigned char **xrec = alocamuc(H, W);
    unsigned char **xrec1 = alocamuc(H, W);
    unsigned char **xrec1bis = alocamuc(H, W);

    double **xrecd = alocamd(H, W); // image reconstruite

    lecture_pgm(nom, x); // x contient l'image initiale

    int **err = alocami(H, W); // allocation de la matrice des erreurs de prediction
    int **err1 = alocami(H, W);
    int **err1bis = alocami(H, W);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err[i][j] = (int)(x[i][j]);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            err1bis[i][j] = (int)(x[i][j]);

    fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*

    // CODE DCT

    double **tdct = alocamd(H, W);  // image transformee
    double **xd = alocamd(H, W);    // image initiale en double

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            xd[i][j] = (double)(x[i][j]);

    dct2dim(xd, tdct, H, W);
    fprintf(stderr, "OK DCT directe\n");


    ////////////////////////////////////////      TEST 1      ///////////////////////////////////////////////////////

/*
    double **tdct1 = alocamd(H, W); // image transformee test 1
    int a, b;
    for (a = 0; a < W; a++) {
        for (b = 0; b < H; b++)
            for (i = 0; i < H; i++)
                for (j = 0; j < W; j++)
                    tdct1[a][b] += matDct1[i][a] * matDct2[b][j] * xd[i][j];
        fprintf(stderr, "a = %d / %d \n",a,H-1);
    }
    */


    ////////////////////////////////////    TEST 1 bis     ///////////////////////////////////////////////////

    /*

    double **matDct1= alocamd(H, W);

    for ( j = 0; j < W; j++)
        matDct1[0][j] = 1 / (sqrt(H));
    for( i = 1 ; i < H; i++){
        for(j = 0; j < W; j++){
            matDct1[i][j]= (sqrt(2.) / sqrt((double) H)) * cos(M_PI * (2 * j + 1) * i / (2*H));
        }
    }


    double **matDct2 = alocamd(W, H);

    for ( i = 0; i < H; i++){
        matDct2[i][0] = 1/(sqrt(W));
        for ( j = 1; j < W; j++){
            matDct2[i][j]= (sqrt(2.) / sqrt((double) W)) * cos(M_PI * (2 * i + 1) * j / (2*W));
        }
    }

    // matDct2 par transposition
     /*
for ( i = 0; i < H; i++)
    for ( j = 0; j < W; j++)
        matDct2[i][j] = matDct1[j][i];
*/
/*

    double **tdct1b = alocamd(H, W); // image transition test 1bis
    double **tdct1bis = alocamd(H, W);   // image transformee test 1bis

    for (i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            tdct1b[i][j] = 0.0;
            tdct1bis[i][j] = 0.0;
        }

    int aa, bb, ii;

    for (aa = 0; aa < H; aa++)
        for (bb = 0; bb <W; bb++)
            for (ii = 0; ii < H; ii++)
                tdct1b[aa][bb] += matDct1[aa][ii] * xd[ii][bb];

    int jj;

    for (aa = 0; aa < H; aa++)
        for (bb = 0; bb < W; bb++)
            for (jj = 0; jj < W; jj++)
                tdct1bis[aa][bb] += tdct1b[aa][jj] * matDct2[jj][bb];

    fprintf(stderr, "DCT1bis directe effectuée\n");


    //  Exemple de produit matriciel classique C = A*B
    /*  avec C de taille (H, W) , A de taille (H, K) et B de taille (K, W)
    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            int z;
            for (z = 0; z < K; z++)
                C[i][j] += A[i][z] * B[z][j];
        }
    */

    ////////////////////////////////////////////////////   TEST 2  /////////////////////////////////////////////////////
    // Il faut clairement le refaire. Ou le laisser tomber
    /*
    double **tdct2 = alocamd(H, W); // image transformee test 2
    int u, v;
    int x2, y2;
    for ( u = 0; u < H; u++) {
        for(v = 0; v < W; v++)
            for(y2 = 0; y2 < W; y2++)
                for(x2 = 0; x2 < H; x2++) {
                 //   if (u = 0)
                 //      tdct2[u][v] += (sqrt(2/(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                 //   if (v = 0)
                 //       tdct2[u][v] += (sqrt(2/(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                    tdct2[u][v] += (2/(sqrt(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                }
        fprintf(stderr, "u = %d / %d \n",u,H-1);
    }
*/

    ///////////////////////////////////////    Comparaison avec la "vraie" dct  ////////////////////////////////////////
    /*
    fprintf(stderr, "\ntdct[3] = \n[");
    for(j = 0; j < W; j++)
        fprintf(stderr, " %g ", tdct[3][j]);
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

    //  Quantification des coefficients de la DCT ... avec tdct1bis

    for (i = 1; i < H; i++)
        for(j = 1; j < W; j++) {
            err1bis[i][j] = quantiz(tdct1bis[i][j], step);
            tdct1bis[i][j] = (double)(err1bis[i][j]);
        }
    for (j = 0; j < W; j++) {
        err1bis[0][j] =  quantiz(tdct1bis[0][j], 1);
        tdct1bis[0][j] = (double)(err1bis[0][j]);
        }
    for (i = 0; i < H; i++) {
        err1bis[i][0] =  quantiz(tdct1bis[i][0], 1);
        tdct1bis[i][0] = (double)(err1bis[i][0]);
        }

    //  FIN Quantification des coefficients de la DCT

    double entro = calc_entropie(err, H, W);

    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

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

    ////////////////////////////////////////// dct2dim_inv  TEST1bis    ////////////////////////////////////////////////


    double **matDct1Inv = alocamd(W, H);

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++)
            matDct1Inv[i][j] = matDct1[j][i];

    double **matDct2Inv = alocamd(H, W);

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++)
            matDct2Inv[i][j] = matDct2[j][i];


    double **xrecd1b = alocamd(H, W);
    double **xrecd1bis = alocamd(H, W);

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++) {
            xrecd1b[i][j] = 0.0;
            xrecd1bis[i][j] = 0.0;
        }


    for (aa = 0; aa < H; aa++)
        for (bb = 0; bb < W; bb++)
            for (ii = 0; ii < H; ii++)
                xrecd1b[aa][bb] += matDct1Inv[aa][ii] * tdct1bis[ii][bb];

    for (aa = 0; aa < H; aa++)
        for (bb = 0; bb <W; bb++)
            for (jj = 0; jj < W; jj++)
                xrecd1bis[aa][bb] += xrecd1b[aa][jj] * matDct2Inv[jj][bb];

    fprintf(stderr, "DCT1bis inverse effectué\n");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for (i = 0; i < H; i++)
        for (j = 0; j < W; j++) {
            if (xrecd1bis[i][j] < 0.0)
                xrec1bis[i][j] = 0;
            else if (xrecd1bis[i][j] > 255.0)
                xrec1bis[i][j] = 255;
            else
                xrec1bis[i][j] = (unsigned char) (xrecd1bis[i][j]);
        }

    SaveIntImage_pgm_tronc(nom_err1bis, err1bis, H, W);
    ecriture_pgm(nom_out1bis, xrec1bis, W, H);

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    dalocd(xrecd, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    //
    dalocd(xrecd1bis, H);
    dalocd(tdct1bis, H);

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


    int Bx = 64, By = 64;

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

    for(i = 0; i < H; i++)
        for(j = 0; j < W; j++)
            err1[i][j] = (int)tdct1[i][j];



    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

    dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By);
    // dct2dim_bloc_inv(tdct1, xrecd1, H, W, Bx, By);


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
    double **xrecd1 = alocamd(H, W); // image reconstruite


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

    SaveIntImage_pgm_tronc(nom_err1, err1, H, W);
    ecriture_pgm(nom_out1, xrec1, W, H);

    dalocd(xrecd, H);
    dalocd(xrecd0,H);
    dalocd(xrecd1, H);
    dalocd(xd, H);
    dalocd(tdct, H);
    dalocd(tdct0, H);
    dalocd(tdct1, H);


    //  FIN DCT PAR BLOCS



    //  PREDICTION

    //MyCodeur(x, H, W, step, err);
    //MyDecodeur(xrec, H, W, err);

    //codeur_adapt(x, err, H, W, step);
    //decodeur_adapt(err, xrec, H, W);


    entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);
    SaveIntImage_pgm(nom_err, err, H, W);

    ecriture_pgm(nom_out, xrec, W, H);

    //FIN PREDICTION

    dalocuc(x, H);
    dalocuc(xrec, H);
    daloci(err, H);

    return 1;
}

