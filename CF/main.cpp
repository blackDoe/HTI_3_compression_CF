#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fichiers.h"
#include "matrix.h"
#include "pred.h"
#include "dct.h"


int SaveIntImage_pgm(char *nom, int **im, int Height, int Width)
{
    int i, j;

    unsigned char ** ima_uc = alocamuc(Height, Width);

    int maxim = 0;

    for(i = 0; i<Height; i++)
        for(j = 0; j<Width; j++)
        {
            if(abs(im[i][j]) > maxim)
                maxim = abs(im[i][j]);
        }
    //  fprintf(stderr, "\nnom = %s  maxim = %d\n", nom, maxim);

    for(i = 0; i<Height; i++)
        for(j = 0; j<Width; j++)
        {
            ima_uc[i][j] = (unsigned char) ( (abs(im[i][j]) *255)/ maxim);
        }

    ecriture_pgm(nom, ima_uc, Width, Height);

    dalocuc(ima_uc, Height);

    return 1;
}


int SaveIntImage_pgm_tronc(char *nom, int **im, int Height, int Width)
{
    int i, j;

    unsigned char ** ima_uc = alocamuc(Height, Width);


    for(i = 0; i<Height; i++)
        for(j = 0; j<Width; j++)
        {
            if( abs(im[i][j]) >= 255.0 )
                ima_uc[i][j] = 255;
            else
                ima_uc[i][j] = (unsigned char) ( (abs(im[i][j])) );
        }

    ecriture_pgm(nom, ima_uc, Width, Height);

    dalocuc(ima_uc, Height);

    return 1;
}


int main (int argc, char *argv[])
{
    char nom[200], nom_out[200], nom_err[300];
    int W, H; // les dimensions de li'image: H = Height, W = Width
    int i, j;

    if (argc != 3)
    {
        fprintf(stderr, "Utilisation: %s <nom_image_pgm> <pas_quantification> \n", argv[0]); exit(0);
    }

    strcpy(nom, argv[1]);
    strcpy(nom_out, argv[1]); strcat(nom_out, ".out");
    strcpy(nom_err, argv[1]); strcat(nom_err, ".err");

    int step = atoi(argv[2]);

    lecture_dim(nom, &W, &H);

    fprintf(stderr, "Width = %d  Height = %d\n", W, H);

    unsigned char **x = alocamuc(H, W);
    unsigned char **xrec = alocamuc(H, W);

    lecture_pgm(nom, x); // x contient l'image initiale

    int **err = alocami(H, W); // allocation de la matric des erreurs de prediction

    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
            err[i][j] = (int)x[i][j];
    fprintf(stderr, "\nentropie initiale = %g [bits/pixel]\n", calc_entropie(err, H, W));



    // CODE DCT
    double **tdct = alocamd(H, W); // image transformee
    double **xd = alocamd(H, W);   // image initiale en double
    double **xrecd = alocamd(H, W); // image reconstruite



    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
            xd[i][j] = (double)x[i][j];


    dct2dim(xd, tdct, H, W);
    fprintf(stderr, "OK DCT directe\n");

    ////////////////////////////////////////      TEST 1      ///////////////////////////////////////////////////////



    double **matDct1= alocamd(H, W);
    for(j=0;j<W;j++)
        matDct1[0][j] = 1/(sqrt(H));
    for(i=1;i<H;i++){
        for(j=0;j<W;j++){
            matDct1[i][j]= (sqrt(2/H))*cos(M_PI*(2*j+1)*i/(2*H));
        }
    }

    double **matDct2 = alocamd(W, H);

    /*
    for(j=0;j<H;j++)
        matDct2[0][j] = 1/(sqrt(W));
    for(i=1;i<W;i++){
        for(j=0;j<H;j++){
            matDct2[i][j]= (sqrt(2/W))*cos(M_PI*(2*j+1)*i/(2*W));
        }
    }

     */

    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
            matDct2[i][j] = matDct1[j][i];


    double **tdct1 = alocamd(H, W); // image transformee test 1
    double **transition = alocamd(H, H);



    int a;
    int b;
    for (a=0;a<W;a++) {
        for (b = 0; b < H; b++) {
            for (i = 0; i < H; i++) {
                for (j = 0; j < W; j++) {
                    tdct1[a][b] += matDct1[i][a] * matDct2[b][j] * xd[i][j];
                }
            }
        }
        fprintf(stderr, "a = %d / %d \n",a,H-1);
    }



    // Exemple de produit matriciel classique C = A*B
    // avec C de taille (H,W) , A de taille (H,K) et B de taille (K,W)

    /*

    for (i=0;i<H;i++)
        for (j=0;j<W;j++)
        {
            int z;
            for (z=0;z<K;z++)
                C[i][j] += A[i][z] * B[z][j];
        }

    */

    /////////////////////////////////////    TEST 2     /////////////////////////////////////////////////////

/*
    double **tdct2 = alocamd(H, W); // image transformee test 2


    int u;
    int v;
    int y2;
    int x2;
    for(u=0;u<H;u++){
        for(v=0;v<W;v++){
            for(y2=0;y2<W;y2++){
                for(x2=0;x2<H;x2++){
                 //   if(u=0){
                 //      tdct2[u][v] += (sqrt(2/(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                 //   }
                 //   if(v=0){
                 //       tdct2[u][v] += (sqrt(2/(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                 //   }
                    tdct2[u][v] += (2/(sqrt(H*W)))*xd[x2][y2]*cos((2*x2+1)*u*M_PI/(2*H))*cos((2*y2+1)*v*M_PI/(2*W));
                }
            }
        }
        fprintf(stderr, "u = %d / %d \n",u,H-1);
    }

*/

    ////////////////////////////     Comparaison avec la "vraie" dct      ////////////////////////////////////

    fprintf(stderr, "tdct[3] = \n[");
    for(j = 0; j < W; j++) {
        fprintf(stderr, " %g ", tdct[3][j]);
    }
    fprintf(stderr, "]\n");



    fprintf(stderr, "tdct1[3] = \n[");
    for(j = 0; j < W; j++) {
        fprintf(stderr, " %g ", tdct1[3][j]);
    }
    fprintf(stderr, "]\n");

     /*

    fprintf(stderr, "tdct2[3] = \n[");
    for(j = 0; j < W; j++) {
        fprintf(stderr, " %g ", tdct2[3][j]);
    }
    fprintf(stderr, "]\n");

*/

    //////////////////////////////////////////////////////////////////////////////////////////////////////



//Quantification des coefficients de la DCT
    for(i=1;i<H;i++)
        for(j=1;j<W;j++)
        {
            err[i][j] = quantiz(tdct[i][j], step);
            tdct[i][j] = (double)err[i][j];
        }
    for(j=0;j<W;j++)
    {err[0][j] =  quantiz(tdct[0][j], 1); tdct[0][j] = (double)err[0][j];}
    for(i=0;i<H;i++)
    {err[i][0] =  quantiz(tdct[i][0], 1); tdct[i][0] = (double)err[i][0];}

    //FIN Quantification des coefficients de la DCT


    double entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

    dct2dim_inv(tdct, xrecd, H, W);
    for(i=0;i<H;i++)
        for(j=0;j<W;j++)
            if(xrecd[i][j] < 0.0)  xrec[i][j] = 0;
            else if(xrecd[i][j] > 255.0)  xrec[i][j] = 255;
            else  xrec[i][j] = (unsigned char)xrecd[i][j];
    SaveIntImage_pgm_tronc(nom_err, err, H, W);
    ecriture_pgm(nom_out, xrec, W, H);

    dalocd(xrecd,H);
    dalocd(xd,H);
    dalocd(tdct,H);


    //FIN CODE DCT


/*
  // DCT PAR BLOCS

  double **tdct = alocamd(H, W); // image transformee
  double **xd = alocamd(H, W);   // image initiale en double
  double **xrecd = alocamd(H, W); // image reconstruite

    for(i=0;i<H;i++)
      for(j=0;j<W;j++)
        xd[i][j] = (double)x[i][j];



  int Bx=4, By=4;

  dct2dim_bloc(xd, tdct, H, W, Bx, By, step);

  for(i=0;i<H;i++)
    for(j=0;j<W;j++)
      err[i][j] = (int)tdct[i][j];

  double entro = calc_entropie(err, H, W);
  fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);

  dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By);

  for(i=0;i<H;i++)
      for(j=0;j<W;j++)
        if(xrecd[i][j] < 0.0)  xrec[i][j] = 0;
        else if(xrecd[i][j] > 255.0)  xrec[i][j] = 255;
        else  xrec[i][j] = (unsigned char)xrecd[i][j];

    SaveIntImage_pgm_tronc(nom_err, err, H, W);
    ecriture_pgm(nom_out, xrec, W, H);
     dalocd(xrecd,H);
    dalocd(xd,H);
    dalocd(tdct,H);


  // FIN DCT PAR BLOCS
*/
//


//PREDICTION



    // MyCodeur(x, H, W, step, err);
    // MyDecodeur(xrec, H, W, err);


    //codeur_adapt(x, err, H, W, step);
    //decodeur_adapt(err, xrec, H, W);


    entro = calc_entropie(err, H, W);
    fprintf(stderr, "\nentropie = %g [bits/pixel]\n", entro);
    SaveIntImage_pgm(nom_err, err, H, W);


    ecriture_pgm(nom_out, xrec, W, H);

//FIN PREDICTION

    dalocuc(x,H);
    dalocuc(xrec,H);
    daloci(err,H);

    return 1;
}




