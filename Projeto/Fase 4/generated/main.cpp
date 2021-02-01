#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "point.h"
#include "pontoText.h"
#include <fstream>
#include <sstream>
//#include "main.cpp"
using namespace std;

//bezier patch variaveis
int patchesNum;
unsigned int *patches;
int controlPointsNum;
float *controlPoints;




void cometas(string filename)
{
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos do plano vão ser guardados
        std::vector<Ponto> pontos;

            int l = 0;
            float r=0.05;
            int slices=20;
            int stacks=20;
            float x,y,z;
            srand((unsigned) time(NULL));
            for(int i=0; i<2; i++){
                x=y=z=0;
                while(x*x + z*z +y*y < 1600){
                    x = ((float)rand()/(float)(RAND_MAX)) * 60 - 20;
                    y = 0;
                    z = ((float)rand()/(float)(RAND_MAX)) * 60 - 20;

                    float alpha = (2 * M_PI) / slices;
                    float beta = (M_PI) / stacks;

                    for (int i = 0; i < slices; i++) {
                        for (int q = 0; q < stacks / 2; q++) {

                            //Parte superior da esfera
                            //Triângulo inferior
                            pontos.push_back(Ponto(x+(r * cos((q+ 1) * beta) * sin((i + 1) * alpha)), y+(r * sin((q + 1) * beta)),
                                                   z+(r * cos((q+ 1) * beta) * cos((i + 1) * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(q * beta) * sin(i * alpha)),y+( r * sin(q * beta)),
                                                  z+( r * cos(q * beta) * cos(i * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(q * beta) * sin((i + 1) * alpha)), y+(r * sin(q * beta)),
                                                   z+(r * cos(q * beta) * cos((i + 1) * alpha))));

                            //Triângulo superior
                            pontos.push_back(Ponto(x+(r * cos(q * beta) * sin(i * alpha)), y+(r * sin(q * beta)),
                                                   z+(r * cos(q* beta) * cos(i * alpha))));
                            pontos.push_back(Ponto(x+(r * cos((q + 1) * beta) * sin((i + 1) * alpha)),y+( r * sin((q + 1) * beta)),
                                                   z+(r * cos((q + 1) * beta) * cos((i + 1) * alpha))));
                            pontos.push_back(Ponto(x+(r * cos((q + 1) * beta) * sin(i * alpha)),y+( r * sin((q + 1) * beta)),
                                                  z+( r * cos((q + 1) * beta) * cos(i * alpha))));

                            //Parte inferior da esfera
                            //Triângulo inferior
                            pontos.push_back(Ponto(x+(r * cos(-(q + 1) * beta) * sin(i * alpha)),y+( r * sin(-(q + 1) * beta)),
                                                   z+(r * cos(-(q + 1) * beta) * cos(i * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(-(q + 1) * beta) * sin((i + 1) * alpha)),y+( r * sin(-(q + 1) * beta)),
                                                   z+(r * cos(-(q+ 1) * beta) * cos((i + 1) * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(-q * beta) * sin((i + 1) * alpha)), y+(r * sin(-q * beta)),
                                                   z+(r * cos(-q * beta) * cos((i + 1) * alpha))));

                            //Triângulo superior
                            pontos.push_back(Ponto(x+(r * cos(-q * beta) * sin((i+ 1) * alpha)),y+( r * sin(-q * beta)),
                                                  z+( r * cos(-q * beta) * cos((i + 1) * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(-q * beta) * sin(i* alpha)),y+( r * sin(-q * beta)),
                                                  z+ ( r * cos(-q * beta) * cos(i * alpha))));
                            pontos.push_back(Ponto(x+(r * cos(-(q + 1) * beta) * sin(i * alpha)),y+( r * sin(-(q + 1) * beta)),
                                                  z+ ( r * cos(-(q + 1) * beta) * cos(i * alpha))));

                        }
                    }
                }
            }

        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }
    }

    fclose(f);

    }

//Função que armazena os pontos para a geração de triângulos que permitem a construção do plano
void plane(float largura, string filename) {
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL){
        // vetor onde os pontos do plano vão ser guardados
        std::vector<Ponto> pontos;
        float lado = largura / 2;

        //primeiro triangulo, virado para cima
        pontos.push_back(Ponto(lado,0.0,lado));
        pontos.push_back(Ponto(lado,0.0,-lado));
        pontos.push_back(Ponto(-lado,0.0,lado));
        //segundo triangulo, virado para cima
        pontos.push_back(Ponto(-lado, 0.0, -lado));
        pontos.push_back(Ponto(-lado, 0.0, lado));
        pontos.push_back(Ponto(lado, 0.0, -lado));

        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }
    }

    fclose(f);
}

//Função que armazena os pontos para a geração de triângulos que permitem a construção da esfera
void sphere(float raio, int slices, int stacks, string filename) {

    FILE *f;
    f = fopen(filename.c_str(), "w");
    float x, y, z;

    if (f != NULL) {

        vector<Ponto> pontos;
        vector<Ponto> normais;
        vector<PontoText> texturas;
        float alpha = 2 * M_PI / slices;
        float beta = M_PI / stacks;
        float itStacks = 1.0f/((float) stacks);
        float itSlices = 1.0f/((float) slices);
        float cSlices = 0;
        float cStacks = 0;

        //metade superior da esfera
        for (int stack = stacks/2; stack >= 0; stack--){
            for (int slice=0; slice <= slices; slice++) {

                pontos.push_back(Ponto(raio * cosf(stack*beta) * sinf(slice*alpha), raio * sinf(stack*beta), raio * cosf(stack*beta) * cosf(slice*alpha)));
                pontos.push_back(Ponto(raio * cosf((stack - 1)*beta) * sinf(slice*alpha), raio * sinf((stack - 1)*beta), raio * cosf((stack - 1)*beta) * cosf(slice*alpha)));

                normais.push_back(Ponto(cosf(stack*beta) * sinf(slice*alpha), sinf(stack*beta), cosf(stack*beta) * cosf(slice*alpha)));
                normais.push_back(Ponto(cosf((stack - 1)*beta) * sinf(slice*alpha), sinf((stack - 1)*beta), cosf((stack - 1)*beta) * cosf(slice*alpha)));

                texturas.push_back(PontoText(cSlices, cStacks));
                texturas.push_back(PontoText(cSlices, cStacks - itStacks));

                cSlices = cSlices + itSlices;
            }
            cStacks = cStacks - itStacks;
            cSlices = 0.0;
        }

        cSlices = 0;
        cStacks = 0.5;

        //metade inferior da esfera
        for (int stack = 0; stack <= stacks/2; stack++){
            for (int slice=0; slice <= slices; slice++) {
                pontos.push_back(Ponto(raio * cosf(-stack * beta) * sinf(slice*alpha), raio * sinf(-stack * beta), raio * cosf(-stack * beta) * cosf(slice*alpha)));
                pontos.push_back(Ponto(raio * cosf(-(stack + 1)*beta) * sinf(slice*alpha), raio * sinf(-(stack + 1)*beta), raio * cosf(-(stack + 1)*beta) * cosf(slice*alpha)));

                normais.push_back(Ponto(cosf(-stack * beta) * sinf(slice*alpha), sinf(-stack * beta), cosf(-stack * beta) * cosf(slice*alpha)));
                normais.push_back(Ponto(cosf(-(stack + 1)*beta) * sinf(slice*alpha), sinf(-(stack + 1)*beta), cosf(-(stack + 1)*beta) * cosf(slice*alpha)));

                texturas.push_back(PontoText(cSlices, cStacks));
                texturas.push_back(PontoText(cSlices, cStacks - itStacks));

                cSlices = cSlices + itSlices;
            }
            cStacks = cStacks - itStacks;
            cSlices = 0.0;
        }
        //imprime pontos
        fprintf(f, "%lu\n", pontos.size()*3);
        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }
        for (unsigned int normal = 0; normal < normais.size(); normal++) {
            fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
        }
        for (unsigned int textura = 0; textura < texturas.size(); textura++){
            fprintf(f, "%f %f\n", texturas[textura].getX(), texturas[textura].getY());
        }
    }

}

//Função que armazena os pontos para a geração de triângulos que permitem a construção da caixa
void box(float x, float y, float z, int nDivisoes, string filename) {

    int i;
    int q;

    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos da plano vão ser guardados
        std::vector<Ponto> pontos;

        float px = x / nDivisoes;
        float py = y / nDivisoes;
        float pz = z / nDivisoes;

        for (i = 0; i < nDivisoes; i++) {
            for (q = 0; q < nDivisoes; q++) {
                //face frontal
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2, z / 2));
                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2, z / 2));
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2 + py, z / 2));

                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2, z / 2));
                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2 + py, z / 2));
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2 + py, z / 2));

                //face trás
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2, -z / 2));
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2 + py, -z / 2));
                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2, -z / 2));

                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2, -z / 2));
                pontos.push_back(Ponto(px * i - x / 2, py * q - y / 2 + py, -z / 2));
                pontos.push_back(Ponto(px * i - x / 2 + px, py * q - y / 2 + py, -z / 2));

                //face de lateral 1
                pontos.push_back(Ponto(x / 2, py * i - y / 2, pz * q - z / 2));
                pontos.push_back(Ponto(x / 2, py * i - y / 2 + py, pz * q - z / 2));
                pontos.push_back(Ponto(x / 2, py * i - y / 2, pz * q - z / 2 + pz));


                pontos.push_back(Ponto(x / 2, py * i - y / 2, pz * q - z / 2 + pz));
                pontos.push_back(Ponto(x / 2, py * i - y / 2 + py, pz * q - z / 2));
                pontos.push_back(Ponto(x / 2, py * i - y / 2 + py, pz * q - z / 2 + pz));

                //face de lateral 2
                pontos.push_back(Ponto(-x / 2, py * i - y / 2, pz * q - z / 2));
                pontos.push_back(Ponto(-x / 2, py * i - y / 2, pz * q - z / 2 + pz));
                pontos.push_back(Ponto(-x / 2, py * i - y / 2 + py, pz * q - z / 2));


                pontos.push_back(Ponto(-x / 2, py * i - y / 2, pz * q - z / 2 + pz));
                pontos.push_back(Ponto(-x / 2, py * i - y / 2 + py, pz * q - z / 2 + pz));
                pontos.push_back(Ponto(-x / 2, py * i - y / 2 + py, pz * q - z / 2));

                //face de cima
                pontos.push_back(Ponto(px * q - x / 2, y / 2, pz * i - z / 2));
                pontos.push_back(Ponto(px * q - x / 2, y / 2, pz * i - z / 2 + pz));
                pontos.push_back(Ponto(px * q - x / 2 + px, y / 2, pz * i - z / 2));

                pontos.push_back(Ponto(px * q - x / 2, y / 2, pz * i - z / 2 + pz));
                pontos.push_back(Ponto(px * q - x / 2 + px, y / 2, pz * i - z / 2 + pz));
                pontos.push_back(Ponto(px * q - x / 2 + px, y / 2, pz * i - z / 2));

                //face de baixo
                pontos.push_back(Ponto(px * q - x / 2, -y / 2, pz * i - z / 2));
                pontos.push_back(Ponto(px * q - x / 2 + px, -y / 2, pz * i - z / 2));
                pontos.push_back(Ponto(px * q - x / 2, -y / 2, pz * i - z / 2 + pz));

                pontos.push_back(Ponto(px * q - x / 2, -y / 2, pz * i - z / 2 + pz));
                pontos.push_back(Ponto(px * q - x / 2 + px, -y / 2, pz * i - z / 2));
                pontos.push_back(Ponto(px * q - x / 2 + px, -y / 2, pz * i - z / 2 + pz));

            }
        }

        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }

    }
    fclose(f);

}

//Função que armazena os pontos para a geração de triângulos que permitem a construção do cone
void cone (float raio,float altura, int slices, int stacks,string filename){
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos do cone vão ser guardados
        std::vector<Ponto> pontos;
    float alpha = (2* M_PI)/slices;
    float beta = M_PI/stacks;


    //base
    for(int i=0; i<slices; i++) {

        float alpha = 2 * M_PI / slices;

        // circulo base


       pontos.push_back(Ponto(0,0,0));
       pontos.push_back(Ponto(raio*cos((i+1)*alpha),0,raio*-sin((i+1)*alpha)));
       pontos.push_back(Ponto(raio*cos(i*alpha),0,raio*-sin(i*alpha)));

        float raioSeguinte, alturaStack, alturaStack2, raioSeguinte2;
        //lados do cone

        for (int q = 0; q < stacks; q++) {
            alturaStack = (altura / stacks) * q;
            raioSeguinte = (altura - alturaStack) * raio / altura;
            alturaStack2 = (altura / stacks) * (q +1);
            raioSeguinte2 = (altura - alturaStack2) * raio / altura;


            //triangulo cima
           pontos.push_back(Ponto(raioSeguinte2 * cos(alpha*i), alturaStack2, raioSeguinte2 * sin(alpha*i)));
           pontos.push_back(Ponto(raioSeguinte2 * cos(alpha*(i + 1)), alturaStack2, raioSeguinte2 * sin(alpha*(i + 1))));
           pontos.push_back(Ponto(raioSeguinte * cos(alpha*i), alturaStack, raioSeguinte * sin(alpha*i)));

            // traingulo baixo
           pontos.push_back(Ponto(raioSeguinte * cos(alpha*i), alturaStack, raioSeguinte * sin(alpha*i)));
           pontos.push_back(Ponto(raioSeguinte2 * cos(alpha *(i+1)), alturaStack2, raioSeguinte2 * sin(alpha *(i + 1))));
           pontos.push_back(Ponto(raioSeguinte * cos(alpha*(i + 1)), alturaStack, raioSeguinte * sin(alpha*(i + 1))));


        }


        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
         }

         }
        fprintf(f, "%lu\n", pontos.size()*3);
     }
        fclose(f);
}
void orbita(float raio, int slices, string filename) {
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos do plano vão ser guardados
        std::vector<Ponto> pontos;
        float x, z;
        // numero de divisões do circulo
        float alpha = 2 * M_PI / slices;

        // pontos da circunferência

            for (float div = 0; div <= 2 * M_PI + alpha; div += alpha) {
                x = raio * sin(div);
                z = raio * cos(div);
                pontos.push_back(Ponto(x, 0, z));
            }

        for (unsigned int i = 0; i < pontos.size(); i++) {
            fprintf(f, "%f %f %f\n", pontos[i].getX(), pontos[i].getY(), pontos[i].getZ());
        }

    }
}

void anel (float raio, int slices, int n,string filename) {
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos do plano vão ser guardados
        std::vector<Ponto> pontos;
        float x, z;
        // numero de divisões do circulo
        float alpha = 2 * M_PI / slices;

        // pontos da circunferência
        for(int j=0;j<n;j++) {
            for (float div = 0; div <= 2 * M_PI + alpha; div += alpha) {
                x = raio * sin(div);
                z = raio * cos(div);
                pontos.push_back(Ponto(x, 0, z));
            }
            raio+=0.5;
        }

        for (unsigned int i = 0; i < pontos.size(); i++) {
            fprintf(f, "%f %f %f\n", pontos[i].getX(), pontos[i].getY(), pontos[i].getZ());
        }

    }

}

void halo (float raioI, float raioE, float slices, float stacks,string filename){
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos do cone vão ser guardados
        std::vector<Ponto> pontos;
        std::vector<Ponto> normais;
        int i, j;
        float alpha = 0, beta = 0, proximo1, proximo2;
        proximo1 = 2 * M_PI / slices;
        proximo2 = 2 * M_PI / stacks;

        vector<PontoText> texturas;
        float itStacks = 1.0f/((float) stacks);
        float itSlices = 1.0f/((float) slices);
        float cSlices = 0;
        float cStacks = 0;

            for (i = 0; i < slices; i++) {
                for (j = 0; j < stacks; j++) {

                    pontos.push_back(
                            Ponto(cos(alpha) * (raioI + raioE * cos(beta)), sin(alpha) * (raioI + raioE * cos(beta)),
                                  raioE * sin(beta)));
                    pontos.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta)),
                                           sin(alpha + proximo1) * (raioI + raioE * cos(beta)), raioE * sin(beta)));
                    pontos.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           sin(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           raioE * sin(beta + proximo2)));

                    normais.push_back(
                            Ponto(cos(alpha) * (raioI + raioE * cos(beta)), sin(alpha) * (raioI + raioE * cos(beta)),
                                  raioE * sin(beta)));
                    normais.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta)),
                                            sin(alpha + proximo1) * (raioI + raioE * cos(beta)), raioE * sin(beta)));
                    normais.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                            sin(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                            raioE * sin(beta + proximo2)));


                    texturas.push_back(PontoText(cSlices, cStacks));
                    texturas.push_back(PontoText(cSlices, cStacks - itStacks));

                    pontos.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           sin(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           raioE * sin(beta + proximo2)));
                    pontos.push_back(Ponto(cos(alpha) * (raioI + raioE * cos(beta + proximo2)),
                                           sin(alpha) * (raioI + raioE * cos(beta + proximo2)),
                                           raioE * sin(beta + proximo2)));
                    pontos.push_back(
                            Ponto(cos(alpha) * (raioI + raioE * cos(beta)), sin(alpha) * (raioI + raioE * cos(beta)),
                                  raioE * sin(beta)));




                    normais.push_back(Ponto(cos(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           sin(alpha + proximo1) * (raioI + raioE * cos(beta + proximo2)),
                                           raioE * sin(beta + proximo2)));
                    normais.push_back(Ponto(cos(alpha) * (raioI + raioE * cos(beta + proximo2)),
                                           sin(alpha) * (raioI + raioE * cos(beta + proximo2)),
                                           raioE * sin(beta + proximo2)));
                    normais.push_back(
                            Ponto(cos(alpha) * (raioI + raioE * cos(beta)), sin(alpha) * (raioI + raioE * cos(beta)),
                                  raioE * sin(beta)));

                    texturas.push_back(PontoText(cSlices, cStacks));
                    texturas.push_back(PontoText(cSlices, cStacks - itStacks));

                    texturas.push_back(PontoText(cSlices, cStacks));
                    texturas.push_back(PontoText(cSlices, cStacks - itStacks));


                    cSlices = cSlices + itSlices;

                    beta = proximo2 * (j + 1);
                }
                alpha = proximo1 * (i + 1);
                cStacks = cStacks - itStacks;
                cSlices = 0.0;
            }
        fprintf(f, "%lu\n", pontos.size()*3+3);
        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
        fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
    }
        for (unsigned int normal = 0; normal < normais.size(); normal++) {
            fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
        }
        for (unsigned int textura = 0; textura < texturas.size(); textura++){
            fprintf(f, "%f %f\n", texturas[textura].getX(), texturas[textura].getY());
        }

}

fclose(f);

}







void cintura(double r,double slices,double stacks,double slices2,double distancia,int numFaixas, string filename){

    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        std::vector<Ponto> pontos;
        double marcas;
        double phi;
        double ondex = 0, ondez = 0;
        double proximo1 = 2 * M_PI / slices2;

        double beta = (M_PI) / stacks;
        double alpha = 2 * M_PI / slices;
        double alteracao;
        double distanciamento = r * 2;

        vector<Ponto> normais;
        vector<PontoText> texturas;
        float itStacks = 1.0f/((float) stacks);
        float itSlices = 1.0f/((float) slices);
        float cSlices = 0;
        float cStacks = 0;



        for (int faixas = 0; faixas < numFaixas; faixas++) {
            marcas = 0;
            distancia = distancia + (numFaixas * distanciamento);

            /* criar alteramento de orbitas*/
            if (faixas % 2 == 0) {
                alteracao = proximo1 + proximo1 / 4;
            } else alteracao = 0;

            while (marcas < slices2) {

                phi = (proximo1 * marcas) + alteracao;
                ondex = cos(phi);
                ondez = sin(phi);
                ondex = ondex * distancia;
                ondez = ondez * distancia;

                for (int i = 0; i < slices; i = i + 1) {
                    for (int q = 0; q < stacks / 2; q++) {

                        //Parte superior da esfera
                        //Triângulo inferior
                        pontos.push_back(Ponto(ondex + r * cos((q + 1) * beta) * sin((i + 1) * alpha), r * sin((q + 1) * beta),
                                               ondez + r * cos((q + 1) * beta) * cos((i + 1) * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(q * beta) * sin(i * alpha), r * sin(q * beta),
                                               ondez + r * cos(q * beta) * cos(i * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(q * beta) * sin((i + 1) * alpha), r * sin(q * beta),
                                               ondez + r * cos(q * beta) * cos((i + 1) * alpha)));

                        normais.push_back(Ponto(ondex + cos((q + 1) * beta) * sin((i + 1) * alpha),  sin((q + 1) * beta),
                                                ondez + cos((q + 1) * beta) * cos((i + 1) * alpha)));
                        normais.push_back(Ponto(ondex + cos(q * beta) * sin(i * alpha), sin(q * beta),
                                                ondez + cos(q * beta) * cos(i * alpha)));
                        normais.push_back(Ponto(ondex + cos(q * beta) * sin((i + 1) * alpha),  sin(q * beta),
                                                ondez + cos(q * beta) * cos((i + 1) * alpha)));

                        //Triângulo superior
                        pontos.push_back(Ponto(ondex + r * cos(q * beta) * sin(i * alpha), r * sin(q * beta),
                                               ondez + r * cos(q * beta) * cos(i * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos((q + 1) * beta) * sin((i + 1) * alpha), r * sin((q + 1) * beta),
                                               ondez + r * cos((q + 1) * beta) * cos((i + 1) * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos((q + 1) * beta) * sin(i * alpha), r * sin((q + 1) * beta),
                                               ondez + r * cos((q + 1) * beta) * cos(i * alpha)));

                        normais.push_back(Ponto(ondex +  cos(q * beta) * sin(i * alpha),  sin(q * beta),
                                                ondez +  cos(q * beta) * cos(i * alpha)));
                        normais.push_back(Ponto(ondex +  cos((q + 1) * beta) * sin((i + 1) * alpha),  sin((q + 1) * beta),
                                                ondez +  cos((q + 1) * beta) * cos((i + 1) * alpha)));
                        normais.push_back(Ponto(ondex + cos((q + 1) * beta) * sin(i * alpha),  sin((q + 1) * beta),
                                                ondez +  cos((q + 1) * beta) * cos(i * alpha)));

                        //Parte inferior da esfera
                        //Triângulo inferior
                        pontos.push_back(Ponto(ondex + r * cos(-(q + 1) * beta) * sin(i * alpha), r * sin(-(q + 1) * beta),
                                               ondez + r * cos(-(q + 1) * beta) * cos(i * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(-(q + 1) * beta) * sin((i + 1) * alpha), r * sin(-(q + 1) * beta),
                                               ondez + r * cos(-(q + 1) * beta) * cos((i + 1) * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(-q * beta) * sin((i + 1) * alpha), r * sin(-q * beta),
                                               ondez + r * cos(-q * beta) * cos((i + 1) * alpha)));

                        normais.push_back(Ponto(ondex +  cos(-(q + 1) * beta) * sin(i * alpha),  sin(-(q + 1) * beta),
                                                ondez +  cos(-(q + 1) * beta) * cos(i * alpha)));
                        normais.push_back(Ponto(ondex +  cos(-(q + 1) * beta) * sin((i + 1) * alpha),  sin(-(q + 1) * beta),
                                                ondez +  cos(-(q + 1) * beta) * cos((i + 1) * alpha)));
                        normais.push_back(Ponto(ondex +  cos(-q * beta) * sin((i + 1) * alpha),  sin(-q * beta),
                                                ondez +  cos(-q * beta) * cos((i + 1) * alpha)));

                        //Triângulo superior
                        pontos.push_back(Ponto(ondex + r * cos(-q * beta) * sin((i + 1) * alpha), r * sin(-q * beta),
                                               ondez + r * cos(-q * beta) * cos((i + 1) * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(-q * beta) * sin(i * alpha), r * sin(-q * beta),
                                               ondez + r * cos(-q * beta) * cos(i * alpha)));
                        pontos.push_back(Ponto(ondex + r * cos(-(q + 1) * beta) * sin(i * alpha), r * sin(-(q + 1) * beta),
                                               ondez + r * cos(-(q + 1) * beta) * cos(i * alpha)));

                        normais.push_back(Ponto(ondex +  cos(-q * beta) * sin((i + 1) * alpha),  sin(-q * beta),
                                               ondez +  cos(-q * beta) * cos((i + 1) * alpha)));
                        normais.push_back(Ponto(ondex +  cos(-q * beta) * sin(i * alpha),  sin(-q * beta),
                                               ondez +  cos(-q * beta) * cos(i * alpha)));
                        normais.push_back(Ponto(ondex +  cos(-(q + 1) * beta) * sin(i * alpha),  sin(-(q + 1) * beta),
                                               ondez +  cos(-(q + 1) * beta) * cos(i * alpha)));



                        texturas.push_back(PontoText(cSlices, cStacks));
                        texturas.push_back(PontoText(cSlices, cStacks - itStacks));

                        cSlices = cSlices + itSlices;
                    }
                    cStacks = cStacks - itStacks;
                    cSlices = 0.0;
                }

                marcas = marcas + 0.5;
            }
        }
        fprintf(f, "%lu\n", pontos.size()*3);

        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), 1.0f, pontos[ponto].getZ());
        }

        for (unsigned int normal = 0; normal < normais.size(); normal++) {
            fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
        }
        for (unsigned int textura = 0; textura < texturas.size(); textura++){
            fprintf(f, "%f %f\n", texturas[textura].getX(), texturas[textura].getY());
        }
    }
    fclose(f);
}

void disco( float inner, float outer,int pts, string filename )
{
    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {

        float itStacks = 0;
        float itSlices = 1.0f/((float) pts);
        float cSlices = 0;
        float cStacks = 0;
        std::vector<Ponto> pontos;
        vector<Ponto> normais;

        for ( int i = 0; i <= pts; ++i) {
            float angle = (i / (float) pts) * 3.14159f * 2.0f;
            pontos.push_back(Ponto(inner * cos(angle), inner * sin(angle),0));
            pontos.push_back(Ponto(outer * cos(angle), outer * sin(angle),0));

            normais.push_back(Ponto(inner * cos(angle), inner * sin(angle),0));
            normais.push_back(Ponto(outer * cos(angle), outer * sin(angle),0));


        }


        for ( int i = 0; i <= pts; ++i) {
            float angle = (i / (float) pts) * 3.14159f * 2.0f;
            pontos.push_back(Ponto(outer * cos(angle), outer * sin(angle),0));
            pontos.push_back(Ponto(inner * cos(angle), inner * sin(angle),0));

            normais.push_back(Ponto(outer * cos(angle), outer * sin(angle),0));
            normais.push_back(Ponto(inner * cos(angle), inner * sin(angle),0));


        }

        fprintf(f,"%lu\n",pontos.size()*3);
        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }
        for (unsigned int normal = 0; normal < normais.size(); normal++) {
            fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
        }
        for (unsigned int i = 0; i < pontos.size(); i++) {
            fprintf(f, "%f %f\n", 0.1f, 0.5f);
            fprintf(f, "%f %f\n", 0.9f, 0.5f);
        }
    }
    fclose(f);
}


void fundo(float raio, int slices, int stacks, string filename) {
    FILE *f;
    f = fopen(filename.c_str(), "w");
    float x, y, z;

    if (f != NULL) {

        vector<Ponto> pontos;
        vector<Ponto> normais;
        vector<PontoText> texturas;
        float alpha = 2 * M_PI / slices;
        float beta = M_PI / stacks;
        float itStacks = 1.0f/((float) stacks);
        float itSlices = 1.0f/((float) slices);
        float cSlices = 0;
        float cStacks = 0;

        //metade superior da esfera
        for (int stack = stacks/2; stack >= 0; stack--){
            for (int slice=0; slice <= slices; slice++) {
                pontos.push_back(Ponto(raio * cosf((stack - 1)*beta) * sinf(slice*alpha), raio * sinf((stack - 1)*beta), raio * cosf((stack - 1)*beta) * cosf(slice*alpha)));
                pontos.push_back(Ponto(raio * cosf(stack*beta) * sinf(slice*alpha), raio * sinf(stack*beta), raio * cosf(stack*beta) * cosf(slice*alpha)));

                normais.push_back(Ponto(cosf((stack - 1)*beta) * sinf(slice*alpha), sinf((stack - 1)*beta), cosf((stack - 1)*beta) * cosf(slice*alpha)));
                normais.push_back(Ponto(cosf(stack*beta) * sinf(slice*alpha), sinf(stack*beta), cosf(stack*beta) * cosf(slice*alpha)));

                texturas.push_back(PontoText(cSlices, cStacks - itStacks));
                texturas.push_back(PontoText(cSlices, cStacks));


                cSlices = cSlices + itSlices;
            }
            cStacks = cStacks - itStacks;
            cSlices = 0.0;
        }

        cSlices = 0;
        cStacks = 0.5;

        //metade inferior da esfera
        for (int stack = 0; stack <= stacks/2; stack++){
            for (int slice=0; slice <= slices; slice++) {
                pontos.push_back(Ponto(raio * cosf(-(stack + 1)*beta) * sinf(slice*alpha), raio * sinf(-(stack + 1)*beta), raio * cosf(-(stack + 1)*beta) * cosf(slice*alpha)));
                pontos.push_back(Ponto(raio * cosf(-stack * beta) * sinf(slice*alpha), raio * sinf(-stack * beta), raio * cosf(-stack * beta) * cosf(slice*alpha)));

                normais.push_back(Ponto(cosf(-stack * beta) * sinf(slice*alpha), sinf(-stack * beta), cosf(-stack * beta) * cosf(slice*alpha)));
                normais.push_back(Ponto(cosf(-(stack + 1)*beta) * sinf(slice*alpha), sinf(-(stack + 1)*beta), cosf(-(stack + 1)*beta) * cosf(slice*alpha)));

                texturas.push_back(PontoText(cSlices, cStacks - itStacks));
                texturas.push_back(PontoText(cSlices, cStacks));


                cSlices = cSlices + itSlices;
            }
            cStacks = cStacks - itStacks;
            cSlices = 0.0;
        }
        //imprime pontos
        fprintf(f, "%lu\n", pontos.size()*3);
        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }
        for (unsigned int normal = 0; normal < normais.size(); normal++) {
            fprintf(f, "%f %f %f\n", normais[normal].getX(), normais[normal].getY(), normais[normal].getZ());
        }
        for (unsigned int textura = 0; textura < texturas.size(); textura++){
            fprintf(f, "%f %f\n", texturas[textura].getX(), texturas[textura].getY());
        }
    }
}

void parsing (string filename) {

    ifstream file(filename.c_str());
    string firstLine, line;

    // Colocar valor dos indices de cada patch no array.
    getline(file, firstLine);
    patchesNum = atoi(firstLine.c_str());
    //C++ always requires cast in mallocs
    patches = (unsigned int*) malloc(sizeof(unsigned int) * patchesNum * 16);
    for (int p = 0; p < patchesNum; p++) {
        getline(file, line);
        istringstream indexes(line);
        string indexP;
        for (int i = 0; i < 16 && getline(indexes, indexP, ','); i++){ //getline(char, streamsize, delimitador)
            patches[p * 16 + i] = atoi(indexP.c_str());
        }
    }

    // Colocar valor dos pontos de controlo no array.
    getline(file, firstLine);
    controlPointsNum = atoi(firstLine.c_str());
    controlPoints = (float *)malloc(sizeof(float) * 3 * controlPointsNum);
    for (int cp = 0; cp < controlPointsNum; cp++) {
        getline(file, line);
        istringstream indexes(line);
        string indexCP;
        for (int i = 0; i < 3 && getline(indexes, indexCP, ','); i++){
            controlPoints[cp * 3 + i] = (float)atof(indexCP.c_str());
        }
    }
}

Ponto pontoBezier(int p, float u, float v) {
    Ponto ponto = Ponto(0.0, 0.0, 0.0);

    // Polinomio de Bernstein
    float bernsteinU[4] = { powf(1-u, 3), 3 * u * powf(1-u, 2), 3 * powf(u, 2) * (1-u), powf(u, 3) };
    float bernsteinV[4] = { powf(1-v, 3), 3 * v * powf(1-v, 2), 3 * powf(v, 2) * (1-v), powf(v, 3) };


    for (int j=0; j<4; j++)
        for (int i=0; i<4; i++) {

            //Indice dentro de um patch p
            int indexPatch = j*4+i;
            //Respetivo indice do ponto de controlo
            int indexCP = patches[p*16 + indexPatch];
            ponto = Ponto(ponto.getX() + controlPoints[indexCP * 3 + 0] * bernsteinU[j] * bernsteinV[i],
                          ponto.getY() + controlPoints[indexCP * 3 + 1] * bernsteinU[j] * bernsteinV[i],
                          ponto.getZ() + controlPoints[indexCP * 3 + 2] * bernsteinU[j] * bernsteinV[i]);


        }

    return ponto;
}

Ponto pontoBezierNormal(int p, float u, float v) {
    Ponto ponto = Ponto(0.0, 0.0, 0.0);

    for (int j=0; j<4; j++){
        for (int i=0; i<4; i++) {

            int indexPatch = j*4+i;
            int indexCP = patches[p*16 + indexPatch];
            ponto = Ponto(ponto.getX() + controlPoints[indexCP * 3 + 0],
                          ponto.getY() + controlPoints[indexCP * 3 + 1],
                          ponto.getZ() + controlPoints[indexCP * 3 + 2]);
        }
    }
    return ponto;
}

// A funcionar para Triangles
void bezier(string patchFile, string filename, int tesselation) {
    parsing(patchFile);
    vector<Ponto> pontos;
    vector<Ponto> normais;
    vector<PontoText> texturas;
    float inc = 1.0f/tesselation;

    for (int p=0; p<patchesNum; p++) {

        for (int tv = 0; tv<tesselation; tv++) {
            float v = (float) tv/tesselation;

            for (int tu = 0; tu < tesselation; tu++) {
                float u = (float) tu/tesselation;

                // triângulo superior
                pontos.push_back(pontoBezier(p, (u+inc), (v+inc)));
                pontos.push_back(pontoBezier(p, u, v));
                pontos.push_back(pontoBezier(p, u, (v+inc)));

                // triângulo inferior
                pontos.push_back(pontoBezier(p, (u+inc), (v+inc)));
                pontos.push_back(pontoBezier(p, (u+inc), v));
                pontos.push_back(pontoBezier(p, u, v));

                normais.push_back(pontoBezierNormal(p, (u+inc), (v+inc)));
                normais.push_back(pontoBezierNormal(p, u, v));
                normais.push_back(pontoBezierNormal(p, (u+inc), v));


                normais.push_back(pontoBezierNormal(p, (u+inc), (v+inc)));
                normais.push_back(pontoBezierNormal(p, u, (v+inc)));
                normais.push_back(pontoBezierNormal(p, u, v));

                texturas.push_back(PontoText((u+inc), (v+inc)));
                texturas.push_back(PontoText(u, v));
                texturas.push_back(PontoText(u, (v+inc)));

                texturas.push_back(PontoText((u+inc), (v+inc)));
                texturas.push_back(PontoText((u+inc), v));
                texturas.push_back(PontoText(u, v));
            }
        }
    }

    FILE *f;
    f = fopen(filename.c_str(), "w");

    fprintf(f, "%lu\n", pontos.size()*3);

    // Printar os pontos no ficheiro .3d
    for (int ponto=0; ponto<pontos.size(); ponto++) {
        fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
    }
    for (int ponto=0; ponto<normais.size(); ponto++) {
        fprintf(f, "%f %f %f\n", normais[ponto].getX(), normais[ponto].getY(), normais[ponto].getZ());
    }
    for (int ponto=0; ponto<texturas.size(); ponto++) {
        fprintf(f, "%f %f\n", texturas[ponto].getX(), texturas[ponto].getY());
    }

    fclose(f);

    free(patches);
    free(controlPoints);
}

int main(int argc, char* argv[]) {

    if (argc < 2) {
        printf("Parametros Inválidos.\nExiting...");
        return -1;
    }


    string tipo = argv[1];
    string generated_path = "../engine/figures/";


    if (tipo == "plane") {
        float largura = std::stof(argv[2]);
        string filename = argv[3];
        filename = generated_path + filename;

        plane(largura, filename);

        return 0;
    }
    if (tipo=="box"){

        float largura=std::stof(argv[2]);
        float comprimento=std::stof(argv[3]);
        float altura=std::stof(argv[4]);
        int nDiv=std::stof(argv[5]);

        string filename = argv[6];
        filename = generated_path + filename;
        box(largura,comprimento,altura, nDiv,filename);

        return 0;
    }
    else if (tipo=="sphere"){

        float radius=std::stof(argv[2]);
        int slices=std::stof(argv[3]);
        int stacks=std::stof(argv[4]);

        string filename = argv[5];
        filename = generated_path + filename;
        sphere (radius, slices , stacks,filename);

        return 0;
    }
    if (tipo == "fundo") {
        float raio = stof(argv[2]);
        int slices = stoi(argv[3]);
        int stacks = stoi(argv[4]);
        string filename = argv[5];
        filename = generated_path + filename;

        //slices e stacks nao podem ser zero porque vamos as usar para dividir
        if (slices > 1 && stacks > 0) {
            fundo(raio, slices, stacks, filename);
            return 0;
        }
        else {
            printf("Parametros Inválidos.\nExiting...\n");
            return -1;
        }
        return 0;
    }

    else if (tipo=="cone") {


        float radius = std::stof(argv[2]);
        float heigth = std::stof(argv[3]);
        int slices = std::stof(argv[4]);
        int stacks = std::stof(argv[5]);

        string filename = argv[6];
        filename = generated_path + filename;
        cone(radius, heigth, slices, stacks, filename);
        return 0;
    }
    else if(tipo =="disco"){
        float inner = std::stof(argv[2]);
        float outer = std::stof(argv[3]);
        int pts = std::stof(argv[4]);

        string filename = argv[5];
        filename = generated_path + filename;
        disco(inner, outer,pts,filename);
        return 0;
    }
    else if (tipo=="halo") {


        float raioInterior = std::stof(argv[2]);
        float raioExterior = std::stof(argv[3]);
        int slices = std::stof(argv[4]);
        int stacks = std::stof(argv[5]);

        string filename = argv[6];
        filename = generated_path + filename;
        halo(raioInterior, raioExterior, slices, stacks, filename);
        return 0;
    }
    else if (tipo=="cintura") {


        double r = std::stof(argv[2]);
        double slices = std::stof(argv[3]);
        double stacks = std::stof(argv[4]);
        double slices2 = std::stof(argv[5]);
        double distancia = std::stof(argv[6]);
        int numFaixas = std::stof(argv[7]);

        string filename = argv[8];
        filename = generated_path + filename;
        cintura( r,  slices,  stacks,  slices2,  distancia, numFaixas, filename);
        return 0;
    }
    else if (tipo=="orbita") {


        float raio = std::stof(argv[2]);
        int slices = std::stof(argv[3]);

        string filename = argv[4];
        filename = generated_path + filename;
        orbita(raio, slices, filename);
        return 0;
    }
    else if (tipo=="anel") {


        float raio = std::stof(argv[2]);
        int slices = std::stof(argv[3]);
        int n = std::stof(argv[4]);

        string filename = argv[5];
        filename = generated_path + filename;
        anel(raio, slices, n, filename);
        return 0;
    }
    else if (tipo=="cometas") {
        string filename = argv[2];
        filename = generated_path + filename;
        cometas( filename);
        return 0;
    }
    else if(tipo == "bezier") {
        string patchFile = argv[2];
        string filename  = argv[3];
        int tesselation  = stoi(argv[4]);
        patchFile = generated_path + patchFile;
        filename = generated_path + filename;

        bezier(patchFile, filename, tesselation);

        return 0;
    }
    // Se chegar aqui significa que os parâmetros são inválidos
    printf("Parametros Inválidos.\nExiting...");
    return -1;
}
