#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include "point.h"
//#include "main.cpp"
using namespace std;


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
void sphere(float r, int slices, int stacks, string filename) {

    FILE *f;
    f = fopen(filename.c_str(), "w");

    if (f != NULL) {
        // vetor onde os pontos da esfera vão ser guardados
        std::vector<Ponto> pontos;

        float alpha = (2 * M_PI) / slices;
        float beta = (M_PI) / stacks;

        for (int i = 0; i < slices; i++) {
            for (int q = 0; q < stacks / 2; q++) {

                //Parte superior da esfera
                //Triângulo inferior
                pontos.push_back(Ponto(r * cos((q+ 1) * beta) * sin((i + 1) * alpha), r * sin((q + 1) * beta),
                           r * cos((q+ 1) * beta) * cos((i + 1) * alpha)));
                pontos.push_back(Ponto(r * cos(q * beta) * sin(i * alpha), r * sin(q * beta),
                           r * cos(q * beta) * cos(i * alpha)));
                pontos.push_back(Ponto(r * cos(q * beta) * sin((i + 1) * alpha), r * sin(q * beta),
                           r * cos(q * beta) * cos((i + 1) * alpha)));

                //Triângulo superior
                pontos.push_back(Ponto(r * cos(q * beta) * sin(i * alpha), r * sin(q * beta),
                           r * cos(q* beta) * cos(i * alpha)));
                pontos.push_back(Ponto(r * cos((q + 1) * beta) * sin((i + 1) * alpha), r * sin((q + 1) * beta),
                           r * cos((q + 1) * beta) * cos((i + 1) * alpha)));
                pontos.push_back(Ponto(r * cos((q + 1) * beta) * sin(i * alpha), r * sin((q + 1) * beta),
                           r * cos((q + 1) * beta) * cos(i * alpha)));

                //Parte inferior da esfera
                //Triângulo inferior
                pontos.push_back(Ponto(r * cos(-(q + 1) * beta) * sin(i * alpha), r * sin(-(q + 1) * beta),
                           r * cos(-(q + 1) * beta) * cos(i * alpha)));
                pontos.push_back(Ponto(r * cos(-(q + 1) * beta) * sin((i + 1) * alpha), r * sin(-(q + 1) * beta),
                           r * cos(-(q+ 1) * beta) * cos((i + 1) * alpha)));
                pontos.push_back(Ponto(r * cos(-q * beta) * sin((i + 1) * alpha), r * sin(-q * beta),
                           r * cos(-q * beta) * cos((i + 1) * alpha)));

                //Triângulo superior
                pontos.push_back(Ponto(r * cos(-q * beta) * sin((i+ 1) * alpha), r * sin(-q * beta),
                           r * cos(-q * beta) * cos((i + 1) * alpha)));
                pontos.push_back(Ponto(r * cos(-q * beta) * sin(i* alpha), r * sin(-q * beta),
                           r * cos(-q * beta) * cos(i * alpha)));
                pontos.push_back(Ponto(r * cos(-(q + 1) * beta) * sin(i * alpha), r * sin(-(q + 1) * beta),
                           r * cos(-(q + 1) * beta) * cos(i * alpha)));

            }
        }
        for (unsigned int ponto = 0; ponto < pontos.size(); ponto++) {
            fprintf(f, "%f %f %f\n", pontos[ponto].getX(), pontos[ponto].getY(), pontos[ponto].getZ());
        }

    }
    fclose(f);

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
     }
        fclose(f);
}


int main(int argc, char* argv[]) {

    if (argc < 2) {
        printf("Parametros Inválidos.\nExiting...");
        return -1;
    }


    string tipo = argv[1];
    string generated_path = "../engine/";


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
    else if (tipo=="box"){

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
    else if (tipo=="cone") {


        float radius = std::stof(argv[2]);
        float heigth = std::stof(argv[3]);
        int slices = std::stof(argv[4]);
        int stacks = std::stof(argv[5]);

        string filename = argv[6];
        filename = generated_path + filename;
        cone(radius, heigth, slices, stacks, filename);
        printf("%d%d",1%3,1%6);
        return 0;
    }

    // Se chegar aqui significa que os parâmetros são inválidos
    printf("Parametros Inválidos.\nExiting...");
    return -1;
}
