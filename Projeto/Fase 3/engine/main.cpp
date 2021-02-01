#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "tinyxml2.cpp"
#include "tinyxml2.h"
#include <stdio.h>
#include <IL/il.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#include <algorithm>

#endif

using namespace tinyxml2;
using namespace std;


class Coordinate {
public:
    bool empty = true;
    float x, y, z;
};

class Color {
public:
    bool empty = true;
    float r, g, b;
};
class Rotate {
public:
    bool empty = true;
    float x, y, z, angle=-1, time=-1;
};
class Translate {
public:
    bool empty = true;
    float time;
    vector<Coordinate> catPontos;
    float catmulls[50][3];
};

class Figura {
public:
    vector<vector<Coordinate>> figura;
    int buffer;
    Rotate rotate;
    Translate translate;
    Coordinate scale;
    Color color;
};

class Tree {
public:
    Figura figure;
    vector<Tree> subtrees;
};
Tree tree;



//variàveis utilizadas para a formatação da câmara
float px = 0;
float py = 0;
float pz = 0;

float gt = 0.0;
float dx = 0;
float dy = 0;
float dz = 0;
float clocks = 0;

//ângulos utilizados para a formatação da câmara
float alpha = 0;
float beta = M_PI/4;
float radium = 20.0;

int axis = 0;
GLenum OPTION = GL_FILL;


Translate vectorToArray(Translate trans) {

    for (unsigned int it = 0; it < trans.catPontos.size(); it++) {
        Coordinate coord = trans.catPontos.at(it);
        trans.catmulls[it][0] = coord.x;
        trans.catmulls[it][1] = coord.y;
        trans.catmulls[it][2] = coord.z;
    }
    return trans;
}

void multMatrixVector(float *m, float *v, float *res) {
    for (int j = 0; j < 4; ++j) {
        res[j] = 0;
        for (int k = 0; k < 4; ++k)
            res[j] += v[k] * m[j * 4 + k];
    }
}

void getCatmullRomPoint(float t, float *p0, float *p1, float *p2, float *p3, float *pos, float *deriv) {

    // catmull-rom matrix
    float m[4][4] = {	{-0.5f,  1.5f, -1.5f,  0.5f},
                         { 1.0f, -2.5f,  2.0f, -0.5f},
                         {-0.5f,  0.0f,  0.5f,  0.0f},
                         { 0.0f,  1.0f,  0.0f,  0.0f} };

    // Compute A = M * P
    float a[3][4];

    float px[4] = { p0[0], p1[0], p2[0], p3[0] };
    float py[4] = { p0[1], p1[1], p2[1], p3[1] };
    float pz[4] = { p0[2], p1[2], p2[2], p3[2] };

    multMatrixVector(*m, px, a[0]);
    multMatrixVector(*m, py, a[1]);
    multMatrixVector(*m, pz, a[2]);

    // Compute pos = T * A
    float tv[4] =  { t*t*t, t*t, t, 1 };
    float tvd[4] = { 3*t*t, 2*t, 1, 0 };

    pos[0] = tv[0]*a[0][0] + tv[1]*a[0][1] + tv[2]*a[0][2] + tv[3]*a[0][3];
    pos[1] = tv[0]*a[1][0] + tv[1]*a[1][1] + tv[2]*a[1][2] + tv[3]*a[1][3];
    pos[2] = tv[0]*a[2][0] + tv[1]*a[2][1] + tv[2]*a[2][2] + tv[3]*a[2][3];

    // compute deriv = T' * A
    deriv[0] = tvd[0]*a[0][0] + tvd[1]*a[0][1] + tvd[2]*a[0][2] + tvd[3]*a[0][3];
    deriv[1] = tvd[0]*a[1][0] + tvd[1]*a[1][1] + tvd[2]*a[1][2] + tvd[3]*a[1][3];
    deriv[2] = tvd[0]*a[2][0] + tvd[1]*a[2][1] + tvd[2]*a[2][2] + tvd[3]*a[2][3];
}

void getGlobalCatmullRomPoint(float gt, float *pos, float *deriv, Translate trans) {

    int pointCount = trans.catPontos.size();
    float t = gt * pointCount; // this is the real global t
    int index = floor(t);  // which segment
    t = t - index; // where within  the segment.

    // indices store the points
    int indices[4];
    indices[0] = (index + pointCount-1)%pointCount;
    indices[1] = (indices[0]+1)%pointCount;
    indices[2] = (indices[1]+1)%pointCount;
    indices[3] = (indices[2]+1)%pointCount;

    getCatmullRomPoint(t, trans.catmulls[indices[0]], trans.catmulls[indices[1]], trans.catmulls[indices[2]], trans.catmulls[indices[3]], pos, deriv);
}

void renderCatmullRomCurve(Translate trans) {
    // desenhar a curva usando segmentos de reta - GL_LINE_LOOP
    float pos[3], deriv[3];
    glBegin(GL_LINE_LOOP);
    for (float gt = 0; gt <= 1; gt += 0.001) {
        getGlobalCatmullRomPoint(gt, pos, deriv, trans);
        glColor3f(0.15,0.15,0.15);
        glVertex3f(pos[0], pos[1], pos[2]);
    }
    glEnd();
}

void renderCatmullTranslate(Translate trans) {
    float pos[3], deriv[3];

    gt = fmod( glutGet(GLUT_ELAPSED_TIME), (float) (trans.time * 1000) ) / (trans.time * 1000);
    getGlobalCatmullRomPoint(gt, pos, deriv, trans);
    glTranslatef(pos[0], pos[1], pos[2]);
}

void renderRotate(Rotate rot){
    if(rot.angle != -1)
        glRotatef(rot.angle, rot.x, rot.y, rot.z);
    else if(rot.time != -1){
        clocks = glutGet(GLUT_ELAPSED_TIME);
        float angle = 360 * (fmod(clocks, (float) (rot.time * 1000) ) / (rot.time * 1000));
        glRotatef(angle, rot.x, rot.y, rot.z);
    }
}


vector<Coordinate> getFigure(string figure, int* i) {
    vector<Coordinate> newFig;
    string line;
    ifstream file;
    file.open(figure);

    // Verifica se é para desenhar com VBO ou triangles e envia valor pelo i para a getGroup.
    int numSpaces, nVertices;
    if (std::getline(file, line)) {
        file.close();
        numSpaces = std::count(line.begin(), line.end(), ' ');
        if (numSpaces == 0) {
            // Criar nova stream para poder voltar atrás e registar o valor da 1ª linha e as coordenadas.
            ifstream fileVBO;
            fileVBO.open(figure);
            fileVBO >> nVertices;
            *i = nVertices;
            /* Guarda os vértices */
            while (!fileVBO.eof()) {
                Coordinate newC;
                fileVBO >> newC.x >> newC.y >> newC.z;
                newFig.push_back(newC);
            }
            fileVBO.close();
            return newFig;
        }
            // Se tiver espaços na 1ª linha é pq são coordenadas, então não é para desenhar com VBO, logo i=1
        else *i = 1;
    }

    // Simplemente regista coordenadas da figura, nos restantes casos (quando é para desenhar com triângulos).
    ifstream fileTriang;
    fileTriang.open(figure);
    while (!fileTriang.eof()) {
        Coordinate newC;
        fileTriang >> newC.x >> newC.y >> newC.z;
        newFig.push_back(newC);
    }
    fileTriang.close();
    return newFig;
}

Tree getGroup(XMLElement* node)  {

    Tree t;
    string figuresPath =  "../figures/";

    t.figure.buffer = 0;
    XMLElement* child = node->FirstChildElement();
    for (; child != nullptr; child = child->NextSiblingElement()) {
        string tag = child->Value();

        if (strcmp(tag.c_str(), "translate") == 0) {
            t.figure.translate.empty = false;
            t.figure.translate.time = child->DoubleAttribute("time");

            XMLElement* translateNode = child->FirstChildElement();
            for (; translateNode != nullptr; translateNode = translateNode->NextSiblingElement()) {
                Coordinate coord;
                coord.x = translateNode->DoubleAttribute("X");
                coord.y = translateNode->DoubleAttribute("Y");
                coord.z = translateNode->DoubleAttribute("Z");
                t.figure.translate.catPontos.push_back(coord);
            }
            t.figure.translate = vectorToArray(t.figure.translate);
        }
        else if (strcmp(tag.c_str(), "rotate") == 0) {
            t.figure.rotate.empty = false;
            t.figure.rotate.x = child->DoubleAttribute("axisX");
            t.figure.rotate.y = child->DoubleAttribute("axisY");
            t.figure.rotate.z = child->DoubleAttribute("axisZ");
            if (child->DoubleAttribute("angle")) t.figure.rotate.angle = child->DoubleAttribute("angle");
            else if (child->DoubleAttribute("time"))  t.figure.rotate.time = child->DoubleAttribute("time");
        }
        else if (strcmp(tag.c_str(), "scale") == 0) {
            t.figure.scale.empty = false;
            t.figure.scale.x = child->DoubleAttribute("X");
            t.figure.scale.y = child->DoubleAttribute("Y");
            t.figure.scale.z = child->DoubleAttribute("Z");
        }
        else if (strcmp(tag.c_str(), "color") == 0) {
            t.figure.color.empty = false;
            t.figure.color.r = child->DoubleAttribute("R");
            t.figure.color.g = child->DoubleAttribute("G");
            t.figure.color.b = child->DoubleAttribute("B");
        }
        else if (strcmp(tag.c_str(), "models") == 0) {
            XMLElement* modelsNode = child->FirstChildElement();

            for (; modelsNode != NULL; modelsNode = modelsNode->NextSiblingElement()) {
                string figureName = figuresPath + modelsNode->Attribute("file");
                int i;
                vector<Coordinate> newFig = getFigure(figureName, &i);
                t.figure.buffer = i;
                t.figure.figura.push_back(newFig);
            }
        }
        else if (strcmp(tag.c_str(), "group") == 0) {
            t.subtrees.push_back(getGroup(child));
        }
    }

    return t;
}

void getScene() {
    XMLDocument doc;
    //XMLError load = doc.LoadFile("../figures/teapotDemo.xml"); //abre ficheiro XML
    XMLError load = doc.LoadFile("../figures/scene.xml"); //abre ficheiro XML
    //XMLError load = doc.LoadFile("../figures/planetDemo.xml"); //abre ficheiro XML
    if (load != XML_SUCCESS) {
        printf("Erro no ficheiro XML\n");
        return;
    }

    XMLElement *pRoot = doc.FirstChildElement("scene");
    if (pRoot == nullptr)
        return;

    tree = getGroup(pRoot);
}

// Funçoes de desenho da informaçao carregada
void drawVBOs(vector<vector<Coordinate> > figures, int numPontos) {

    float *array = new float[numPontos];

    vector<vector<Coordinate> >::iterator fig;
    glEnableClientState(GL_VERTEX_ARRAY);
    unsigned int buffers;
    glGenBuffers(1, &buffers);
    glBindBuffer(GL_ARRAY_BUFFER, buffers);



    for (fig = figures.begin(); fig != figures.end(); fig++) {
        vector<Coordinate>::iterator it_coords;
        int it = 0;
        for (it_coords = fig->begin(); it_coords != fig->end(); it_coords++) {
            array[it++] = it_coords->x;
            array[it++] = it_coords->y;
            array[it++] = it_coords->z;
        }
        glBufferData(GL_ARRAY_BUFFER, it*sizeof(float), array, GL_STATIC_DRAW);
        glVertexPointer(3,GL_FLOAT,0,0);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, numPontos/3);
    }
    delete[] array;
    glDeleteBuffers(1, &buffers);
}

void drawFiguras(vector<vector<Coordinate> > figuras) {

    vector<vector<Coordinate> >::iterator fig;
    for ( fig = figuras.begin(); fig != figuras.end(); fig++){

        glBegin(GL_TRIANGLES);
        vector<Coordinate>::iterator coordenadas;

        for (coordenadas= fig->begin(); coordenadas != fig->end(); coordenadas++) {
            glVertex3f(coordenadas->x, coordenadas->y, coordenadas->z);
        }
        glEnd();
    }
}


void drawGroup(Tree t) {
    glPushMatrix();
    if (!t.figure.translate.empty){
        renderCatmullRomCurve(t.figure.translate);
        renderCatmullTranslate(t.figure.translate);
        t.figure.translate.catPontos.clear();
    }
    if (!t.figure.rotate.empty){
        renderRotate(t.figure.rotate);
    }
    if (!t.figure.scale.empty){
        glScalef(t.figure.scale.x, t.figure.scale.y, t.figure.scale.z);
    }
    if (!t.figure.color.empty){
        glColor3f(t.figure.color.r, t.figure.color.g, t.figure.color.b);
    }
    if (t.figure.buffer > 1){
        drawVBOs(t.figure.figura, t.figure.buffer);
    }
    else if (t.figure.buffer == 1){
        drawFiguras(t.figure.figura);
    }
    vector<Tree>::iterator it;
    for (it = t.subtrees.begin(); it != t.subtrees.end(); it++)
        drawGroup(*it);

    glPopMatrix();
}

void changeSize(int w, int h) {

    // Prevent a divide by zero, when window is too short
    // (you cant make a window with zero width).
    if(h == 0)
        h = 1;

    // compute window's aspect ratio
    float ratio = w * 1.0 / h;

    // Set the projection matrix as current
    glMatrixMode(GL_PROJECTION);
    // Load Identity Matrix
    glLoadIdentity();

    // Set the viewport to be the entire window
    glViewport(0, 0, w, h);

    // Set perspective
    gluPerspective(45.0f ,ratio, 1.0f ,1000.0f);

    // return to the model view matrix mode
    glMatrixMode(GL_MODELVIEW);
}


void renderScene(void) {

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPolygonMode(GL_FRONT_AND_BACK, OPTION);

    glLoadIdentity();
    px = radium*cosf(beta)*sinf(alpha);
    py = radium*sinf(beta);
    pz = radium*cosf(alpha)*cosf(beta);

    gluLookAt(px + dx, py + dy, pz + dz,
              dx, dy, dz,
              0.0f, 1.0f, 0.0f);


    drawGroup(tree);



    // End of frame
    glutSwapBuffers();
}

//Função que permite a interação com o teclado
void processSpecialKeys(int key_code, int xx, int yy) {
    switch (key_code) {
        case GLUT_KEY_RIGHT: alpha = alpha + 0.05;
            break;
        case GLUT_KEY_LEFT: alpha = alpha - 0.05;
            break;
        case GLUT_KEY_UP:
            if (beta + 0.05 >= 1.5)
                beta = 1.5;
            else
                beta = beta + 0.05;
            break;
        case GLUT_KEY_DOWN:
            if (beta - 0.05 <= -1.5)
                beta = -1.5;
            else
                beta = beta - 0.05;
            break;
        default:
            break;
    }
    glutPostRedisplay();
}

void processKeys(unsigned char c, int xx, int yy) {
    switch (c) {
        case 'd': dx += 0.1;
            break;
        case 'a': dx -= 0.1;
            break;
        case 'w': dy += 0.1;
            break;
        case 's': dy -= 0.1;
            break;
        case 'q': dz -= 0.1;
            break;
        case 'e': dz += 0.1;
            break;
        case 'z':
            if(axis==0) axis=1;
            else axis=0;
            break;

        case 'i':
            OPTION = GL_FILL;
            break;
        case 'o':
            OPTION = GL_LINE;
            break;
        case 'p':
            OPTION = GL_POINT;
            break;
        default:
            break;
    }
    glutPostRedisplay();
}

void Mouse_Func(int button, int state, int x, int y) {
    switch (button) {
        case GLUT_LEFT_BUTTON:
            if (radium - 1 == 0) break;
            else radium --;
            break;
        case GLUT_RIGHT_BUTTON:
            radium ++;
            break;
    }
    glutPostRedisplay();
}


int main(int argc, char **argv) {

// init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(0,50);
    glutInitWindowSize(3000,3000);
    glutCreateWindow("CG@DI-UM");

    getScene();

// Required callback registry
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutReshapeFunc(changeSize);


#ifndef __APPLE__
    glewInit();
#endif

    glutKeyboardFunc(processKeys);
    glutSpecialFunc(processSpecialKeys);
    glutMouseFunc(Mouse_Func);


//  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}
