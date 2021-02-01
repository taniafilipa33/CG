#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <math.h>
#include <vector>
#include <cstdio>
#include <string>
#include <fstream>
#include "tinyxml2.h"
#include "tinyxml2.cpp"

using namespace tinyxml2;
using namespace std;

class Coordinate {
public:
    float x, y, z;
};

class Figura {
public:
    std::vector<Coordinate> figura;
};

std::vector<Figura> figuras;

//variàveis utilizadas para a formatação da câmara
float px = 0;
float py = 0;
float pz = 0;

float dx = 0;
float dy = 0;
float dz = 0;

//ângulos utilizados para a formatação da câmara
float alpha = M_PI/4;
float beta = M_PI/4;
float radium = 17.0;
GLenum OPTION = GL_FILL;

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

void getFiguras() {

    std::vector<string> figurasToLoad; //vector com os nomes das figuras presentes no ficheiro XML
    string generated_path = "../";

    XMLDocument doc;
    XMLError load = doc.LoadFile("../scene.xml"); //abre ficheiro XML
    if (load != XML_SUCCESS) {
        printf("Erro no ficheiro XML\n");
        return;
    }

    XMLNode *pRoot = doc.FirstChildElement("scene");
    if (pRoot == nullptr) return;
    XMLElement *scenefiguras = pRoot->FirstChildElement("model");
    for (; scenefiguras != nullptr; scenefiguras = scenefiguras->NextSiblingElement("model")) {
        string newfigura = generated_path + scenefiguras->Attribute("file");
        figurasToLoad.push_back(newfigura);
    }

    for (auto i: figurasToLoad) {
        ifstream file;
        Figura newFig;
        file.open(i);

        while (!file.eof()) {
            Coordinate newC;
            file >> newC.x >> newC.y >> newC.z;
            newFig.figura.push_back(newC);
        }

        figuras.push_back(newFig);
    }
}
void drawFiguras() {

    std::vector<Figura>::iterator it_fig;

    for ( it_fig = figuras.begin(); it_fig != figuras.end(); it_fig++){
        std::vector<Coordinate> crdnt = it_fig->figura; //coordenadas dos pontos da figura
        std::vector<Coordinate>::iterator it;


        glBegin(GL_TRIANGLES);
        int color;

        for (it = crdnt.begin(), color = 0; it != crdnt.end() ; it++, color++) {
             if (color< 3 ) {glColor3f(1.0, 0.2, 0.4);} //Rosa
              else {
                 if (color < 6) { glColor3f(0.0, 1.0, 1.0); } //Teal
                 else {color = 0;
                     glColor3f(1.0, 0.2, 0.4);} //rosa
             }

            glVertex3f(it->x, it->y, it->z); // vertice(x,y,z)
        }
        glEnd();
    }
}

void renderScene(void) {

    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    px = radium*cosf(beta)*sinf(alpha);
    py = radium*sinf(beta);
    pz = radium*cosf(alpha)*cosf(beta);
    gluLookAt(px, py, pz,
              dx, dy, dz,
              0.0f, 1.0f, 0.0f);

// put drawing instructions here
    glBegin(GL_LINES);
// X axis in red
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-100.0f, 0.0f, 0.0f);
    glVertex3f( 100.0f, 0.0f, 0.0f);
// Y Axis in Green
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, -100.0f, 0.0f);
    glVertex3f(0.0f, 100.0f, 0.0f);
// Z Axis in Blue
    glColor3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, -100.0f);
    glVertex3f(0.0f, 0.0f, 100.0f);
    glEnd();

    glPolygonMode(GL_FRONT_AND_BACK, OPTION);

    glBegin(GL_TRIANGLES);

    getFiguras();
    drawFiguras();



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

int main(int argc, char **argv) {

// init GLUT and the window
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGBA);
    glutInitWindowPosition(100,100);
    glutInitWindowSize(800,800);
    glutCreateWindow("CG@DI-UM");

// Required callback registry
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutSpecialFunc(processSpecialKeys);


//  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}
