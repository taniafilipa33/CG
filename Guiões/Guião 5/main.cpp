#include <stdio.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <ctime>

float alfa = 0.0f, beta = 0.5f, radius = 100.0f, angle = 0;
float camX, camY, camZ;
float *translacoes;

void arvoresPosicao(){

    translacoes = (float*) malloc(200*2*sizeof(float));
    float x,z;
    int l = 0;
    srand((unsigned) time(NULL));
    for(int i=0; i<200; i++){
        x = 0;
        z = 0;

        while(x*x + z*z < 1600){
            x = ((float)rand()/(float)(RAND_MAX)) * 200 - 100;
            z = ((float)rand()/(float)(RAND_MAX)) * 200 - 100;
        }

        translacoes[l++] = x;
        translacoes[l++] = z;
    }
}

void arvores(){
    int pos = 0;
    float x,z;
    float g1=0.4;
    float g2=0.6;

    for(int i=0; i<100; i++){
        glPushMatrix();
        x = translacoes[pos++];
        z = translacoes[pos++];
        glTranslatef(x,0,z);

        glColor3f(0.4, 0.2, 0);
        glRotatef(-90,1,0,0);
        glutSolidCone(1,3,4,4);

        glTranslatef(0, 0, 2);
        if (i % 2 == 0){
            glColor3f(0, g1, 0);
        }
        else {
            glColor3f(0, g2, 0);
        }

        glutSolidCone(3,5,8,8);
        glPopMatrix();
    }
}

void desenhaTeapots(int nTeapots){
    for(int teapot=0; teapot<nTeapots; teapot++) {
        float alpha = teapot * 360/nTeapots;
        glPushMatrix();
        glTranslatef(0,3,0);
        glRotatef(alpha,0,1,0);
        glTranslatef(30,0,0);
        glColor3f(0.9,0,0);
        glutSolidTeapot(2);
        glPopMatrix();
    }
}

void desenhaTeapotsCentro(int nTeapots){
    for(int teapot = 0; teapot<nTeapots; teapot++) {
        float alpha = teapot * 360/nTeapots;
        glPushMatrix();
        glTranslatef(0,3,0);
        glRotatef(alpha,0,1,0);
        glTranslatef(15,0,0);
        glColor3f(0.3,0,1);
        glutSolidTeapot(2);
        glPopMatrix();
    }
}

void spherical2Cartesian() {

    camX = radius * cos(beta) * sin(alfa);
    camY = radius * sin(beta);
    camZ = radius * cos(beta) * cos(alfa);
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

    // set the camera
    glLoadIdentity();
    gluLookAt(camX, camY, camZ,
              0.0, 0.0, 0.0,
              0.0f, 1.0f, 0.0f);

    glColor3f(0.2f, 0.8f, 0.2f);
    glBegin(GL_TRIANGLES);
    glVertex3f(100.0f, 0, -100.0f);
    glVertex3f(-100.0f, 0, -100.0f);
    glVertex3f(-100.0f, 0, 100.0f);

    glVertex3f(100.0f, 0, -100.0f);
    glVertex3f(-100.0f, 0, 100.0f);
    glVertex3f(100.0f, 0, 100.0f);
    glEnd();
    arvores();
    glColor3f(0.08627, 0.04705, 0.15686);
    glutSolidTorus(2,5,10,10);
    desenhaTeapots(10);
    desenhaTeapotsCentro(10);
    // End of frame
    glutSwapBuffers();
}


void processKeys(unsigned char c, int xx, int yy) {

// put code to process regular keys in here

}


void processSpecialKeys(int key, int xx, int yy) {

    switch (key) {

        case GLUT_KEY_RIGHT:
            alfa -= 0.1; break;

        case GLUT_KEY_LEFT:
            alfa += 0.1; break;

        case GLUT_KEY_UP:
            beta += 0.1f;
            if (beta > 1.5f)
                beta = 1.5f;
            break;

        case GLUT_KEY_DOWN:
            beta -= 0.1f;
            if (beta < -1.5f)
                beta = -1.5f;
            break;

        case GLUT_KEY_PAGE_DOWN: radius -= 1.0f;
            if (radius < 1.0f)
                radius = 1.0f;
            break;

        case GLUT_KEY_PAGE_UP: radius += 1.0f; break;
    }
    spherical2Cartesian();
    glutPostRedisplay();

}


void printInfo() {

    printf("Vendor: %s\n", glGetString(GL_VENDOR));
    printf("Renderer: %s\n", glGetString(GL_RENDERER));
    printf("Version: %s\n", glGetString(GL_VERSION));

    printf("\nUse Arrows to move the camera up/down and left/right\n");
    printf("Home and End control the distance from the camera to the origin");
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

// Callback registration for keyboard processing
    glutKeyboardFunc(processKeys);
    glutSpecialFunc(processSpecialKeys);

//  OpenGL settings
    glEnable(GL_DEPTH_TEST);
    //glEnable(GL_CULL_FACE);
    arvoresPosicao();
    spherical2Cartesian();

    printInfo();

// enter GLUT's main cycle
    glutMainLoop();

    return 1;
}