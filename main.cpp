#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <GL/glut.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>

int DIM = 12;

struct data
{
  int index;
  float open;
  float high;
  float low;
  float close;
  void print()
  {
    std::cout << index << "\t" << open << "\t" << high << "\t" << low << "\t" << close << std::endl;
  }
};

struct Point
{
  std::vector<float> dat;
};

std::vector<data> vec;

std::vector<float> ema26;
std::vector<float> ema12;
std::vector<float> macd;
std::vector<Point> projection;

// exponential moving average
void EMA(std::vector<data> const & vec, std::vector<float> & ema , int win , int offset )
{
  if(win>offset)
  {
    std::cout << "win>offset" << std::endl;
    exit(1);
  }
  std::vector<float> W(win);
  float factor = -3.0f/(win*win);
  float normalization = 0.0f;
  for(int i=0;i<win;i++)
  {
    W[i] = exp(factor*i*i);
    normalization += W[i];
  }
  normalization = 1.0f/normalization;
  for(int i=0;i<win;i++)
  {
    W[i] *= normalization;
  }
  ema.resize(vec.size()-offset);
  for(int i=0;i<ema.size();i++)
  {
    ema[i] = 0;
    for(int k=0;k<win;k++)
    {
      ema[i] += vec[k+i].close*W[k];
    }
  }
}

// exponential moving standard deviation
void EMSD(std::vector<data> const & vec, std::vector<float> const & ema,std::vector<float> & emsd,int win,int offset)
{
  std::vector<float> W(win);
  float factor = -3.0f/(win*win);
  float normalization = 0.0f;
  for(int i=0;i<win;i++)
  {
    W[i] = exp(factor*i*i);
    normalization += W[i];
  }
  normalization = 1.0f/normalization;
  for(int i=0;i<win;i++)
  {
    W[i] *= normalization;
  }
  emsd.resize(vec.size()-offset);
  for(int i=0;i<emsd.size();i++)
  {
    emsd[i] = 0;
    for(int k=0;k<win;k++)
    {
      emsd[i] += (vec[k+i].close-ema[k+i])*(vec[k+i].close-ema[k+i])*W[k];
    }
    emsd[i] = sqrt(emsd[i]);
    emsd[i] *= sqrt(win);
  }
}

int dim_x = 0;
int dim_y = 1;

void draw()
{
  glBegin(GL_LINES);
  int n = std::min(std::min(vec.size(),ema26.size()),projection.size());
  float factor = 2.0f/n;
  float vfactor = 2.0f/25000;
  float mfactor = 10.0f;
  float pfactor = 5.0f;
  for(int i=2;i<n;i++)
  {
    glColor3f(1,1,1);
    glVertex3f(1.0f-i*factor,-1.0f+vfactor*vec[i  ].close,0);
    glVertex3f(1.0f-i*factor,-1.0f+vfactor*vec[i-1].close,0);
    glColor3f(1,0,0);
    glVertex3f(1.0f-i*factor,-1.0f+vfactor*ema26[i  ],0);
    glVertex3f(1.0f-i*factor,-1.0f+vfactor*ema26[i-1],0);
    glColor3f(0,1,0);
    glVertex3f(1.0f-i*factor,mfactor*macd[i  ],0);
    glVertex3f(1.0f-i*factor,mfactor*macd[i-1],0);
    glColor3f(1,1,0);
    glVertex3f(pfactor*projection[i  ].dat[dim_x],pfactor*projection[i  ].dat[dim_y],0);
    glVertex3f(pfactor*projection[i-1].dat[dim_x],pfactor*projection[i-1].dat[dim_y],0);
  }
  glEnd();
}

void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  draw();
  glutSwapBuffers();
}

void init()
{
  /* Use depth buffering for hidden surface elimination. */
  glEnable(GL_DEPTH_TEST);

  /* Setup the view of the cube. */
  glMatrixMode(GL_PROJECTION);
  gluPerspective( /* field of view in degree */ 40.0,
    /* aspect ratio */ 1.0,
    /* Z near */ 1.0, /* Z far */ 10.0);
  glMatrixMode(GL_MODELVIEW);
  gluLookAt(0.0, 0.0, 1.8,  /* eye is at (0,0,5) */
    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    0.0, 1.0, 0.);      /* up is in positive Y direction */

  /* Adjust cube position to be asthetic angle. */
  glTranslatef(0.0, 0.0, -1.0);
  glRotatef(0, 1.0, 0.0, 0.0);
  glRotatef(0, 0.0, 0.0, 1.0);
}

void read_data(std::string filename)
{
  std::ifstream infile(filename.c_str());
  std::string line;
  int i=0;
  while (std::getline(infile, line))
  {
    data D;
    std::stringstream iss(line);
    std::string token;
    // date
    iss >> token;
    boost::erase_all(token,"/");
    boost::erase_all(token,",");
    D.index = atoi(token.c_str());
    // open
    iss >> token;
    boost::erase_all(token,",");
    D.open = atof(token.c_str());
    // high
    iss >> token;
    boost::erase_all(token,",");
    D.high = atof(token.c_str());
    // low
    iss >> token;
    boost::erase_all(token,",");
    D.low = atof(token.c_str());
    // close
    iss >> token;
    boost::erase_all(token,",");
    D.close = atof(token.c_str());
    vec.push_back(D);
  }
  infile.close();
}

void keyboard(unsigned char key,int x,int y)
{
  switch(key)
  {
    case 27:exit(1);break;
    case 'q':dim_x++;if(dim_x>=DIM)dim_x=DIM-1;break;
    case 'a':dim_x--;if(dim_x<0)dim_x=0;break;
    case 'w':dim_y++;if(dim_y>=DIM)dim_y=DIM-1;break;
    case 's':dim_y--;if(dim_y<0)dim_y=0;break;
    default:break;
  }
}

void MACD(std::vector<float>const& v12,std::vector<float>const& v26,std::vector<float>& macd)
{
  if(v12.size()!=v26.size())
  {
    std::cout << "|v12|!=|v26|" << std::endl;
    exit(1);
  }
  macd.resize(v12.size());
  for(int i=0;i<macd.size();i++)
  {
    macd[i] = (v12[i]-v26[i])/v26[i];
  }
}

void compute_projection(int delta,int dim,std::vector<float>const& in,std::vector<Point>& out)
{
  out.resize(in.size()-delta*dim);
  for(int i=0;i<out.size();i++)
  {
    Point pt;
    for(int k=0;k<dim;k++)
    {
      pt.dat.push_back(in[i+delta*k]);
    }
    out[i] = pt;
  }
}

void idle()
{
  glutPostRedisplay();
  usleep(100);
}

int main(int argc, char **argv)
{
  std::cout << "Welcome to Stock Strange Attractor" << std::endl;
  read_data("HistoricalPrices.csv");
  EMA(vec,ema26,26,26);
  EMA(vec,ema12,12,26);
  MACD(ema12,ema26,macd);
  compute_projection(1,DIM,macd,projection);
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutCreateWindow("red 3D lighted cube");
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
  return 0;
}

