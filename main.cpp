#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <GL/glut.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <armadillo>

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
std::vector<Point> projection_laplacian;

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
  int n = std::min(std::min(vec.size(),ema26.size()),projection_laplacian.size());
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
    glVertex3f(pfactor*projection_laplacian[i  ].dat[dim_x],pfactor*projection_laplacian[i  ].dat[dim_y],0);
    glVertex3f(pfactor*projection_laplacian[i-1].dat[dim_x],pfactor*projection_laplacian[i-1].dat[dim_y],0);
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

float L2_norm(Point const & pt1,Point const & pt2)
{
  float ret = 0;
  for(int i=0;i<pt1.dat.size();i++)
  {
    ret += (pt1.dat[i] - pt2.dat[i])*(pt1.dat[i] - pt2.dat[i]);
  }
  return sqrtf(ret);
}

// nx = p.size()^2
void compute_metrics ( std::vector < Point > const & p , float * d , int nx )
{
  for(int x=0,i=0;x<nx;x++)
  for(int y=0;y<nx;y++,i++)
  {
    d[i] = L2_norm(p[x],p[y]);
  }
}

void compute_metrics_1d ( std::vector < Point > const & p , Point const & pt, float * d , int nx )
{
  for(int x=0,i=0;x<nx;x++)
  {
    d[i] = L2_norm(p[x],pt);
  }
}

// d = [nx x nx]
// W_ij = exp(-|x_i - x_j|^2/t^2)
void compute_distance_matrix ( float * W , float * d , int nx , float t )
{
  float inv_sigma_sqr = -1.0f / (t*t);
  for(int j=0,k=0;j<nx;j++,k++)
  {
    // L^2 norm
    W[k] = exp(d[k]*d[k]*inv_sigma_sqr);
  }
}


// D_ii = sum_j W_ij
void compute_weight_matrix ( float * D , float * W , int nx )
{
  for(int i=0,k=0;i<nx;i++)
  {
    D[i] = 0;
    for(int j=0;j<nx;j++,k++)
    {
      D[i] += W[k];
    }
  }
}

// solve generalized eigen-vector problem
// L f = lambda D f
// L = D - W, where D is diagonal
// so D^-1 L f = lambda f
// D^-1 ( D - W ) f = lambda f
// ( I - D^-1 W ) f = lambda f
void solve_eigen_problem ( double * A , float * D , float * W , float * eig_val , float * eig_vec , int nx , float epsilon )
{
  for(int i=0,k=0;i<nx;i++)
    for(int j=0;j<nx;j++,k++)
    {
      A[k] = (i==j)?1.0f:0.0f;
    }
  for(int i=0,k=0;i<nx;i++)
    for(int j=0;j<nx;j++,k++)
    {
      A[k] -= W[k] / (fabs(D[j]) + epsilon);
    }

  arma::mat Amat = arma::mat(A,nx,nx);

  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym ( eigval, eigvec, Amat );

  for(int i=0;i<nx;i++)
  {
    eig_val[i] = eigval[i];
  }

  for(int i=0,k=0;i<nx;i++)
    for(int j=0;j<nx;j++,k++)
    {
      eig_vec[k] = eigvec(i,j);
    }

  std::cout << "eigval:" << std::endl;
  for(int i=0;i<30;i++)
  {
    std::cout << eigval[i] << std::endl;
  }
}

float * Laplacian_eigen_maps ( std::vector < Point > const & p , int n , float sigma , int out_n )
{

  float * out = new float[n*out_n];

  float epsilon = 1e-5;

  if(n>=p.size())
  {
    n = p.size();
  }

  int nx = n*n;

  float * d = new float[nx];

  // nx = p.size()^2
  compute_metrics ( p , d , n );
 
  float * W = new float[nx];

  // d = [nx x nx]
  // W_ij = exp(-|x_i - x_j|^2/t^2)
  compute_distance_matrix ( W , d , nx , sigma );
 
  float * D = new float[n];

  // D_ii = sum_j W_ij
  compute_weight_matrix ( D , W , n );
 
  double * A = new double[nx];

  float * eig_val = new float[n];
  float * eig_vec = new float[nx];

  // solve generalized eigen-vector problem
  // L f = lambda D f
  // L = D - W, where D is diagonal
  // so D^-1 L f = lambda f
  // D^-1 ( D - W ) f = lambda f
  // ( I - D^-1 W ) f = lambda f
  solve_eigen_problem ( A , D , W , eig_val , eig_vec , n , epsilon );

  for(int x=0,k=0;x<n;x++)
    for(int y=0;y<out_n;y++,k++)
    {
      out[k] = eig_vec[k];
    }

  delete [] eig_val;
  delete [] eig_vec;
  delete [] A;
  delete [] D;
  delete [] W;
  delete [] d;

  return out;

}

void calculate_eigen_map_projection(std::vector<Point>const& p,int n,Point const& pt,float const* eig_vec,Point & proj_pt,int n_out,float sigma)
{
  if(n>=p.size())
  {
    n = p.size();
  }

  int nx = n*1;

  float * d = new float[nx];

  compute_metrics_1d ( p , pt , d , n );
 
  float * W = new float[nx];

  // d = [nx x nx]
  // W_ij = exp(-|x_i - x_j|^2/t^2)
  compute_distance_matrix ( W , d , nx , sigma );

  proj_pt.dat.clear();
  for(int v=0,k=0;v<n_out;v++)
  {
    float proj = 0;
    for(int i=0;i<n;i++,k++)
    {
      proj += d[i] * eig_vec[k];
    }
    proj_pt.dat.push_back(proj);
  }

}

void normalize(std::vector<Point>& p)
{
  std::vector<float> norm_vec(p[0].dat.size());
  for(int i=0;i<p.size();i++)
  {
    for(int k=0;k<norm_vec.size();k++)
    {
      if(fabs(p[i].dat[k])>norm_vec[k]&&fabs(p[i].dat[k])<1e3)
      {
        norm_vec[k] = fabs(p[i].dat[k]);
      }
    }
  }
  for(int k=0;k<norm_vec.size();k++)
  {
    std::cout << k << "\t" << norm_vec[k] << std::endl;
  }
  for(int i=0;i<p.size();i++)
  {
    for(int k=0;k<norm_vec.size();k++)
    {
      p[i].dat[k] /= norm_vec[k] + 0.001f;
    }
  }
}

int main(int argc, char **argv)
{
  std::cout << "Welcome to Stock Strange Attractor" << std::endl;
  read_data("HistoricalPrices.csv");
  EMA(vec,ema26,26,26);
  EMA(vec,ema12,12,26);
  MACD(ema12,ema26,macd);
  compute_projection(1,DIM,macd,projection);
  float sigma = 0.01f;
  int num_ref_pts = 100;
  float * eig_vec = Laplacian_eigen_maps ( projection , num_ref_pts , sigma , DIM );
  for(int i=0;i<projection.size();i++)
  {
    Point pt;
    calculate_eigen_map_projection(projection,num_ref_pts,projection[i],eig_vec,pt,DIM,sigma);
    projection_laplacian.push_back(pt);
  }
  normalize(projection_laplacian);
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

