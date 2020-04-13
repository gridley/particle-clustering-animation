#include <SDL2/SDL.h>
#include <vector>
#include <iostream>
#include <cmath>
#include "openmc/position.h"
#include "openmc/distribution_multi.h"

#define WIDTH 1920
#define HEIGHT 1080
#define NEUT_SIZE 3

// Cross sections for the whole problem
double nu = 2.5;
double sig_s = 0.27;
double sig_c = 0.02;
double sig_f = sig_c / (nu-1.0);
double sig_t = sig_s + sig_c + sig_f;
double vel = 20000.0 * 100.0; // cm/s (epithermal-ish)

uint64_t dir_seed = 0;
openmc::Isotropic dir_distr;

inline double randf() { return (double)rand() / (double)RAND_MAX; }

struct Color
{
  unsigned char r;
  unsigned char g;
  unsigned char b;
};

// generate a random saturated color
Color random_color() {
  std::array<unsigned char, 3> rgb;
  rgb[0] = rand() % 256;
  rgb[1] = rand() % 256;
  rgb[2] = rand() % 256;

  if (rgb[0] > 200 and rgb[1] > 200 and rgb[2] > 200) {
    unsigned char indx = rand() % 3;
    rgb[indx] = 10;
  }

  Color result;
  result.r = rgb[0];
  result.g = rgb[1];
  result.b = rgb[2];

  // Make sure no nearly white colors that don't contrast well
  // get sampled

  return result;
}

struct Neutron
{
  openmc::Position r;
  openmc::Direction u;
  double distance_to_collision;
  bool alive;
  Color c; //coloring for progeny tracking
};

class WindowManager
{
  public:
    int window_width, window_height;
    SDL_Window* window;
    WindowManager();
    ~WindowManager();
};
WindowManager::WindowManager() : window_width(WIDTH), window_height(HEIGHT),
  window(nullptr)
{
  SDL_Init(SDL_INIT_VIDEO);
  window = SDL_CreateWindow("Particle clustering demo", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
      window_width, window_height, SDL_WINDOW_SHOWN);

  if (window == nullptr)
    std::cerr << "Failed to create window instance." << std::endl;
}
WindowManager::~WindowManager()
{
  SDL_DestroyWindow(window);
  SDL_Quit();
}

void initialize_neutrons(std::vector<Neutron>& neuts) {
  int nrows = std::sqrt(neuts.size()/2);
  int ncols = 2 * nrows;
  neuts.resize(nrows * ncols); // round a bit

  double dr = HEIGHT / (nrows + 1);
  double dc = WIDTH / (ncols + 1);

  // initialize neutron positions
  for (unsigned row=0; row<nrows; ++row) {
    for (unsigned col=0; col<ncols; ++col) {
      auto& r = neuts[row * ncols + col].r;
      r.x = dc/2.0 + col * dc;
      r.y = dr/2.0 + row * dr;
      r.z = 0;
      neuts[row * ncols + col].u = dir_distr.sample(&dir_seed);
      neuts[row * ncols + col].alive = true;
      neuts[row * ncols + col].c = random_color();
      neuts[row * ncols + col].distance_to_collision = 0;
    }
  }

  // just to be extra sure...
  int oldsize = neuts.size();
  neuts.resize(oldsize * 2);
  for (int i=0; i<oldsize; ++i) {
    neuts[oldsize+i].alive = false;
  }
}

void draw_neutrons(std::vector<Neutron>& neuts, SDL_Renderer* renderer) {
  for (auto& n: neuts) {
    if (not n.alive) continue;
    int posx = n.r.x;
    int posy = n.r.y;
    if (posx < NEUT_SIZE/2) n.r.x = WIDTH-NEUT_SIZE;
    if (posy < NEUT_SIZE/2) n.r.y = HEIGHT-NEUT_SIZE;
    if (posx > WIDTH-NEUT_SIZE/2) n.r.x = NEUT_SIZE/2;
    if (posy > HEIGHT-NEUT_SIZE/2) n.r.y = NEUT_SIZE/2;
    SDL_SetRenderDrawColor(renderer, n.c.r, n.c.g, n.c.b, 180);
    posx--;
    posy--;
    #pragma unroll
    for (unsigned char xx=0; xx<NEUT_SIZE; ++xx) {
      #pragma unroll
      for (unsigned char yy=0; yy<NEUT_SIZE; ++yy) {
        SDL_RenderDrawPoint(renderer, posx+xx, posy+yy);
      }
    }
  }
}

void move_neutrons(std::vector<Neutron>& neuts, double dt) {
  for (int nn=0; nn<neuts.size(); ++nn) {
    Neutron& n = neuts[nn];
    if (not n.alive) continue;
    if (n.distance_to_collision < 0) {
      double xi = randf();
      if (xi < sig_c / sig_t) {
        // captured!
        n.alive = false;
      } else if (xi < (sig_c + sig_f)/sig_t) {
        // fission, first calculate number of new ones to make
        int n_new = (int)nu;
        if (randf() < nu-n_new) n_new++;

        // Create fission from self:
        n_new--;

        // Now loop over the vector, either find a dead neutron or push_back one
        while(n_new > 0) {
          bool found = false;
          for (int nnn=0; nnn<neuts.size(); ++nnn) {
            Neutron& nprime = neuts[nnn];
            if (not nprime.alive) {
              nprime = n;
              nprime.u = dir_distr.sample(&dir_seed); // different direction though
              nprime.distance_to_collision = -std::log(randf()) / sig_t;
              found = true;
              break;
            }
          }
          if (not found) {
            neuts.emplace_back(n);
            neuts[neuts.size()-1].u = dir_distr.sample(&dir_seed);
            neuts[neuts.size()-1].distance_to_collision = -std::log(1-randf()) / sig_t;
          }
          n_new--;
        }
      }
      // scatter or fission does this:
      n.distance_to_collision = -std::log(randf()) / sig_t;
      n.u = dir_distr.sample(&dir_seed);
    }

    // If it's not a collision, just move some
    double dist = vel * dt;
    n.r += n.u * dist;
    n.distance_to_collision -= dist;
  }
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cerr << "Must have one command line argument: the number of neutrons to simulate." << std::endl;
    exit(1);
  }
  srand(1);
  int n_neuts = std::stoi(argv[1]);

  std::vector<Neutron> neutrons(n_neuts);
  initialize_neutrons(neutrons);
  WindowManager window;

  SDL_Renderer* renderer;
  SDL_Event event;

  renderer = SDL_CreateRenderer(window.window, -1, SDL_RENDERER_ACCELERATED || SDL_HINT_RENDER_VSYNC);
  if (renderer == nullptr)
  {
    std::cerr<< "unable to make renderer instance:"
         << std::string(SDL_GetError()) << std::endl;
  }

  bool running = true;
  while (running) {

    // Handle keyboard queue
    while(SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_QUIT:
          running=false;
          break;
      }
    }

    SDL_SetRenderDrawColor(renderer, 240, 240, 250, 255);
    SDL_RenderClear(renderer);
    move_neutrons(neutrons, 1e-6);
    draw_neutrons(neutrons, renderer);
    SDL_RenderPresent(renderer);
  }

  SDL_DestroyRenderer(renderer);
}
