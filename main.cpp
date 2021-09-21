#include <SDL2/SDL.h>

#include <iostream>
#include <vector>
#include <chrono>

#include <cstdio>
#include <math.h>

//#define use_SSE


#include "physics.h"
#ifdef use_SSE
	#include "physics_sse.h"
	#define Sim ParticleSimSSE 
#else
	#define Sim ParticleSim
#endif



int WIDTH = 1280;
int HEIGHT = 720;

SDL_Color light_blue {109, 198, 208, 0};
SDL_Color pink {244, 139, 201};
SDL_Color black {0, 0, 0, 0};

class Window
{
public:
	Window()
	{
		if (SDL_Init(SDL_INIT_EVERYTHING)) 
		{
			fprintf(stderr, "Unable to initialize SDL:  %s\n", SDL_GetError());
		}
		window = SDL_CreateWindow("Physics",SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,WIDTH,HEIGHT,SDL_WINDOW_RESIZABLE);
		if(window==nullptr) 
		{
			SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR,"ERROR",SDL_GetError(),nullptr); 
		}
		renderer = SDL_CreateRenderer(window,-1,SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
		if(renderer==nullptr) 
		{
			SDL_ShowSimpleMessageBox(SDL_MESSAGEBOX_ERROR,"ERROR",SDL_GetError(),window); 
		}
	}

	~Window()
	{
		SDL_DestroyRenderer(renderer);
		SDL_DestroyWindow(window);
		
		SDL_Quit();	
	}

	bool alive;
	bool paused = false;
	bool L_clicked = false;
	bool R_clicked = false;
	int holding_id = -1;
	bool step_frame = false;
	
	struct Dimensions { int width; int height; };
	Dimensions get_size() const
	{
		int w,h;
		SDL_GetWindowSize(window, &w, &h);
		return {w,h};
	}

	void SetRenderDrawColor(SDL_Color c)
	{
		auto [r,g,b,a] = c;
		SDL_SetRenderDrawColor(renderer, r, g, b, a);
	}
	SDL_Window* window;
	SDL_Renderer* renderer;
};

void handle_SDL_events(Window& w)
{
	SDL_Event e;
	while (SDL_PollEvent(&e))
	{
		switch(e.type)
		{
			case SDL_QUIT:
			{
				w.alive = false;
			} break;
			case SDL_KEYDOWN:
			{
				switch(e.key.keysym.sym)
				{
					case(SDLK_ESCAPE):
					{
						SDL_Event q;
						q.type = SDL_QUIT;
						SDL_PushEvent(&q);
					} break;

					case(SDLK_SPACE):
					{
						w.paused = !w.paused;
					} break;
					case(SDLK_RIGHT):
					{
						if(w.paused)
						{
							w.step_frame = true;
						}	
					} break;

				}
			} break;
			case SDL_MOUSEBUTTONUP:
			{
				if(e.button.button==SDL_BUTTON_LEFT)
				{
					w.L_clicked = false;
				}
				else if(e.button.button==SDL_BUTTON_RIGHT)
				{
					w.R_clicked = false;
				}
			} break;
			case SDL_MOUSEBUTTONDOWN:
			{
				if(e.button.button==SDL_BUTTON_LEFT)
				{
					w.L_clicked = true;
				}
				else if(e.button.button==SDL_BUTTON_RIGHT)
				{
					w.R_clicked = true;
				}
			} break;
		}
	}
}

float timestep = 1.0/60;//seconds

int main(int argv, char** args)
{
	Window w;
	w.alive = true;

	Vector origin = {(float) WIDTH/4, 50};
	const int n_grid = 64;
	Sim sim(n_grid*n_grid);

	float grid_size = 640/n_grid;

	for(int j = 0; j<n_grid; ++j)
	{
		for(int i = 0; i<n_grid; ++i)
		{
			Vector pos = origin + Vector{i*grid_size, j*grid_size};

			sim.add_particle(pos, j==0);

			if(i>0)
			{
				Bond b;
				b.particle_1 = n_grid*j+i-1;
				b.particle_2 = n_grid*j+i;
				b.length = grid_size;
				sim.bonds.push_back(b);
			}
			if(j>0)
			{
				Bond b;
				b.particle_1 = n_grid*(j-1)+i;
				b.particle_2 = n_grid*j+i;
				b.length = grid_size;
				sim.bonds.push_back(b);
			}
			
		}
	}
	sim.calc_forces();



	int mouse_x, mouse_y;
	auto current_time =  std::chrono::high_resolution_clock::now();
	float frame_time = 0;
	float time = 0;
	while (w.alive)
	{
		w.SetRenderDrawColor(black);
		SDL_RenderClear(w.renderer);
		handle_SDL_events(w);
		if(w.L_clicked)
		{
			SDL_GetMouseState(&mouse_x, &mouse_y);
			if(w.holding_id==-1)
			{
				// Look for new rect to hold
				SDL_Point q = {mouse_x, mouse_y};
				for(int i = 0; i<sim.num_particles(); ++i)
				{
					
				#ifdef use_SSE
					Particle p = {sim.particles->x[i], sim.particles->y[i]};
				#else
					Particle p = {sim.particles[i].position[0], sim.particles[i].position[1]};
				#endif

					SDL_Rect r = {(int)p.position[0], (int)p.position[1], 10, 10};
					if(SDL_PointInRect(&q, &r))
					{
						w.holding_id = i;
						break;
					}
				}
			}
			else
			{
				#ifdef use_SSE
					sim.particles->x[w.holding_id] = (float)mouse_x;
					sim.particles->y[w.holding_id] = (float)mouse_y;
				#else
					sim.particles[w.holding_id].position = {(float)mouse_x, (float)mouse_y};
				#endif
			}
		}
		else
		{
			w.holding_id = -1;
		}

		auto new_time = std::chrono::high_resolution_clock::now();
		frame_time += std::chrono::duration_cast<std::chrono::microseconds>(new_time - current_time).count() / 1000000.0;
		current_time =  new_time;

		if(!w.paused | w.step_frame)
		{
			w.step_frame = false;
			int step_count = 0;
			while(frame_time > timestep)
			{
				time+=timestep;
				++step_count;
				frame_time -= timestep;
				sim.step(timestep);
			}
		}
		else
		{
			frame_time = 0;
		}

		w.SetRenderDrawColor(light_blue);
		for(int i = 0; i<sim.num_particles(); ++i)
		{
			#ifdef use_SSE
				Particle p = {sim.particles->x[i], sim.particles->y[i]};
			#else
				Particle p = {sim.particles[i].position[0], sim.particles[i].position[1]};
			#endif
			SDL_Rect r = {(int)p.position[0], (int)p.position[1], 2, 2};
			SDL_RenderFillRect(w.renderer, &r);
		}

		for(const auto& b : sim.bonds)
		{
		#ifdef use_SSE
			SDL_RenderDrawLine(w.renderer, 
			 sim.particles->x[b.particle_1],
			 sim.particles->y[b.particle_1],
			 sim.particles->x[b.particle_2],
			 sim.particles->y[b.particle_2]);
		#else
			SDL_RenderDrawLine(w.renderer, 
			 sim.particles[b.particle_1].position[0],
			 sim.particles[b.particle_1].position[1],
			 sim.particles[b.particle_2].position[0],
			 sim.particles[b.particle_2].position[1]);
		#endif
		}
		
		SDL_RenderPresent(w.renderer);
	}
	return 0;
}
