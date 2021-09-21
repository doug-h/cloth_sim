
#include <cassert>

const int PACKING = SDL_SIMDGetAlignment()/sizeof(float);

struct ParticleData
{
	// Structure of Arrays model for SSE
	float *x,*y;
	float *x_old,*y_old;
	float *force_x, *force_y;
	const int max_particles;
	int num_particles = 0;

	ParticleData(const int n) : max_particles{n}
	{
		std::cout << "Alignment: " << SDL_SIMDGetAlignment() << '\n';

		x = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));
		y = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));

		x_old = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));
		y_old = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));

		force_x = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));
		force_y = (float*)SDL_SIMDAlloc(max_particles*sizeof(float));
	}

	~ParticleData()
	{
		SDL_SIMDFree(x);
		SDL_SIMDFree(y);

		SDL_SIMDFree(x_old);
		SDL_SIMDFree(y_old);

		SDL_SIMDFree(force_x);
		SDL_SIMDFree(force_y);
	}
};

struct ParticleSimSSE
{
	
	const float gravity = 100;
	const __m256 ZERO = _mm256_setzero_ps();
	const __m256 HALF = _mm256_set1_ps(0.5);
	const __m256 GRAV = _mm256_set1_ps(gravity);

	float *floatbuffer;

	ParticleData* particles;
	std::vector<FixedPoint> fixed_points;
	std::vector<Bond> bonds;

	ParticleSimSSE(int num_particles)
	{
		particles = new ParticleData(num_particles);
		floatbuffer = (float*)SDL_SIMDAlloc(PACKING*sizeof(float));
	}
	~ParticleSimSSE()
	{
		delete particles;
		SDL_SIMDFree(floatbuffer);
	}

	void add_particle(Vector position = {0,0}, bool fixed=false)
	{
		particles->x[particles->num_particles] = position[0];
		particles->y[particles->num_particles] = position[1];
		particles->x_old[particles->num_particles] = position[0];
		particles->y_old[particles->num_particles] = position[1];

		if(fixed)
		{
			FixedPoint f = {position, particles->num_particles};
			fixed_points.push_back(f);
		}

		++particles->num_particles;
	}

	void calc_forces()
	{
		for(int i = 0; i<particles->max_particles; i+=PACKING)
		{
			_mm256_store_ps(&particles->force_x[i], ZERO);
			_mm256_store_ps(&particles->force_y[i], GRAV);
		}
	}
	void integrate(float timestep)
	{
		for(int i = 0; i<particles->max_particles; i+=PACKING)
		{
			//TODO: test performance of timestep multiplied by intrinsic
			__m256 dt2 = _mm256_set1_ps(timestep*timestep);
			__m256 rx = _mm256_load_ps(&particles->x[i]);
			__m256 rx_old = _mm256_load_ps(&particles->x_old[i]);
			__m256 fx = _mm256_load_ps(&particles->force_x[i]);
			__m256 ry = _mm256_load_ps(&particles->y[i]);
			__m256 ry_old = _mm256_load_ps(&particles->y_old[i]);
			__m256 fy = _mm256_load_ps(&particles->force_y[i]);

			//r += (r-r_old) + a*(dt*dt);
			_mm256_store_ps(&particles->x[i],
				_mm256_add_ps(
					rx,
					_mm256_add_ps(
						_mm256_sub_ps(rx, rx_old),
						_mm256_mul_ps(fx,dt2))));
			_mm256_store_ps(&particles->x_old[i], rx);

			_mm256_store_ps(&particles->y[i],
				_mm256_add_ps(
					ry,
					_mm256_add_ps(
						_mm256_sub_ps(ry, ry_old),
						_mm256_mul_ps(fy,dt2))));
			_mm256_store_ps(&particles->y_old[i], ry);
		}
	}

	void apply_bonds()
	{
		for(int j = 0; j<NUM_ITERATIONS; ++j)
		{
			for(int i = 0; i<bonds.size(); i+=PACKING)
			{
				__m256 x1 = _mm256_set_ps(
					particles->x[bonds[i].particle_1],
					particles->x[bonds[i+1].particle_1],
					particles->x[bonds[i+2].particle_1],
					particles->x[bonds[i+3].particle_1],
					particles->x[bonds[i+4].particle_1],
					particles->x[bonds[i+5].particle_1],
					particles->x[bonds[i+6].particle_1],
					particles->x[bonds[i+7].particle_1]
					);
				__m256 y1 = _mm256_set_ps(
					particles->y[bonds[i].particle_1],
					particles->y[bonds[i+1].particle_1],
					particles->y[bonds[i+2].particle_1],
					particles->y[bonds[i+3].particle_1],
					particles->y[bonds[i+4].particle_1],
					particles->y[bonds[i+5].particle_1],
					particles->y[bonds[i+6].particle_1],
					particles->y[bonds[i+7].particle_1]
					);
				__m256 x2 = _mm256_set_ps(
					particles->x[bonds[i].particle_2],
					particles->x[bonds[i+1].particle_2],
					particles->x[bonds[i+2].particle_2],
					particles->x[bonds[i+3].particle_2],
					particles->x[bonds[i+4].particle_2],
					particles->x[bonds[i+5].particle_2],
					particles->x[bonds[i+6].particle_2],
					particles->x[bonds[i+7].particle_2]
					);
				__m256 y2 = _mm256_set_ps(
					particles->y[bonds[i].particle_2],
					particles->y[bonds[i+1].particle_2],
					particles->y[bonds[i+2].particle_2],
					particles->y[bonds[i+3].particle_2],
					particles->y[bonds[i+4].particle_2],
					particles->y[bonds[i+5].particle_2],
					particles->y[bonds[i+6].particle_2],
					particles->y[bonds[i+7].particle_2]
					);
				__m256 length = _mm256_set_ps(
					bonds[i].length,
					bonds[i+1].length,
					bonds[i+2].length,
					bonds[i+3].length,
					bonds[i+4].length,
					bonds[i+5].length,
					bonds[i+6].length,
					bonds[i+7].length
					);

				__m256 delta_x = _mm256_sub_ps(x2,x1);
				__m256 delta_y = _mm256_sub_ps(y2,y1);

				__m256 deltalength = _mm256_sqrt_ps(_mm256_add_ps(_mm256_mul_ps(delta_x,delta_x),_mm256_mul_ps(delta_y,delta_y)));
				__m256 delta = _mm256_mul_ps(HALF, _mm256_div_ps(_mm256_sub_ps(deltalength, length),deltalength));
				delta_x = _mm256_mul_ps(delta_x,delta);
				delta_y = _mm256_mul_ps(delta_y,delta);

				_mm256_store_ps(floatbuffer, delta_x);
				float* f = floatbuffer;
				for(int k = PACKING-1; k>=0; --k)
				{
					particles->x[bonds[i+k].particle_1] += *f;
					particles->x[bonds[i+k].particle_2] -= *f;
					f++;
				}

				_mm256_store_ps(floatbuffer, delta_y);
				f = floatbuffer;
				for(int k = PACKING-1; k>=0; --k)
				{
					particles->y[bonds[i+k].particle_1] += *f;
					particles->y[bonds[i+k].particle_2] -= *f;
					f++;
				}
			}

			for(const auto& p : fixed_points)
			{
				particles->x[p.particle] = p.point[0];
				particles->y[p.particle] = p.point[1];
			}
		}
	}

	void step(float timestep)
	{
		integrate(timestep);
		apply_bonds();
	}

	int num_particles() const
	{
		return particles->num_particles;
	}
};
