#include <SDL2/SDL.h>
#include <vector>
#include <algorithm>

//#define APPROX_SQRT

template<typename T>
struct Vector_T
{
	T values[2];
	#define x values[0]
	#define y values[1]

	Vector_T operator+(const Vector_T& other)
	{
		return {x + other.x, y + other.y};
	}
	Vector_T& operator+=(const Vector_T& other)
	{
		x += other.x;
		y += other.y;
		return *this;
	}
	Vector_T operator-(const Vector_T& other)
	{
		return {x - other.x, y - other.y};
	}
	Vector_T& operator-=(const Vector_T& other)
	{
		x -= other.x;
		y -= other.y;
		return *this;
	}
	Vector_T operator*(const T c) const
	{
		return {x*c, y*c};
	}
	Vector_T& operator*=(const T c)
	{
		x*=c;
		y*=c;
		return *this;
	}
	Vector_T operator/(const T c) const
	{
		return {x/c, y/c};
	}
	Vector_T& operator/=(const T c)
	{
		x/=c;
		y/=c;
		return *this;
	}
	float dot(const Vector_T other) const
	{
		return x*other.x + y*other.y;
	}
	float sqr_norm() const
	{
		return x*x + y*y;
	}
	T& operator[](int i)
	{
		return values[i];
	}
	const T& operator[](int i) const
	{
		return values[i];
	}

	#undef x
	#undef y
};


template<typename T>
std::ostream& operator<<(std::ostream& os, const Vector_T<T> v) 
{
	os << '(' << v[0] << ',' << v[1] << ')';
	return os;
}

using Vector = Vector_T<float>;

#define NUM_ITERATIONS 10
constexpr float gravity = 200;


struct Line
{
	float start_x;
	float start_y;

	float end_x;
	float end_y;
};
struct Particle
{
	Vector position;
	Vector old_position;
	Vector force;
};
struct Bond
{
	int particle_1;
	int particle_2;
	float length;
	float length_squared;
};
struct FixedPoint
{
	Vector point;
	int particle;
};

struct ParticleSim
{
	std::vector<Particle> particles;
	std::vector<Bond> bonds;
	std::vector<FixedPoint> fixed_points;

	ParticleSim(int n)
	{
	}

	void calc_forces()
	{
		for(auto& p : particles)
		{
			p.force = {0, gravity};
		}
	}

	void integrate(float dt)
	{
		for(auto& p : particles)
		{
			Vector& r = p.position;
			Vector& r_old = p.old_position;
			Vector& a = p.force;

			Vector temp = r;
			r += (r-r_old) + a*(dt*dt);
			r_old = temp;
		}
	}

	void apply_bonds()
	{
		for(int j = 0; j<NUM_ITERATIONS; ++j)
		{
			for(const auto& b : bonds)
			{
				Vector& x1 = (particles[b.particle_1]).position;
				Vector& x2 = (particles[b.particle_2]).position;

				// Taylor approx. for sqrt
				Vector delta = x2-x1;

				#ifdef APPROX_SQRT
					delta *= -(b.length_squared/(delta.dot(delta)+b.length_squared)-0.5);
				#else
					float deltalength = sqrt(delta.dot(delta));
					float diff=(deltalength-b.length)/deltalength;
					delta *= 0.5*diff;
				#endif
				x1 += delta;
				x2 -= delta;
			}

			for(const auto& p : fixed_points)
			{
				particles[p.particle].position = p.point;
			}
		}
	}

	void step(float timestep)
	{
		calc_forces();
		integrate(timestep);
		apply_bonds();
	}

	Line bond_as_points(const Bond& b) const
	{
		const Particle& p1 = particles[b.particle_1];
		const Particle& p2 = particles[b.particle_2];
		Vector start = p1.position;
		Vector end   = p2.position;
		return {start[0], start[1], end[0], end[1]};
	}

	void add_particle(Vector position = {0,0}, bool fixed=false)
	{
		Particle p;
		p.position = position;
		p.old_position = position;
		p.force = {0,0};

		particles.push_back(p);
		if(fixed)
		{
			FixedPoint f = {position, (int)particles.size()-1};
			fixed_points.push_back(f);
		}
	}

	void add_bond(int p1_id, int p2_id, float length)
	{
		Bond b = {p1_id, p2_id, length, length*length};
		bonds.push_back(b);
	}

	void fix_point(int particle_id)
	{
		if(!point_is_fixed(particle_id))
		{
			FixedPoint f = {particles[particle_id].position, particle_id};
			fixed_points.push_back(f);
		}
	}

	bool point_is_fixed(int particle_id)
	{
		return std::any_of(fixed_points.begin(),fixed_points.end(),
				[&](FixedPoint f){ return f.particle==particle_id; });
	}
	void remove_fixed_point(int particle_id)
	{
		fixed_points.erase(std::remove_if(fixed_points.begin(), fixed_points.end(), 
			[&](FixedPoint f){ return f.particle==particle_id; }), fixed_points.end());
	}

	int num_particles() const
	{
		return particles.size();
	}
};