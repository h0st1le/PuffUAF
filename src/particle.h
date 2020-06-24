#ifndef PARTICLE_H_
#define PARTICLE_H_

#ifdef LOW_MEMORY
#define PUFFdouble float
#define PUFFncDouble ncFloat
#else
#define PUFFdouble double
#define PUFFncDouble ncDouble
#endif

class Particle {
    
public:
    PUFFdouble   x, y, z;  // 3-dimensional location
    PUFFdouble   size;  // radius in meters
    long   startTime;
    PUFFdouble   mass_fraction;
    bool     grounded; // is ash on the ground?
    bool     exists;   // has ash left the boundary?
#ifdef PUFF_STATISTICS
    PUFFdouble    dif_x, dif_y, dif_z;
    PUFFdouble    adv_x, adv_y, adv_z;
#endif
    
    Particle();
    Particle(double xx, double yy, double zz);
    ~Particle();
    
    Particle & operator=(const Particle &pnt);
    Particle & operator+(Particle &pnt); 
    Particle & operator-(Particle &pnt); 
		bool operator<(Particle pnt) const;
		bool operator<=(Particle pnt) const;
		Particle & operator*(float f);
    

};

#endif 
