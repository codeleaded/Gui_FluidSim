#ifndef FLUID_H
#define FLUID_H

#include "/home/codeleaded/System/Static/Library/Circle.h"
#include "/home/codeleaded/System/Static/Library/Random.h"
#include "/home/codeleaded/System/Static/Library/Pixel.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"

/*
#define FLUID_SIZE_X		10.0f
#define FLUID_SIZE_Y		10.0f
#define FLUID_COUNT_X		10
#define FLUID_COUNT_Y		10

#define DENSITY_H					2.0f
#define MASS				3141.5f
#define DENSITY				1000.0f
#define K					2000.0f
#define VISCOSITY			0.1f
#define MU_WATER 			100.0f // 0.001f

float DENSITY_H = 1.0f;
float MASS = 1.0f;
float DENSITY = 1000.0f;
float K = 2000.0f;
float VISCOSITY = 0.1f;
float MU_WATER = 0.001f; // 0.001f
float GRAVITY = 12.0f;

typedef struct FluidPoint {
	Vec2 p;
	Vec2 v;
	Vec2 a;
	Vec2 pp;
	float ps;
	float dy;
	Pixel c;
} FluidPoint;

FluidPoint FluidPoint_New(Vec2 p){
	FluidPoint f;
	f.p = p;
	f.v = (Vec2){ 0.0f,0.0f };
	f.a = (Vec2){ 0.0f,0.0f };
	f.ps = 0.0f;
	f.dy = 0.0f;
	f.c = BLUE;
	return f;
}
void FluidPoint_PredictUpdate(FluidPoint* f,float t){
	//f->v = Vec2_Mulf(f->a,t);
	f->pp = Vec2_Add(f->p,Vec2_Mulf(f->v,t));
	
	if(f->pp.x < 0.0f) 		f->pp.x = 0.0f;
	if(f->pp.x > 10.0f)  	f->pp.x = 10.0f;
	if(f->pp.y < 0.0f) 		f->pp.y = 0.0f;
	if(f->pp.y > 10.0f)  	f->pp.y = 10.0f;
}
void FluidPoint_Update(FluidPoint* f,float t){
	//f->v = Vec2_Mulf(f->a,t);
	f->v = Vec2_Add(f->v,Vec2_Mulf(f->a,t));
	f->p = Vec2_Add(f->p,Vec2_Mulf(f->v,t));
	f->c = Pixel_Mulf(WHITE,F32_Clamp(f->dy / DENSITY,0.0f,1.0f));

	if(f->p.x < 0.0f){
		f->p.x = 0.0f;
		f->v.x *= -1.0f;
	}
	if(f->p.x > 10.0f){
		f->p.x = 10.0f;
		f->v.x *= -1.0f;
	}
	if(f->p.y < 0.0f){
		f->p.y = 0.0f;
		f->v.y *= -1.0f;
	}
	if(f->p.y > 10.0f){
		f->p.y = 10.0f;
		f->v.y *= -1.0f;
	}
}
void FluidPoint_Free(FluidPoint* f){
	
}


typedef struct Fluid {
	Vector ps;// Vector<FluidPoint>
	float pressure;
} Fluid;

float poly6(float r, float h) {
    if (r >= 0 && r <= h) {
        float hr2 = h * h - r * r;
        return (315.0f / (64.0f * F32_PI * pow(h, 9))) * pow(hr2, 3);
    }
    return 0.0f;
}
Vec2 spiky_kernel_gradient(Vec2 r_vec,float h){
	float r = Vec2_Mag(r_vec);
    if(0 < r && r <= h){
		float coeff = -45.0f / (F32_PI * powf(h, 6)) * powf(h - r, 2);
        return Vec2_Mulf(Vec2_Norm(r_vec),coeff);
	}
    return (Vec2){ 0.0f,0.0f };
} 
// float estimate_density_at_point(Vec2 pos,Fluid* f, float h){
// 	float rho = 0.0f;
// 	for(int i = 0;i<f->ps.size;i++){
// 		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
// 		float r = Vec2_Mag(Vec2_Sub(pos,p->p));
//         rho += MASS * poly6_kernel(r, h);
// 	}
//     return rho;
// }
// float pressure_from_density(float rho,float rho0,float k){
// 	return k * (rho - rho0);
// }
// float pressure_at_point(Vec2 pos,Fluid* f,float h,float rho0,float k){
// 	float rho = estimate_density_at_point(pos,f,h);
//     return pressure_from_density(rho,rho0,k);
// }
// Vec2 pressure_force_at_point(Vec2 pos,Fluid* f,float h,float rho0,float k){
// 	Vec2 force = (Vec2){ 0.0f,0.0f };
// 	float rho_x = estimate_density_at_point(pos,f,h);
//     float p_x = pressure_from_density(rho_x,rho0,k);
// 	for(int i = 0;i<f->ps.size;i++){
// 		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
// 		Vec2 r_vec = Vec2_Sub(pos,p->p);
//         float r = Vec2_Mag(r_vec);
//         if(0.0f < r && r <= h){
// 			float p_avg = (p->ps + p_x) / 2.0f;
//             Vec2 gradW = spiky_kernel_gradient(r_vec, h);
//             force = Vec2_Add(force,Vec2_Mulf(gradW,-MASS * p_avg / DENSITY));
// 		}
// 	}
//     return force;
// }

float Fluid_Poly6(float r,float h){
	if(r >= h) return 0.0f;
	float factor = h * h - r * r;
	return (315.0f / (64.0f * F32_PI * F32_Pow(h,9.0f))) * factor * factor * factor;
}
float Fluid_Density(FluidPoint* f1,FluidPoint* f2){
	float dx = f2->p.x - f1->p.x;
	float dy = f2->p.y - f1->p.y;
	float r = F32_Sqrt(dx * dx + dy * dy);
	return MASS * Fluid_Poly6(r,DENSITY_H);
}
float Fluid_Pressure(float dy){
	return K * (dy - DENSITY);
}
float Fluid_Spiky_Gradient(float r,float h){
	if(r==0 || r>=h) return 0.0f;
	return (-45.0f / (F32_PI * F32_Pow(h,6.0f))) * (h - r) * (h - r);
}
float Fluid_Laplacian(float r,float h){
	if(r>=h) return 0.0f;
	return (45.0f / (F32_PI * F32_Pow(h,6.0f))) * (h - r);
}

Vec2 Fluid_Force(FluidPoint* f1,FluidPoint* f2){
	float dx = f2->p.x - f1->p.x;
	float dy = f2->p.y - f1->p.y;
	float r = F32_Sqrt(dx * dx + dy * dy);
	if(r < DENSITY_H && r > 0){
		float grad = Fluid_Spiky_Gradient(r,DENSITY_H);
		float pressure = (f1->ps + f2->ps) / (2.0f * f2->dy);
		return (Vec2){
			.x = -MASS * pressure * grad * (dx / r),
			.y = -MASS * pressure * grad * (dy / r)
		};
	}
	return (Vec2){ 0.0f,0.0f };
}
Vec2 Fluid_ViscosityForce(FluidPoint* f1,FluidPoint* f2){
	float dx = f2->p.x - f1->p.x;
	float dy = f2->p.y - f1->p.y;
	float r = F32_Sqrt(dx * dx + dy * dy);
	if(r < DENSITY_H){
		float lap = Fluid_Laplacian(r,DENSITY_H);
		return (Vec2){
			.x = MU_WATER * MASS * (f2->v.x - f1->v.x) / f2->dy * lap,
			.y = MU_WATER * MASS * (f2->v.y - f1->v.y) / f2->dy * lap
		};
	}
	return (Vec2){ 0.0f,0.0f };
}

float Fluid_Kernel(float r,float h){
	//float volume = F32_PI * F32_Pow(h,8.0f) / 4.0f;
	//float v = F32_Max(0.0f,h * h - r * r);
	//return v * v * v / volume;

	if(r >= h) return 0.0f;
	float v = F32_PI * F32_Pow(h,4.0f) / 6.0f;
	return (h - r) * (h - r) / v;
}
float Fluid_Kernel_Derif(float r,float h){
	// if(r >= h) return 0.0f;
	// float f = h * h - r * r;
	// float s = -24.0f / F32_PI * F32_Pow(h,8.0f);
	// return s * r * f * f;

	if(r >= h) return 0.0f;
	float v = 12.0f / (F32_PI * F32_Pow(h,4.0f));
	return (h - r) * v;
}
float Fluid_Kernel_Density(Fluid* f, int index){
    FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps, index);
    float dy = 0.0f;
    for(int i = 0; i < f->ps.size; i++){
        FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps, i);
        float r = Vec2_Mag(Vec2_Sub(p2->p, p1->p));
        dy += MASS * poly6_kernel(r, DENSITY_H);
    }
    return dy;
}
float Fluid_Kernel_Pressure(float dy){
	return K * (dy - DENSITY);
}
float Fluid_Kernel_Visc(float r,float h){
	float volume = F32_PI * F32_Pow(h,8.0f) / 4.0f;
	float v = F32_Max(0.0f,h * h - r * r);
	return v * v * v / volume;
}

Vec2 Fluid_Kernel_Gradient(Fluid* f,int index){
	FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,index);
	
	Vec2 prop = { 0.0f,0.0f };
	for(int i = 0;i<f->ps.size;i++){
		if(i==index) continue;

		FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,i);
		Vec2 dir = Vec2_Sub(p2->p,p1->p);
		float r = Vec2_Mag(dir);
		if(r==0.0f) dir = (Vec2){ Random_f64_New(),Random_f64_New() };
		dir = Vec2_Norm(dir);
		float slope = Fluid_Kernel_Derif(r,DENSITY_H);
		float dy = p2->dy;
		prop = Vec2_Add(prop,Vec2_Mulf(dir,-p1->ps * slope * MASS / dy));
	}
	return prop;
}
Vec2 Fluid_Kernel_Force(Fluid* f,int index){
	FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,index);
	
	Vec2 pressure = { 0.0f,0.0f };
	for(int i = 0;i<f->ps.size;i++){
		if(i==index) continue;
		
		FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,i);
		
		Vec2 dir = Vec2_Sub(p2->p,p1->p);
		float r = Vec2_Mag(dir);
		if(r==0.0f) dir = (Vec2){ Random_f64_New(),Random_f64_New() };
		dir = Vec2_Norm(dir);
		
		float slope = Fluid_Kernel_Derif(r,DENSITY_H);
		float dy = p2->dy;
		float sharedpress = (Fluid_Kernel_Pressure(p1->dy) + Fluid_Kernel_Pressure(p2->dy)) / 2;
		pressure = Vec2_Add(pressure,Vec2_Mulf(dir,-sharedpress * slope * MASS / dy));
	}
	return pressure;
}
Vec2 Fluid_Kernel_Viscosity(Fluid* f,int index){
	FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,index);
	
	Vec2 pressure = { 0.0f,0.0f };
	for(int i = 0;i<f->ps.size;i++){
		if(i==index) continue;
		
		FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,i);
		
		Vec2 dir = Vec2_Sub(p2->p,p1->p);
		float r = Vec2_Mag(dir);
		if(r==0.0f) dir = (Vec2){ Random_f64_New(),Random_f64_New() };
		dir = Vec2_Norm(dir);
		
		float infl = Fluid_Kernel_Visc(r,DENSITY_H);
		pressure = Vec2_Add(pressure,Vec2_Mulf(Vec2_Sub(p2->v,p1->v),infl));
	}
	return pressure;
}


void compute_density(Fluid* f, float h, float mass) {
    for (int i = 0; i < f->ps.size; i++) {
		FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,i);
        
		float density = 0.0f;
        for (int j = 0; j < f->ps.size; j++) {
			//if(i == j) continue;

			FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,j);
            float r = Vec2_Mag(Vec2_Sub(p1->p,p2->p));
            density += mass * poly6(r, h);
        }
        p1->dy = density;
    }
}
void compute_pressure(Fluid* f, float k, float rest_density) {
    for (int i = 0; i < f->ps.size; i++) {
		FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,i);
        p1->ps = k * (p1->dy - rest_density);
		p1->ps = F32_Clamp(p1->ps,-1000.0f,5000.0f); // Optional
    }
}
void compute_pressure_force(Fluid* f, float h, float mass) {
    for (int i = 0; i < f->ps.size; i++) {
		FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,i);
        Vec2 pressure_force = {0.0f, 0.0f};

        for (int j = 0; j < f->ps.size; j++) {
        	//if (i == j) continue;

			FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,j);
            Vec2 rij = Vec2_Sub(p1->p,p2->p);
            float dist = Vec2_Mag(rij);
            if (dist > 0 && dist <= h) {
                Vec2 gradW = spiky_kernel_gradient(rij, h);
                float pressure_term = (p1->ps + p2->ps) / (2.0f * p2->dy);
                pressure_force = Vec2_Add(pressure_force,Vec2_Mulf(gradW,-mass * pressure_term));
            }
        }

        p1->a = Vec2_Add(p1->a,Vec2_Divf(pressure_force,p1->dy));
    }
}



Fluid Fluid_New(){
	Fluid f;
	f.ps = Vector_New(sizeof(FluidPoint));
	f.pressure = 0.0f;

	for(int i = 0;i<200;i++){
		Vec2 p = { Random_f64_MinMax(0.0f,10.0f),Random_f64_MinMax(0.0f,10.0f) };
		Vector_Push(&f.ps,(FluidPoint[]){ FluidPoint_New(p) });
	}

	return f;
}
// Vec2 Fluid_Nearest(Fluid* f,FluidPoint* fp){
// 	return pressure_force_at_point(fp->p,f,1.0f,1000.0f,10000.0f);
// }
void Fluid_Update(Fluid* f,float t){
	compute_density(f,DENSITY_H,MASS);
    compute_pressure(f,K,DENSITY);
    compute_pressure_force(f,DENSITY_H,MASS);
	
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,i);
		//FluidPoint_PredictUpdate(p1,t);
		//p1->dy = Fluid_Kernel_Density(f,i);
		//p1->ps = Fluid_Kernel_Pressure(p1->dy);
		//printf("dy[%d] = %f\n", i, p1->dy);
	}
	
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p1 = (FluidPoint*)Vector_Get(&f->ps,i);
		
		//p1->a = Fluid_Kernel_Force(f,i);
		//p1->a = Vec2_Add(p1->a,Fluid_Kernel_Force(f,i));
		//p1->a = Vec2_Add(p1->a,Fluid_Kernel_Viscosity(f,i));

		//p1->a.x /= p1->dy;
		//p1->a.y /= p1->dy;

		//p1->a.y += GRAVITY;
		
		// p1->a.x = 0.0f;
		// p1->a.y = 0.0f;
		// for(int j = 0;j<f->ps.size;j++){
		// 	if(i==j) continue;
		// 	FluidPoint* p2 = (FluidPoint*)Vector_Get(&f->ps,j);
			
		// 	Vec2 ff = Fluid_Force(p1,p2);
		// 	p1->a.x += ff.x;
		// 	p1->a.y += ff.y;
			
		// 	Vec2 vf = Fluid_ViscosityForce(p1,p2);
		// 	p1->a.x += vf.x;
		// 	p1->a.y += vf.y;
		// }
		//p1->a = Vec2_Add(p1->a,(Vec2){ 0.0f,p1->dy * -9.81f });
		
		//printf("A: %f,%f\n",p1->a.x,p1->a.y);
		// if(p1->dy > 1e-6f){
		// 	p1->a.x /= p1->dy;
		// 	p1->a.y /= p1->dy;
		// }
		
		//printf("-------\n");
	}
	//exit(0);
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
		FluidPoint_Update(p,t);
		p->a = (Vec2){ 0.0f,0.0f };
	}
}
void Fluid_Render(Fluid* f,TransformedView* tv){
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);

		Vec2 sp = TransformedView_WorldScreenPos(tv,p->p);
		float r = TransformedView_WorldScreenLX(tv,0.1f);
		RenderCircle(sp,r,p->c);
	}
}
void Fluid_Free(Fluid* f){
	Vector_Free(&f->ps);
}
*/

#define BORDER_X    				50.0f
#define BORDER_Y    				50.0f
#define RADIUS      				1.0f
#define RADIUS_TERM      			(RADIUS * RADIUS)

// #define DENSITY_H           		1.0f//(2.0f * RADIUS)
// #define DENSITY_WATER      			1000.0f
// #define MASS_PARTICLE				(DENSITY_WATER * powf(RADIUS * 1.5f, 2.0f))
// #define DENSITY_K       			20000.0f

float DENSITY_H;
float DENSITY_WATER;
float MASS_PARTICLE;
float DENSITY_K;
float SOUND_SPEED;

//#define SOUND_SPEED     			50.0f//sqrtf(DENSITY_K / DENSITY_WATER)
#define DENSITY_KINEMATIC_WATER		1.0E-3f
#define VISCOSITY_ALPHA				1.5f
#define VISCOSITY_BETA				2.5f

#define GAMMA_AIR 					1.4f
#define GAMMA_WATER 				7.0f


typedef struct FluidPoint {
	Vec2 p;
	Vec2 v;
	Vec2 a;
	Vec2 pp;
	Vec2 aex;
	float ro;
	float pres;
	float cs;
	float m;
	Pixel c;
} FluidPoint;

FluidPoint FluidPoint_New(Vec2 p){
	FluidPoint f;
	f.p = p;
	f.v = (Vec2){ 0.0f,0.0f };
	f.a = (Vec2){ 0.0f,15.0f };
	f.pp = (Vec2){ 0.0f,0.0f };
	f.aex = (Vec2){ 0.0f,0.0f };
	f.ro = 0.0f;
	f.pres = 0.0f;
	f.cs = 0.0f;
	f.m = 0.0f;
	f.c = BLUE;
	return f;
}
void FluidPoint_Update(FluidPoint* f,float t){
	//f->v = Vec2_Mulf(f->a,t);
	//f->v = Vec2_Add(f->v,Vec2_Mulf(f->a,t));
	//f->v = Vec2_Add(f->v,Vec2_Mulf(f->pp,t));
	//f->v = Vec2_Mulf(Vec2_Add(f->pp,f->a),t);
	
	if(isnan(f->pp.x) || isnan(f->pp.y)) f->pp = (Vec2){ 0.0f,0.0f };
	f->v = Vec2_Add(f->v,Vec2_Mulf(Vec2_Add(f->pp,Vec2_Add(f->a,f->aex)),t));
	f->aex = (Vec2){ 0.0f,0.0f };

    //f->v = Vec2_Mulf(f->v,0.9999f);
	f->p = Vec2_Add(f->p,Vec2_Mulf(f->v,t));

	if(f->p.x < 0.0f){
		f->p.x = 0.0f;
		f->v.x *= -1.0f;
	}
	if(f->p.x > BORDER_X){
		f->p.x = BORDER_X;
		f->v.x *= -1.0f;
	}
	if(f->p.y < 0.0f){
		f->p.y = 0.0f;
		f->v.y *= -1.0f;
	}
	if(f->p.y > BORDER_Y){
		f->p.y = BORDER_Y;
		f->v.y *= -0.9f;
	}

	if(isnan(f->p.x) || isnan(f->p.y) || isnan(f->v.x) || isnan(f->v.y)){
	    f->v = (Vec2){0.0f, 0.0f};
	    f->pp = (Vec2){0.0f, 0.0f};
	    if (isnan(f->p.x) || isnan(f->p.y))
			f->p = (Vec2){ F32_Clamp(f->p.x,0.0f,BORDER_X),F32_Clamp(f->p.y,0.0f,BORDER_Y) };
	}
	
	//f->c = Pixel_Mulf(WHITE,F32_Clamp(1.0f,0.0f,1.0f));
	f->c = Pixel_toRGBA(
		F32_Sin_Sq(f->ro / DENSITY_WATER) * 0.9f + 0.1f,
		F32_Sin_Sq(f->ro / DENSITY_WATER) * 0.7f + 0.3f,
		F32_Sin_Sq(f->ro / DENSITY_WATER) * 0.8f + 0.2f,
		1.0f
	);
}

float Gradient_Function(float a){
    const float border = 0.25f;
    
    float ret;
    if(a > border)  ret = 1.0f / (a * a);
    else            ret = 1.0f / (border * border);
    return ret * 0.5f;

    //const float t = 10.0f;
    //const float z = 1.0f;
    //return t * (-1.0f / (1.0f + powf(2.718f,-z * (a - 1.0f))) + 1.0f);
}
float Fluid_Kernel_W(float r,float h){
    const float q = r / h;
    const float sigma = 10.0f / (7.0f * F32_PI * h * h);
    if(q >= 0.0f && q < 1.0f){
        return sigma * (1.0f - 1.5f * q * q + 0.75f * q * q * q);
    } else if(q < 2.0f){
        const float v = (2.0f - q);
        return sigma * (0.25f * v * v * v);
    } else {
        return 0.0f;
    }
}
float Fluid_Kernel_dWdr(float r, float h){
    if(r <= 0.0f) r = 1e-12f; // avoid div by zero in direction calc
    const float q = r / h;
    const float sigma = 10.0f / (7.0f * F32_PI * h * h);
    if(q >= 0.0f && q < 1.0f){
        // d/dq of (1 - 1.5 q^2 + 0.75 q^3) = -3 q + 2.25 q^2
        const float dW_dq = (-3.0f * q + 2.25f * q * q);
        return sigma * (dW_dq / h);
    } else if(q < 2.0f){
        // d/dq of 0.25*(2-q)^3 = 0.25 * -3*(2-q)^2
        const float v = (2.0f - q);
        const float dW_dq = -0.75f * v * v;
        return sigma * (dW_dq / h);
    } else {
        return 0.0f;
    }
}
float Fluid_Kernel_P(float k,float roi,float ro0){
    return k * (roi - ro0);
}
float Fluid_Kernel_CMF_A(FluidPoint* p1, FluidPoint* p2){
    const Vec2 rij = Vec2_Sub(p1->p, p2->p);
    const Vec2 vij = Vec2_Sub(p1->v, p2->v);
    const float rij2 = Vec2_Mag2(rij);
    const float vijdotr = Vec2_Dot(vij, rij);
    const float divzero2 = 0.01f * RADIUS * RADIUS;
    const float cbar = 0.5f * (p1->cs + p2->cs);
    const float rhobar = 0.5f * (p1->ro + p2->ro);
    if(vijdotr < 0.0f){
        const float mu = (RADIUS * vijdotr) / (rij2 + divzero2);
        return (-VISCOSITY_ALPHA * cbar * mu + VISCOSITY_BETA * mu * mu) / (rhobar);
    } else {
        return 0.0f;
    }
}
float Fluid_Kernel_CMF_Phy(FluidPoint* p1,FluidPoint* p2){
	const Vec2 rij = Vec2_Sub(p1->p, p2->p);
    const Vec2 vij = Vec2_Sub(p1->v, p2->v);
    const float rij2 = Vec2_Mag2(rij);
    const float vijdotr = Vec2_Dot(vij,rij);
    const float divzero2 = 0.01f * DENSITY_H * DENSITY_H;
	const float mi = (DENSITY_H * vijdotr) / (rij2 + divzero2);
    const float nu = mi / p1->ro;
    return (2.0f * nu * vijdotr) / ((p1->ro + p2->ro) * (rij2 + divzero2));
}


typedef struct Fluid {
	Vector ps;// Vector<FluidPoint>
	float pressure;
} Fluid;

Fluid Fluid_New(){
	Fluid f;
	f.ps = Vector_New(sizeof(FluidPoint));
	f.pressure = 0.0f;

	for(float i = 0.0f;i<BORDER_Y;i+=RADIUS*2.0f){
		for(float j = 0.0f;j<BORDER_X;j+=RADIUS*2.0f){
			Vec2 p = { j,i };
			Vector_Push(&f.ps,(FluidPoint[]){ FluidPoint_New(
				Vec2_Add(p,(Vec2){ Random_f64_MinMax(-0.5f*RADIUS,0.5f*RADIUS),Random_f64_MinMax(-0.5f*RADIUS,0.5f*RADIUS) }))});
		}
	}

	return f;
}
void Fluid_ApplyForce(Fluid* f,Vec2 pos,float radius,float force){
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);

        const Vec2 dir = Vec2_Sub(pos,p->p);
        const float mag2 = Vec2_Mag2(dir);
	    if(mag2 < radius * radius){
            const Vec2 norm = Vec2_Norm(dir);
            p->aex = Vec2_Mulf(norm,F32_Max(sqrtf(mag2),0.1f) * force);
        }
	}
}
void Fluid_CalcFP(Fluid* f,FluidPoint* p){
	p->ro = 0.0f;

    for(int j = 0;j<f->ps.size;j++){
		FluidPoint* o = (FluidPoint*)Vector_Get(&f->ps,j);
		if(p == o) continue;
		
		const Vec2 dir = Vec2_Sub(o->p,p->p);
        const float mag = F32_Max(Vec2_Mag(dir),0.01f);
		
        const float w = Fluid_Kernel_W(mag,DENSITY_H);
		p->ro += o->m * w;
	}
	p->ro = F32_Clamp(p->ro,DENSITY_WATER * 0.8f,DENSITY_WATER * 5.0f);
	
	//p->pres = Fluid_Kernel_P(DENSITY_K,p->ro,DENSITY_WATER);
	//p->pres = fmaxf(p->pres,0.0f);
	
	float rho_ratio = fmaxf(p->ro / DENSITY_WATER,0.1f);
	p->pres = (SOUND_SPEED*SOUND_SPEED*DENSITY_WATER/GAMMA_WATER)*(powf(rho_ratio,GAMMA_WATER)-1.0f);
	//p->pres = (SOUND_SPEED * SOUND_SPEED * DENSITY_WATER) * (rho_ratio - 1.0f);
	
	p->pres = F32_Clamp(p->pres,0.0f,DENSITY_WATER * 5.0f * DENSITY_K);
	
	p->cs = SOUND_SPEED;//sqrtf(p->pres / p->ro);
}
void Fluid_Calc(Fluid* f){
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
		p->m = MASS_PARTICLE;
		Fluid_CalcFP(f,p);
	}
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
        p->pp = (Vec2){ 0.0f,0.0f };

        for(int j = 0;j<f->ps.size;j++){
	    	if(i == j) continue;
			FluidPoint* o = (FluidPoint*)Vector_Get(&f->ps,j);

			// const Vec2 dir = Vec2_Sub(o->p,p->p);
			// const Vec2 grad = Vec2_Norm(dir);
            // const float mag = Vec2_Mag(dir);

			// p->pp = Vec2_Add(
			// 	p->pp,
			// 	Vec2_Mulf(
			// 		grad,
			// 		(
			// 			p->pres / (p->ro * p->ro)
			// 			+ o->pres / (o->ro * o->ro)
			// 			//+ Fluid_Kernel_CMF_A(p,o)
			// 		) * o->m
			// 	)
			// );

			const Vec2 rij = Vec2_Sub(p->p,o->p); // r_i - r_j
        	const float rij_len = F32_Max(Vec2_Mag(rij),0.01f);
        	if(rij_len > 2.0f * DENSITY_H) continue;

        	const Vec2 rhat = Vec2_Mulf(rij,1.0f / rij_len);
        	const float dWdr = Fluid_Kernel_dWdr(rij_len,DENSITY_H);
        	const float pres_term = (p->pres / (p->ro * p->ro) + o->pres / (o->ro * o->ro)) * o->m;
        	const float visc = Fluid_Kernel_CMF_A(p,o) * o->m;
        	const float scalar = (pres_term + visc) * dWdr;

			p->pp = Vec2_Add(p->pp,Vec2_Mulf(rhat,-scalar));
        	//p->pp = Vec2_Add(p->pp,Vec2_Mulf(rhat,-0.5f * scalar));
			//o->pp = Vec2_Add(o->pp,Vec2_Mulf(rhat,0.5f * scalar));

			//printf("%f = %f,%f,%f\n",pres_term,p->pres,p->ro,o->m);
			//printf("%f = %f * %f\n",visc,Fluid_Kernel_CMF_A(p,o),o->m);
			//printf("%f = (%f + %f) * %f\n",scalar,pres_term,visc,dWdr);
	    }
        //p->pp = Vec2_Neg(p->pp);
	}
}
void Fluid_Update(Fluid* f,float t){
    for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
        FluidPoint_Update(p,t);
	}
}
void Fluid_Render(unsigned int *Target,int Target_Width,int Target_Height,Fluid* f,TransformedView* tv){
	for(int i = 0;i<f->ps.size;i++){
		FluidPoint* p = (FluidPoint*)Vector_Get(&f->ps,i);
		if(isnan(p->p.x) || isnan(p->p.y)) continue;

		Vec2 sp = TransformedView_WorldScreenPos(tv,p->p);
		float r = TransformedView_WorldScreenLX(tv,RADIUS);
		Circle_RenderX(Target,Target_Width,Target_Height,sp,r,p->c);
	}
}
void Fluid_Free(Fluid* f){
	Vector_Free(&f->ps);
}

#endif