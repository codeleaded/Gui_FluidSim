#ifndef FLUID_H
#define FLUID_H

#include "/home/codeleaded/System/Static/Library/Circle.h"
#include "/home/codeleaded/System/Static/Library/Random.h"
#include "/home/codeleaded/System/Static/Library/Pixel.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"

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

#define GRID_X						5.0f
#define GRID_Y						5.0f

#ifdef FLUID_ITER

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
void FluidPoint_Update(FluidPoint* p,float t){
	//p->v = Vec2_Mulf(p->a,t);
	//p->v = Vec2_Add(p->v,Vec2_Mulf(p->a,t));
	//p->v = Vec2_Add(p->v,Vec2_Mulf(p->pp,t));
	//p->v = Vec2_Mulf(Vec2_Add(p->pp,p->a),t);
	
	if(isnan(p->pp.x) || isnan(p->pp.y)) p->pp = (Vec2){ 0.0f,0.0f };
	p->v = Vec2_Add(p->v,Vec2_Mulf(Vec2_Add(p->pp,Vec2_Add(p->a,p->aex)),t));
	p->aex = (Vec2){ 0.0f,0.0f };

    //p->v = Vec2_Mulf(p->v,0.9999f);
	p->p = Vec2_Add(p->p,Vec2_Mulf(p->v,t));

	if(p->p.x < 0.0f){
		p->p.x = 0.0f;
		p->v.x *= -0.9f;
	}
	if(p->p.x > BORDER_X){
		p->p.x = BORDER_X;
		p->v.x *= -0.9f;
	}
	if(p->p.y < 0.0f){
		p->p.y = 0.0f;
		p->v.y *= -0.9f;
	}
	if(p->p.y > BORDER_Y){
		p->p.y = BORDER_Y;
		p->v.y *= -0.9f;
	}

	if(isnan(p->p.x) || isnan(p->p.y) || isnan(p->v.x) || isnan(p->v.y)){
	    p->v = (Vec2){0.0f, 0.0f};
	    p->pp = (Vec2){0.0f, 0.0f};
	    if (isnan(p->p.x) || isnan(p->p.y))
			p->p = (Vec2){ F32_Clamp(p->p.x,0.0f,BORDER_X),F32_Clamp(p->p.y,0.0f,BORDER_Y) };
	}
	
	//p->c = Pixel_Mulf(WHITE,F32_Clamp(1.0f,0.0f,1.0f));
	p->c = Pixel_toRGBA(
		F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.9f + 0.1f,
		F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.7f + 0.3f,
		F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.8f + 0.2f,
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
void Fluid_Insert(Fluid* f,FluidPoint* p){
	Vector_Push(&f->ps,p);
}
int Fluid_Size(Fluid* f){
	return f->ps.size;	
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
        const float mag = F32_Max(Vec2_Mag(dir),0.1f);
		if(mag > 2.0f * DENSITY_H) continue;
		
        const float w = Fluid_Kernel_W(mag,DENSITY_H);
		p->ro += o->m * w;
	}
	p->ro = F32_Clamp(p->ro,DENSITY_WATER * 0.8f,DENSITY_WATER * 3.0f);
	
	float rho_ratio = fmaxf(p->ro / DENSITY_WATER,0.1f);
	p->pres = (SOUND_SPEED*SOUND_SPEED*DENSITY_WATER/GAMMA_WATER)*(powf(rho_ratio,GAMMA_WATER)-1.0f);
	
	p->pres = F32_Clamp(p->pres,0.0f,DENSITY_WATER * 3.0f * DENSITY_K);
	
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

			const Vec2 rij = Vec2_Sub(p->p,o->p); // r_i - r_j
        	const float rij_len = F32_Max(Vec2_Mag(rij),0.1f);
        	if(rij_len > 2.0f * DENSITY_H) continue;

        	const Vec2 rhat = Vec2_Norm(rij);
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

#else

Vec2 Fluid_Pos(Vec2 p){
	return (Vec2){
		floorf(p.x / GRID_X),
		floorf(p.y / GRID_Y),
	};
}
Vec2 Fluid_toPos(Vec2 p){
	return (Vec2){
		p.x * GRID_X,
		p.y * GRID_Y,
	};
}

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
	f.a = (Vec2){ 0.0f,0.0f };
	f.pp = (Vec2){ 0.0f,0.0f };
	f.aex = (Vec2){ 0.0f,0.0f };
	f.ro = 0.0f;
	f.pres = 0.0f;
	f.cs = 0.0f;
	f.m = 0.0f;
	f.c = BLUE;
	return f;
}
void FluidPoint_Update(FluidPoint* p,int cell_x,int cell_y,float t){
	//p->v = Vec2_Mulf(p->a,t);
	//p->v = Vec2_Add(p->v,Vec2_Mulf(p->a,t));
	//p->v = Vec2_Add(p->v,Vec2_Mulf(p->pp,t));
	//p->v = Vec2_Mulf(Vec2_Add(p->pp,p->a),t);
	
	if(isnan(p->pp.x) || isnan(p->pp.y)) p->pp = (Vec2){ 0.0f,0.0f };
	p->v = Vec2_Add(p->v,Vec2_Mulf(Vec2_Add(p->pp,Vec2_Add(p->a,p->aex)),t));
	p->aex = (Vec2){ 0.0f,0.0f };

    //p->v = Vec2_Mulf(p->v,0.9999f);
	p->p = Vec2_Add(p->p,Vec2_Mulf(p->v,t));

	if(p->p.x < 0.0f){
		p->p.x = 0.0f;
		p->v.x *= -0.9f;
	}
	if(p->p.x > BORDER_X){
		p->p.x = BORDER_X;
		p->v.x *= -0.9f;
	}
	if(p->p.y < 0.0f){
		p->p.y = 0.0f;
		p->v.y *= -0.9f;
	}
	if(p->p.y > BORDER_Y){
		p->p.y = BORDER_Y;
		p->v.y *= -0.9f;
	}

	if(isnan(p->p.x) || isnan(p->p.y) || isnan(p->v.x) || isnan(p->v.y)){
	    p->v = (Vec2){0.0f, 0.0f};
	    p->pp = (Vec2){0.0f, 0.0f};
	    if (isnan(p->p.x) || isnan(p->p.y))
			p->p = (Vec2){ F32_Clamp(p->p.x,0.0f,BORDER_X),F32_Clamp(p->p.y,0.0f,BORDER_Y) };
	}
	
	//p->c = Pixel_Mulf(WHITE,F32_Clamp(1.0f,0.0f,1.0f));
	//p->c = Pixel_toRGBA(
	//	F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.9f + 0.1f,
	//	F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.7f + 0.3f,
	//	F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.8f + 0.2f,
	//	1.0f
	//);

	// Vec2 pos = Fluid_Pos(p->p);
	// unsigned int posx = (unsigned int)pos.x;
	// unsigned int posy = (unsigned int)pos.y;
	// float fposx = (float)pos.x / (BORDER_X / GRID_X);
	// float fposy = (float)pos.y / (BORDER_Y / GRID_Y);
	
	float fposx = (float)cell_x / (BORDER_X / GRID_X);
	float fposy = (float)cell_y / (BORDER_Y / GRID_Y);

	p->c = Pixel_toRGBA(
		(F32_Sin_Sq(fposx * F32_PI2) + 1.0f) * (0.5f * 0.9f) + 0.1f,
		(F32_Sin_Sq(fposy * F32_PI2) + 1.0f) * (0.5f * 0.9f) + 0.1f,
		F32_Sin_Sq(p->ro / DENSITY_WATER) * 0.8f + 0.2f,
		1.0f
	);
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
	Vector* points;
	unsigned int width;
	unsigned int height;
} Fluid;

void Fluid_Insert(Fluid* f,FluidPoint* p){
	Vec2 pos = Fluid_Pos(p->p);
	unsigned int posx = (unsigned int)pos.x;
	unsigned int posy = (unsigned int)pos.y;
	
	if(posx < f->width && posy < f->height){
		Vector* points = f->points + posy * f->width + posx;
		Vector_Push(points,p);
	}
}
int Fluid_Size(Fluid* f){
	int size = 0;
	for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;
		size += points->size;
	}
	return size;
}
Fluid Fluid_New(){
	Fluid f;
	Vec2 dim = Fluid_Pos((Vec2){ BORDER_X,BORDER_Y });
	f.width = (unsigned int)(dim.x + 1.0f);
	f.height = (unsigned int)(dim.y + 1.0f);
	f.points = (Vector*)malloc(sizeof(Vector) * f.width * f.height);

	for(int i = 0;i<f.width * f.height;i++){
		f.points[i] = Vector_New(sizeof(FluidPoint));
	}

	for(float i = 0.0f;i<BORDER_Y;i+=RADIUS*2.0f){
		for(float j = 0.0f;j<BORDER_X;j+=RADIUS*2.0f){
			Vec2 p = { j,i };
			Fluid_Insert(&f,(FluidPoint[]){ FluidPoint_New(
				Vec2_Add(
					p,
					(Vec2){ Random_f64_MinMax(-0.5f*RADIUS,0.5f*RADIUS),Random_f64_MinMax(-0.5f*RADIUS,0.5f*RADIUS) }
				))}
			);
		}
	}

	return f;
}
void Fluid_ApplyForce(Fluid* f,Vec2 pos,float radius,float force){
	Vec2 min_p = Fluid_Pos(Vec2_Subf(pos,radius));
	Vec2 max_p = Fluid_Pos(Vec2_Addf(pos,radius));

	const int min_x = (int)F32_Clamp(min_p.x,0.0f,f->width);
	const int min_y = (int)F32_Clamp(min_p.y,0.0f,f->height);
	const int max_x = (int)F32_Clamp(max_p.x,0.0f,f->width);
	const int max_y = (int)F32_Clamp(max_p.y,0.0f,f->height);
	
	for(int y = min_y;y<max_y;y++){
		for(int x = min_x;x<max_x;x++){
			Vector* points = f->points + y * f->width + x;
			
			for(int j = 0;j<points->size;j++){
				FluidPoint* p = (FluidPoint*)Vector_Get(points,j);

				const Vec2 dir = Vec2_Sub(pos,p->p);
    	    	const float mag2 = Vec2_Mag2(dir);
		    	if(mag2 < radius * radius){
    	    	    const Vec2 norm = Vec2_Norm(dir);
    	    	    p->aex = Vec2_Mulf(norm,F32_Max(sqrtf(mag2),0.1f) * force);
    	    	}
			}
		}
	}
}
void Fluid_CalcFP(Fluid* f,FluidPoint* p){
	p->ro = 0.0f;

	Vec2 min_p = Fluid_Pos(Vec2_Subf(p->p,DENSITY_H));
	Vec2 max_p = Fluid_Pos(Vec2_Addf(p->p,DENSITY_H));
	const int min_x = (int)F32_Clamp(min_p.x,0.0f,f->width);
	const int min_y = (int)F32_Clamp(min_p.y,0.0f,f->height);
	const int max_x = (int)F32_Clamp(max_p.x + 1.0f,0.0f,f->width);
	const int max_y = (int)F32_Clamp(max_p.y + 1.0f,0.0f,f->height);

	//printf("--> %f,%f: %d,%d -> %d,%d\n",p->p.x,p->p.y,min_x,min_y,max_x,max_y);

	for(int y = min_y;y<max_y;y++){
		for(int x = min_x;x<max_x;x++){
			Vector* points = f->points + y * f->width + x;
			
			for(int j = 0;j<points->size;j++){
				FluidPoint* o = (FluidPoint*)Vector_Get(points,j);

				if(p == o) continue;
		
				const Vec2 dir = Vec2_Sub(o->p,p->p);
        		const float mag = F32_Max(Vec2_Mag(dir),0.1f);
        		if(mag > 2.0f * DENSITY_H) continue;

        		const float w = Fluid_Kernel_W(mag,DENSITY_H);
				p->ro += o->m * w;
			}
		}
	}

	//printf("%f - ",p->ro);
	p->ro = F32_Clamp(p->ro,DENSITY_WATER * 0.8f,DENSITY_WATER * 3.0f);
	
	float rho_ratio = fmaxf(p->ro / DENSITY_WATER,0.1f);
	p->pres = (SOUND_SPEED * SOUND_SPEED * DENSITY_WATER / GAMMA_WATER) * (powf(rho_ratio,GAMMA_WATER) - 1.0f);
	//printf("%f %f %f\n",p->ro,p->pres,rho_ratio);
	
	p->pres = F32_Clamp(p->pres,0.0f,DENSITY_WATER * 3.0f * DENSITY_K);
	
	p->cs = SOUND_SPEED;//sqrtf(p->pres / p->ro);
}
void Fluid_Calc(Fluid* f){
	for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;
		
		for(int j = 0;j<points->size;j++){
			FluidPoint* p = (FluidPoint*)Vector_Get(points,j);
			p->m = MASS_PARTICLE;
			Fluid_CalcFP(f,p);
		}
	}
	
	for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;

		for(int j = 0;j<points->size;j++){
			FluidPoint* p = (FluidPoint*)Vector_Get(points,j);
			
			Vec2 min_p = Fluid_Pos(Vec2_Subf(p->p,DENSITY_H));
			Vec2 max_p = Fluid_Pos(Vec2_Addf(p->p,DENSITY_H));
			const int min_x = (int)F32_Clamp(min_p.x,0.0f,f->width);
			const int min_y = (int)F32_Clamp(min_p.y,0.0f,f->height);
			const int max_x = (int)F32_Clamp(max_p.x + 1.0f,0.0f,f->width);
			const int max_y = (int)F32_Clamp(max_p.y + 1.0f,0.0f,f->height);
			
			for(int y = min_y;y<max_y;y++){
				for(int x = min_x;x<max_x;x++){
					Vector* cell_points = f->points + y * f->width + x;

					for(int k = 0;k<cell_points->size;k++){
						FluidPoint* o = (FluidPoint*)Vector_Get(cell_points,k);

						if(p == o) continue;
					
						const Vec2 rij = Vec2_Sub(p->p,o->p); // r_i - r_j
        				const float rij_len = F32_Max(Vec2_Mag(rij),0.1f);
        				if(rij_len > 2.0f * DENSITY_H) continue;

        				const Vec2 rhat = Vec2_Norm(rij);
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
				}
			}
		}
	}
}
void Fluid_Update(Fluid* f,float t){
    for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;
		
		for(int j = 0;j<points->size;j++){
			FluidPoint* p = (FluidPoint*)Vector_Get(points,j);
			FluidPoint_Update(p,i % f->width,i / f->width,t);
		}
	}
	for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;
		
		for(int j = points->size - 1;j>=0;j--){
			FluidPoint* p = (FluidPoint*)Vector_Get(points,j);

			Vec2 pos = Fluid_Pos(p->p);
			unsigned int posx = (unsigned int)pos.x;
			unsigned int posy = (unsigned int)pos.y;
			unsigned int index = posy * f->width + posx;

			if(index != i){
				//printf("[%d,%d]: %f,%f -> %d,%d(%d)\n",i,j,p->p.x,p->p.y,posx,posy,index);
				Vector* new_points = f->points + index;
				Vector_Push(new_points,p);
				Vector_Remove(points,j);
			}
		}
	}
}
void Fluid_Render(unsigned int *Target,int Target_Width,int Target_Height,Fluid* f,TransformedView* tv){
	for(int i = 0;i<f->width * f->height;i++){
		Vector* points = f->points + i;
		
		for(int j = 0;j<points->size;j++){
			FluidPoint* p = (FluidPoint*)Vector_Get(points,j);
			if(isnan(p->p.x) || isnan(p->p.y)) continue;

			Vec2 sp = TransformedView_WorldScreenPos(tv,p->p);
			float r = TransformedView_WorldScreenLX(tv,RADIUS);
			Circle_RenderX(Target,Target_Width,Target_Height,sp,r,p->c);

			sp = TransformedView_WorldScreenPos(tv,p->p);
			r = TransformedView_WorldScreenLX(tv,DENSITY_H);
			Circle_RenderXWire(Target,Target_Width,Target_Height,sp,r,YELLOW,1.0f);

			Vec2 min_p = Fluid_Pos(Vec2_Subf(p->p,DENSITY_H));
			Vec2 max_p = Fluid_Pos(Vec2_Addf(p->p,DENSITY_H));
			const int min_x = (int)F32_Clamp(min_p.x,0.0f,f->width);
			const int min_y = (int)F32_Clamp(min_p.y,0.0f,f->height);
			const int max_x = (int)F32_Clamp(max_p.x + 1.0f,0.0f,f->width);
			const int max_y = (int)F32_Clamp(max_p.y + 1.0f,0.0f,f->height);
			Vec2 cp = TransformedView_WorldScreenPos(tv,Fluid_toPos((Vec2){ min_x,min_y }));
			Vec2 ct = TransformedView_WorldScreenPos(tv,Fluid_toPos((Vec2){ max_x,max_y }));
			Rect_RenderXWire(Target,Target_Width,Target_Height,cp,Vec2_Sub(ct,cp),GREEN,1.0f);
		}
	}
}
void Fluid_Free(Fluid* f){
	if(f->points){
		for(int i = 0;i<f->width * f->height;i++){
			Vector* points = f->points + i;
			Vector_Free(points);
		}
		free(f->points);
		f->points = NULL;
	}
	f->width = 0U;
	f->height = 0U;
}

#endif

#endif