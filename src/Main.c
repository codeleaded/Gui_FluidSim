#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/Random.h"
#include "/home/codeleaded/System/Static/Library/PerlinNoise.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"

#include "Fluid.h"

TransformedView tv;
Fluid fluid;

void Setup(AlxWindow* w){
	ResizeAlxFont(25,25);

	tv = TransformedView_Make(
		(Vec2){ GetWidth(),GetHeight() },
		(Vec2){ 0.0f,0.0f },
		(Vec2){ 0.05f,0.05f },
		(float)GetWidth() / (float)GetHeight()
	);
	fluid = Fluid_New();

	DENSITY_H = (1.0f * RADIUS);
	DENSITY_WATER = 3000.0f;
	MASS_PARTICLE = 12000.0f;//(DENSITY_WATER * RADIUS_TERM);
	DENSITY_K = 5000.0f;
	SOUND_SPEED = sqrtf(DENSITY_K / DENSITY_WATER);
}

void Update(AlxWindow* w){
	TransformedView_Output(&tv,(Vec2){ GetWidth(),GetHeight() });
	TransformedView_HandlePanZoom(&tv,window.Strokes,GetMouse());
	Vec2 mouse = TransformedView_ScreenWorldPos(&tv,GetMouse());

	if(Stroke(ALX_KEY_Y).DOWN){
		Vector_Push(&fluid.ps,(FluidPoint[]){FluidPoint_New(Vec2_Add(
			mouse,
			(Vec2){ Random_f64_MinMax(-2.0f * RADIUS,2.0f * RADIUS),Random_f64_MinMax(-2.0f * RADIUS,2.0f * RADIUS) }
		)) });
	}

	if(Stroke(ALX_KEY_W).DOWN) DENSITY_H *= 1.01f;
	if(Stroke(ALX_KEY_S).DOWN) DENSITY_H *= 0.99f;
	if(Stroke(ALX_KEY_E).DOWN) MASS_PARTICLE *= 1.01f;
	if(Stroke(ALX_KEY_D).DOWN) MASS_PARTICLE *= 0.99f;
	if(Stroke(ALX_KEY_R).DOWN) DENSITY_WATER *= 1.01f;
	if(Stroke(ALX_KEY_F).DOWN) DENSITY_WATER *= 0.99f;
	if(Stroke(ALX_KEY_T).DOWN) DENSITY_K *= 1.01f;
	if(Stroke(ALX_KEY_G).DOWN) DENSITY_K *= 0.99f;

	//MASS_PARTICLE = (DENSITY_WATER * RADIUS_TERM);
	SOUND_SPEED = sqrtf(DENSITY_K / DENSITY_WATER);

	// if(Stroke(ALX_KEY_Z).DOWN) VISCOSITY *= 1.01f;
	// if(Stroke(ALX_KEY_H).DOWN) VISCOSITY *= 0.99f;
	// if(Stroke(ALX_KEY_U).DOWN) GRAVITY *= 1.01f;
	// if(Stroke(ALX_KEY_J).DOWN) GRAVITY *= 0.99f;
	//if(Stroke(ALX_KEY_U).DOWN) MU_WATER *= 1.01f;
	//if(Stroke(ALX_KEY_J).DOWN) MU_WATER *= 0.99f;

	Fluid_Calc(&fluid);

	if(Stroke(ALX_MOUSE_L).DOWN)
		Fluid_ApplyForce(&fluid,TransformedView_ScreenWorldPos(&tv,GetMouse()),10.0f,5.0f);
	
	if(Stroke(ALX_MOUSE_R).DOWN)
		Fluid_ApplyForce(&fluid,TransformedView_ScreenWorldPos(&tv,GetMouse()),10.0f,-5.0f);

	FluidPoint fp = FluidPoint_New(mouse);
	Fluid_CalcFP(&fluid,&fp);

	//Fluid_Update(&fluid,w->ElapsedTime);
	Fluid_Update(&fluid,0.01f);

	Clear(BLACK);

	Rect r = {
		.p = TransformedView_WorldScreenPos(&tv,(Vec2){ 0.0f,0.0f }),
		.d = TransformedView_WorldScreenLength(&tv,(Vec2){ BORDER_X,BORDER_Y })
	};
	RenderRectWire(r.p.x,r.p.y,r.d.x,r.d.y,WHITE,1.0f);

	Fluid_Render(WINDOW_STD_ARGS,&fluid,&tv);

	String str = String_Format("S: %d, H: %f,M: %f,D: %f,K: %f",fluid.ps.size,DENSITY_H,MASS_PARTICLE,DENSITY_WATER,DENSITY_K);
	RenderCStrSize(str.Memory,str.size,0.0f,0.0f,WHITE);
	String_Free(&str);

	str = String_Format("RO: %f, P: %f",fp.ro,fp.pres);
	RenderCStrSize(str.Memory,str.size,0.0f,GetAlxFont()->CharSizeY,WHITE);
	String_Free(&str);

	// Vec2 p = TransformedView_ScreenWorldPos(&tv,GetMouse());
	// String str = String_Format("P: X: %f, Y: %f",p.x,p.y);
	// RenderCStrSize(str.Memory,str.size,0.0f,0.0f,WHITE);
	// String_Free(&str);
}

void Delete(AlxWindow* w){
	Fluid_Free(&fluid);
}

int main(){
    if(Create("Fluid-Sim",2500,1200,1,1,Setup,Update,Delete))
        Start();
    return 0;
}