#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/Random.h"
#include "/home/codeleaded/System/Static/Library/PerlinNoise.h"
#include "/home/codeleaded/System/Static/Library/TransformedView.h"

#include "Fluid.h"

TransformedView tv;
Fluid fluid;

void Setup(AlxWindow* w){
	ResizeAlxFont(25,25);

	tv = TransformedView_New((Vec2){ GetWidth() * 0.05f,GetHeight() * 0.05f });
	fluid = Fluid_New();
}

void Update(AlxWindow* w){
	TransformedView_HandlePanZoom(&tv,window.Strokes,GetMouse());

	if(Stroke(ALX_MOUSE_L).DOWN){
		Vec2 p = TransformedView_ScreenWorldPos(&tv,GetMouse());
		Vector_Push(&fluid.ps,(FluidPoint[]){ FluidPoint_New(p) });
	}

	if(Stroke(ALX_KEY_W).DOWN) H *= 1.01f;
	if(Stroke(ALX_KEY_S).DOWN) H *= 0.99f;
	if(Stroke(ALX_KEY_E).DOWN) MASS *= 1.01f;
	if(Stroke(ALX_KEY_D).DOWN) MASS *= 0.99f;
	if(Stroke(ALX_KEY_R).DOWN) DENSITY *= 1.01f;
	if(Stroke(ALX_KEY_F).DOWN) DENSITY *= 0.99f;
	if(Stroke(ALX_KEY_T).DOWN) K *= 1.01f;
	if(Stroke(ALX_KEY_G).DOWN) K *= 0.99f;
	if(Stroke(ALX_KEY_Z).DOWN) VISCOSITY *= 1.01f;
	if(Stroke(ALX_KEY_H).DOWN) VISCOSITY *= 0.99f;
	if(Stroke(ALX_KEY_U).DOWN) GRAVITY *= 1.01f;
	if(Stroke(ALX_KEY_J).DOWN) GRAVITY *= 0.99f;
	//if(Stroke(ALX_KEY_U).DOWN) MU_WATER *= 1.01f;
	//if(Stroke(ALX_KEY_J).DOWN) MU_WATER *= 0.99f;

	Fluid_Update(&fluid,w->ElapsedTime);

	Clear(BLACK);

	Fluid_Render(&fluid,&tv);

	Vec2 p = TransformedView_ScreenWorldPos(&tv,GetMouse());
	String str = String_Format("H: %f,M: %f,D: %f,K: %f,V: %f,G: %f",H,MASS,DENSITY,K,VISCOSITY,GRAVITY);
	RenderCStrSize(str.Memory,str.size,0.0f,0.0f,WHITE);
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