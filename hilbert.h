#ifndef WAVE_FUNCTION_H
#define WAVE_FUNCTION_H
#include "lanczos.h"

void build_two_body_jumps(int l);
void generate_wave_function();
void generate_interaction_file(int l); 
#endif
