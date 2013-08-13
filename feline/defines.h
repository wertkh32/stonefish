#pragma once
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define CLAMP(n,a,b) MIN(MAX(a,n),b)
#define TO_RAD(a) ((a)*(M_PI/180.))
#define INF 10000000

#define MAX_ELEMENTS 100000
#define MAX_NODES 20000
#define FPS 30
#define DT (1.0/FPS)
#define GRAVITY 9.81

#define ALPHA 0.3
#define BETA 0.1
