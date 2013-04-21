#pragma once
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define CLAMP(n,a,b) MIN(MAX(a,n),b)
#define TO_RAD(a) ((a)*(M_PI/180.))
#define INF 10000000

#define MAX_ELEMENTS 500
#define MAX_NODES 100
#define FPS 60
#define GRAVITY 9.81