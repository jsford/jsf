///usr/bin/tcc -run $0 $@ ; exit

#define JSF_DUBINS_IMPLEMENTATION
#include "jsf_dubins.h"
#include <stdio.h>

static void print_path(DubinsPath *path) {
  printf("%s\n", DubinsPathTypeString[path->type]);
  printf("\tARC0: (%0.3f, %0.3f)\n", path->curvatures[0], path->arclengths[0]);
  printf("\tARC1: (%0.3f, %0.3f)\n", path->curvatures[1], path->arclengths[1]);
  printf("\tARC2: (%0.3f, %0.3f)\n", path->curvatures[2], path->arclengths[2]);
}

int main(int argc, char **argv) {
  (void)argc; (void)argv;

  DubinsState s0;
  s0.x = -2.0;
  s0.y = -1.0;
  s0.h = 0.0;

  DubinsState s1;
  s1.x = 2.0;
  s1.y = 1.0;
  s1.h = M_PI/2.0;

  DubinsPath path;
  if( !jsf_dubins(&s0, &s1, 1.0, &path) ) {
    printf("Failed to generate a dubins path.\n");
    return 1;
  }

  print_path(&path);
  printf("\n");

  double path_length = jsf_dubins_path_length(&path);
  for(double a=0.0; a<=path_length; a+=path_length/10.0) {
    DubinsState sample;
    sample = jsf_dubins_path_sample(&s0, &path, a);

    printf("%0.3f, %0.3f, %0.3f\n", sample.x, sample.y, sample.h);
  }

  return 0;
}
