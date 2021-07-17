/* jsf_dubins.h - v0.0 - public domain dubins path generator - http://
  
  Do this:
      #define JSF_DUBINS_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the implementation.

   // i.e. it should look like this:
   #include ...
   #include ...
   #include ...
   #define JSF_DUBINS_IMPLEMENTATION
   #include "jsf_dubins.h" 

   AUTHOR - Jordan Ford <jsford94@gmail.com>
   LICENSE - See end of file for license information.
*/

#ifndef JSF_DUBINS_H
#define JSF_DUBINS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

//////////////////////////////////////////////////////////////////
//
// PRIMARY API - generate and manipulate dubins paths
//

typedef enum DubinsPathType {
  DUBINS_LRL = 0,
  DUBINS_RLR = 1,
  DUBINS_RSR = 2,
  DUBINS_RSL = 3,
  DUBINS_LSR = 4,
  DUBINS_LSL = 5
} DubinsPathType;

static const char *DubinsPathTypeString[] = {"LRL", "RLR", "RSR",
                                             "RSL", "LSR", "LSL"};

typedef struct DubinsState {
  double x;
  double y;
  double h;
} DubinsState;

typedef struct DubinsPath {
  double curvatures[3];
  double arclengths[3];
  DubinsPathType type;
} DubinsPath;

bool jsf_dubins(const DubinsState *initial_state,
                const DubinsState *final_state, double radius,
                DubinsPath *path);

double jsf_dubins_path_length(const DubinsPath *path);

DubinsState jsf_dubins_path_sample(const DubinsState *initial_state,
                                   const DubinsPath *path,
                                   const double arclength);
//
//
////  end header file  /////////////////////////////////////////////
#endif // JSF_DUBINS_H

#ifdef JSF_DUBINS_IMPLEMENTATION

#include <math.h>
#include <stdlib.h>

#define DUBINS_EPSILON 1e-8

typedef enum DubinsDirection {
  DUBINS_LEFT = 0,
  DUBINS_RIGHT = 1
} DubinsDirection;

typedef struct DubinsCircle {
  double x;
  double y;
  double r;
  DubinsDirection direction;
} DubinsCircle;

typedef struct DubinsLine {
  double x0, y0;
  double x1, y1;
} DubinsLine;

static double jsf__norm2(double x, double y) { return x * x + y * y; }

static double jsf__norm(double x, double y) { return sqrt(jsf__norm2(x, y)); }

static void jsf__state_tangent_circle(const DubinsState *state, double radius,
                                 const DubinsDirection direction,
                                 DubinsCircle *tangent_circle) {
  double alpha = (direction == DUBINS_LEFT) ? M_PI / 2.0 : -M_PI / 2.0;
  tangent_circle->x = state->x + radius * cos(alpha + state->h);
  tangent_circle->y = state->y + radius * sin(alpha + state->h);
  tangent_circle->r = radius;
  tangent_circle->direction = direction;
}

static bool jsf__tangent_line(const DubinsCircle *c0, const DubinsCircle *c1,
                              DubinsLine *tangent_line) {
  // This method requires both circles to have the same radius.
  if (fabs(c0->r - c1->r) > DUBINS_EPSILON) {
    return false;
  }

  // Compute the distance from c0 to c1.
  double c0c1_x = c1->x - c0->x;
  double c0c1_y = c1->y - c0->y;
  double D = jsf__norm(c0c1_x, c0c1_y);

  // Normalize the vector from the center of c0 to c1.
  double c0c1_hat_x = c0c1_x / D;
  double c0c1_hat_y = c0c1_y / D;

  if (c0->direction == c1->direction) {
    // Construct an exterior tangent line.

    // If the circles are identical, there is no unique tangent line.
    if (fabs(D) < DUBINS_EPSILON) {
      return false;
    }

    // Connect two LEFT oriented circles by rotating clockwise 90 degrees.
    // Connect two RIGHT oriented circles by rotating counter-clockwise 90
    // degrees.
    double rot_x, rot_y;
    if (c0->direction == DUBINS_LEFT) {
      rot_x = c0c1_hat_y;
      rot_y = -c0c1_hat_x;
    } else {
      rot_x = -c0c1_hat_y;
      rot_y = c0c1_hat_x;
    }

    // Find the tangent points on c0 and c1 that create the tangent line.
    tangent_line->x0 = c0->x + c0->r * rot_x;
    tangent_line->y0 = c0->y + c0->r * rot_y;
    tangent_line->x1 = tangent_line->x0 + c0c1_x;
    tangent_line->y1 = tangent_line->y0 + c0c1_y;
  } else {
    // Construct an interior tangent line.

    // If the circles overlap, there is no interior tangent line.
    if (fabs(D) < 2 * c0->r) {
      return false;
    }

    // Compute alpha, the angle between the c0c1 vector and the tangent points.
    double alpha = acos(c0->r / (0.5 * D));

    // Rotate by alpha. CW for LEFT to RIGHT and CCW for RIGHT to LEFT.
    double rot_x, rot_y;
    double c = cos(alpha);
    double s = sin(alpha);
    if (c0->direction == DUBINS_LEFT) {
      rot_x = c * c0c1_hat_x + s * c0c1_hat_y;
      rot_y = -s * c0c1_hat_x + c * c0c1_hat_y;
    } else {
      rot_x = c * c0c1_hat_x - s * c0c1_hat_y;
      rot_y = s * c0c1_hat_x + c * c0c1_hat_y;
    }

    tangent_line->x0 = c0->x + c0->r * rot_x;
    tangent_line->y0 = c0->y + c0->r * rot_y;

    tangent_line->x1 = c1->x - c1->r * rot_x;
    tangent_line->y1 = c1->y - c1->r * rot_y;
  }
  return true;
}

static bool jsf__circle_tangent_circle(const DubinsCircle *c0, const DubinsCircle *c1,
                                       DubinsCircle *tangent_circle) {
  // Input circle radii must match.
  if (fabs(c0->r - c1->r) > DUBINS_EPSILON) {
    return false;
  }

  // Input circle directions must match.
  if (c0->direction != c1->direction) {
    return false;
  }

  // Compute the distance from c0 to c1.
  double c0c1_x = c1->x - c0->x;
  double c0c1_y = c1->y - c0->y;
  double D = jsf__norm(c0c1_x, c0c1_y);

  // No tangent circle (with same radius as c0 and c1)
  // is possible because c0 and c1 are too far apart.
  if (D > 4 * c0->r) {
    return false;
  }

  // Normalize the vector from the center of c0 to c1.
  double c0c1_hat_x = c0c1_x / D;
  double c0c1_hat_y = c0c1_y / D;

  // Compute the angle by which to rotate.
  double alpha = acos(D / 4 * (c0->r));
  alpha *= (c0->direction == DUBINS_RIGHT) ? -1.0 : 1.0;

  // Rotate c0c1_hat by alpha.
  double c = cos(alpha);
  double s = sin(alpha);
  double rot_x = c * c0c1_hat_x - s * c0c1_hat_y;
  double rot_y = s * c0c1_hat_x + c * c0c1_hat_y;

  // Compute the center of the tangent circle.
  tangent_circle->x = c0->x + 2 * c0->r * rot_x;
  tangent_circle->y = c0->y + 2 * c0->r * rot_y;

  // Fill in everything else.
  if (c0->direction == DUBINS_RIGHT) {
    tangent_circle->direction = DUBINS_LEFT;
  } else {
    tangent_circle->direction = DUBINS_RIGHT;
  }
  tangent_circle->r = c0->r;

  return true;
}

static double jsf__angle_between(double x0, double y0, double x1, double y1,
                                 DubinsDirection direction) {
  double n0 = jsf__norm(x0, y0);
  double n1 = jsf__norm(x1, y1);

  double arg = (x0 * x1 + y0 * y1) / (n0 * n1);

  // Clamp to (-1, +1) to avoid domain errors caused by rounding.
  arg = (arg > 1.0) ? 1.0 : arg;
  arg = (arg < -1.0) ? -1.0 : arg;

  double angle = acos(arg);

  double cross2d = x0 * y1 - x1 * y0;
  if (direction == DUBINS_LEFT && cross2d < 0) {
    angle = 2 * M_PI - angle;
  }
  if (direction == DUBINS_RIGHT && cross2d > 0) {
    angle = 2 * M_PI - angle;
  }
  // Anything close to 2pi counts as 0.
  // This avoids doing full circles when you don't need to.
  if (fabs(angle - 2 * M_PI) < DUBINS_EPSILON) {
    angle = 0.0;
  }
  return angle;
}

static void jsf__construct_csc(const DubinsState *s0, const DubinsState *s1,
                               const DubinsCircle *c0, const DubinsLine *line,
                               const DubinsCircle *c1, DubinsPath *path) {

  // Arc from s0, around c0, to line0
  {
    double v0x = s0->x - c0->x;
    double v0y = s0->y - c0->y;
    double v1x = line->x0 - c0->x;
    double v1y = line->y0 - c0->y;
    path->arclengths[0] =
        jsf__angle_between(v0x, v0y, v1x, v1y, c0->direction) * c0->r;
    path->curvatures[0] =
        (c0->direction == DUBINS_LEFT) ? 1.0 / c0->r : -1.0 / c0->r;
  }

  // Straight line from line0 to line1
  path->arclengths[1] = jsf__norm(line->x1 - line->x0, line->y1 - line->y0);
  path->curvatures[1] = 0.0;

  // Arc from line1, around c1, to s1
  {
    double v0x = line->x1 - c1->x;
    double v0y = line->y1 - c1->y;
    double v1x = s1->x - c1->x;
    double v1y = s1->y - c1->y;
    path->arclengths[2] =
        jsf__angle_between(v0x, v0y, v1x, v1y, c1->direction) * c1->r;
    path->curvatures[2] =
        (c1->direction == DUBINS_LEFT) ? 1.0 / c1->r : -1.0 / c1->r;
  }
}

static void jsf__construct_ccc(const DubinsState *s0, const DubinsState *s1,
                               const DubinsCircle *c0, const DubinsCircle *c1,
                               const DubinsCircle *c2, DubinsPath *path) {
  // Tangent point 0 is the average of the centers of c0 and c1.
  double tp0x = 0.5 * (c0->x + c1->x);
  double tp0y = 0.5 * (c0->y + c1->y);

  // Tangent point 1 is the average of the centers of c1 and c2.
  double tp1x = 0.5 * (c1->x + c2->x);
  double tp1y = 0.5 * (c1->y + c2->y);

  // Arc from s0, around c0, to tp0
  {
    double v0x = s0->x - c0->x;
    double v0y = s0->y - c0->y;
    double v1x = tp0x - c0->x;
    double v1y = tp0y - c0->y;
    path->arclengths[0] =
        jsf__angle_between(v0x, v0y, v1x, v1y, c0->direction) * c0->r;
    path->curvatures[0] =
        (c0->direction == DUBINS_LEFT) ? 1.0 / c0->r : -1.0 / c0->r;
  }

  // Arc from tp0, around c1, to tp1
  {
    double v0x = tp0x - c1->x;
    double v0y = tp0y - c1->y;
    double v1x = tp1x - c1->x;
    double v1y = tp1y - c1->y;
    path->arclengths[1] =
        jsf__angle_between(v0x, v0y, v1x, v1y, c1->direction) * c1->r;
    path->curvatures[1] =
        (c1->direction == DUBINS_LEFT) ? 1.0 / c1->r : -1.0 / c1->r;
  }

  // Arc from tp1, around c2, to s1
  {
    double v0x = tp1x - c2->x;
    double v0y = tp1y - c2->y;
    double v1x = s1->x - c2->x;
    double v1y = s1->y - c2->y;
    path->arclengths[2] =
        jsf__angle_between(v0x, v0y, v1x, v1y, c2->direction) * c2->r;
    path->curvatures[2] =
        (c2->direction == DUBINS_LEFT) ? 1.0 / c2->r : -1.0 / c2->r;
  }
}

bool jsf_dubins(const DubinsState *initial_state,
                const DubinsState *final_state, double radius,
                DubinsPath *path) {
  // Radius must be non-negative.
  radius = fabs(radius);

  // Construct the tangent circles to the initial and final states.
  DubinsCircle cl0, cr0, cl1, cr1;
  {
    jsf__state_tangent_circle(initial_state, radius, DUBINS_LEFT, &cl0);
    jsf__state_tangent_circle(initial_state, radius, DUBINS_RIGHT, &cr0);
    jsf__state_tangent_circle(final_state, radius, DUBINS_LEFT, &cl1);
    jsf__state_tangent_circle(final_state, radius, DUBINS_RIGHT, &cr1);
  }

  // Construct the tangent lines connecting the initial and final states.
  DubinsLine l2r, r2l, r2r, l2l;
  bool l2r_valid, r2l_valid, r2r_valid, l2l_valid;
  {
    l2r_valid = jsf__tangent_line(&cl0, &cr1, &l2r);
    r2l_valid = jsf__tangent_line(&cr0, &cl1, &r2l);
    r2r_valid = jsf__tangent_line(&cr0, &cr1, &r2r);
    l2l_valid = jsf__tangent_line(&cl0, &cl1, &l2l);
  }

  // Construct the tangent circles connecting the initial and final states.
  DubinsCircle crc, clc;
  bool crc_valid, clc_valid;
  {
    crc_valid = jsf__circle_tangent_circle(&cl0, &cl1, &crc);
    clc_valid = jsf__circle_tangent_circle(&cr0, &cr1, &clc);
  }

  DubinsPath rsr, lsl, lsr, rsl, lrl, rlr;

  double min_len = INFINITY;
  DubinsPath *shortest_path = NULL;

  // Construct the RSR path.
  if (r2r_valid) {
    jsf__construct_csc(initial_state, final_state, &cr0, &r2r, &cr1, &rsr);
    double len = jsf_dubins_path_length(&rsr);
    rsr.type = DUBINS_RSR;

    min_len = len;
    shortest_path = &rsr;
  }

  // Construct the LSL path.
  if (l2l_valid) {
    jsf__construct_csc(initial_state, final_state, &cl0, &l2l, &cl1, &lsl);
    double len = jsf_dubins_path_length(&lsl);
    lsl.type = DUBINS_LSL;

    if (len < min_len - DUBINS_EPSILON) {
      min_len = len;
      shortest_path = &lsl;
    }
  }

  // Construct the LSR path.
  if (l2r_valid) {
    jsf__construct_csc(initial_state, final_state, &cl0, &l2r, &cr1, &lsr);
    double len = jsf_dubins_path_length(&lsr);
    lsr.type = DUBINS_LSR;

    if (len < min_len - DUBINS_EPSILON) {
      min_len = len;
      shortest_path = &lsr;
    }
  }

  // Construct the RSL path.
  if (r2l_valid) {
    jsf__construct_csc(initial_state, final_state, &cr0, &r2l, &cl1, &rsl);
    double len = jsf_dubins_path_length(&rsl);
    rsl.type = DUBINS_RSL;

    if (len < min_len - DUBINS_EPSILON) {
      min_len = len;
      shortest_path = &rsl;
    }
  }

  // Construct the LRL path.
  if (crc_valid) {
    jsf__construct_ccc(initial_state, final_state, &cl0, &crc, &cl1, &lrl);
    double len = jsf_dubins_path_length(&lrl);
    lrl.type = DUBINS_LRL;

    if (len < min_len - DUBINS_EPSILON) {
      min_len = len;
      shortest_path = &lrl;
    }
  }

  // Construct the RLR path.
  if (clc_valid) {
    jsf__construct_ccc(initial_state, final_state, &cr0, &clc, &cr1, &rlr);
    double len = jsf_dubins_path_length(&rlr);
    rlr.type = DUBINS_RLR;

    if (len < min_len - DUBINS_EPSILON) {
      min_len = len;
      shortest_path = &rlr;
    }
  }

  path->arclengths[0] = shortest_path->arclengths[0];
  path->arclengths[1] = shortest_path->arclengths[1];
  path->arclengths[2] = shortest_path->arclengths[2];

  path->curvatures[0] = shortest_path->curvatures[0];
  path->curvatures[1] = shortest_path->curvatures[1];
  path->curvatures[2] = shortest_path->curvatures[2];
  path->type = shortest_path->type;

  return true;
}

double jsf_dubins_path_length(const DubinsPath *path) {
  return path->arclengths[0] + path->arclengths[1] + path->arclengths[2];
}

DubinsState jsf_dubins_path_sample(const DubinsState *initial_state,
                                   const DubinsPath *path,
                                   const double arclength) {

  const int num_arcs = sizeof(path->arclengths) / sizeof(path->arclengths[0]);

  int arc_idx = 0;
  double arc_prefix_length = 0.0;
  while (arc_idx < num_arcs &&
         arc_prefix_length + path->arclengths[arc_idx] < arclength) {
    arc_prefix_length += path->arclengths[arc_idx];
    arc_idx++;
  }

  // Find the state at the start of the arc containing the sample.
  DubinsState state = *initial_state;
  for (int i = 0; i < num_arcs && i <= arc_idx; ++i) {
    double k = path->curvatures[i];
    double s = path->arclengths[i];

    // This arc contains the sample.
    if (i == arc_idx) {
      s = arclength - arc_prefix_length;
    }

    s = fmax(0.0, s);

    // Arc is a straight line.
    if (k == 0.0) {
      state.x += s * cos(state.h);
      state.y += s * sin(state.h);
    } else {
      // Compute the radius of this arc.
      double r = 1.0 / k;

      // Find the tangent circle to the current state in the direction of this
      // arc.
      DubinsCircle circle;
      {
        DubinsDirection dir = (k > 0.0) ? DUBINS_LEFT : DUBINS_RIGHT;
        jsf__state_tangent_circle(&state, fabs(r), dir, &circle);
      }

      double dh = s * k;
      state.x = circle.x + r * cos(state.h - M_PI / 2.0 + dh);
      state.y = circle.y + r * sin(state.h - M_PI / 2.0 + dh);
      state.h += dh;
    }
  }
  return state;
}

#endif // JSF_DUBINS_IMPLEMENTATION

#ifdef __cplusplus
}
#endif

/*
------------------------------------------------------------------------------
MIT License
Copyright (c) 2021 Jordan Ford
------------------------------------------------------------------------------
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
*/
