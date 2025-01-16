#include "include/Custom_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Convex_hull_traits_adapter_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/property_map.h>
#include <boost/json.hpp>
#include <ctime>
#include <fstream>
#include <iostream>
#include <queue>
#include <utility>
#include <vector>
#include <chrono>
#include <random>
#define CYCLE_MAX 50

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_triangulation_plus_2<Custom_Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::Exact_intersections_tag>> CDT;
typedef CGAL::Polygon_2<K> Polygon;
typedef K::Point_2 Point;
typedef K::Segment_2 Edge;
typedef K::Vector_2 Vector;
typedef CDT::Face_handle Face_handle;

namespace json = boost::json;

std::vector<std::string> steiner_points_x; // Global vector to store x-coordinates of Steiner points
std::vector<std::string> steiner_points_y; // Global vector to store y-coordinates of Steiner points
Polygon boundary_polygon;                  // Polygon for the region boundary
std::vector<K::FT> steiner_points_x_2;
std::vector<K::FT> steiner_points_y_2;
bool randomization = false;

// Fuction that help us to write the resuls on output file (convert a FT number to a string)
std::string FT_to_string(const K::FT &ft) {
  std::ostringstream oss;
  oss << ft;
  return oss.str();
}

std::string return_rational(const K::FT &coord) {
  const auto exact_coord = CGAL::exact(coord);
  return exact_coord.get_num().get_str() + "/" +
         exact_coord.get_den().get_str();
}

// Function to check if a specific triangle is obtuse with there edge length
int is_obtuse(const Point &p1, const Point &p2, const Point &p3) {
  // Calculate squared edge lengths
  K::FT a = CGAL::squared_distance(p2, p3); // Squared distance p2 and p3
  K::FT b = CGAL::squared_distance(p1, p3); // Squared distanc p1 and p3
  K::FT c = CGAL::squared_distance(p1, p2); // quared distance p1 and p2

  // if square lenght of the side opposite that angle is greater than the sum of
  // the squares of the other two sides.
  if (b + c < a) // Case obtuse angle is in opposite side of c
    return 1;
  else if (a + c < b) // Case obtuse angle is in opposite side of b
    return 2;
  else if (a + b < c) // Case obtuse angle is in opposite side of a
    return 3;
  else
    return 0; // If none true then the trianle is not obtuse
}

// Fuction that returns the count of obtuse faces
int return_obtuse(CDT &cdt, Polygon &pol) {
  int count = 0;
  for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
       ++face) {
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();
    if (pol.bounded_side(CGAL::centroid(p1, p2, p3)) !=
        CGAL::ON_UNBOUNDED_SIDE) {

      if ((is_obtuse(face->vertex(0)->point(), face->vertex(1)->point(),
                     face->vertex(2)->point())))
        count++;
    }
  }
  return count;
}

// Function to flip each possible flipable obtuse edge in a single pass
void flip_edges(CDT &cdt, Polygon boundary) {
  // Iterate all  edges
  for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
       ++face) {
    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();
    if (boundary.bounded_side(CGAL::centroid(p1, p2, p3)) !=
        CGAL::ON_UNBOUNDED_SIDE) {
      if (is_obtuse(p1, p2, p3) != 0) {
        for (int i = 0; i < 3; i++) {
          Point newp1 = face->vertex(i)->point();
          Point newp2 = face->vertex((i + 1) % 3)->point();
          Point newp3 = face->vertex((i + 2) % 3)->point();
          if (face->is_constrained(i) == false) {

            if (cdt.is_infinite(face->neighbor(i)) == false) {
              int count_bef = 1;
              int count_aft = 0;

              Face_handle neighbor = face->neighbor(i);
              int index_n = neighbor->index(face);
              Point p4 = neighbor->vertex(index_n)->point();
              Polygon pol;
              pol.push_back(newp1);
              pol.push_back(newp2);
              pol.push_back(p4);
              pol.push_back(newp3);
              if (boundary.bounded_side(CGAL::centroid(p1, p2, p3)) !=
                  CGAL::ON_UNBOUNDED_SIDE) {

                if (pol.is_convex() == true) {
                  if (is_obtuse(newp2, newp3, p4) != 0) {
                    count_bef = count_bef + 1;
                  }
                  if (is_obtuse(newp1, newp2, p4) != 0) {
                    count_aft = count_aft + 1;
                  }
                  if (is_obtuse(newp1, newp3, p4) != 0) {
                    count_aft = count_aft + 1;
                  }
                  if (count_aft < count_bef) {
                    Face_handle dupface;
                    dupface = face;
                    cdt.flip(dupface, i);
                    flip_edges(cdt, pol);
                    return;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

// Project point p onto the edge of p1 and p2
Point project_edge(const Point &p, const Point &p1, const Point &p2) {
  // Create vectors from the points
  K::Vector_2 v1 = p2 - p1;
  K::Vector_2 v2 = p - p1;

  K::FT dot_product = v1 * v2;         // Compute the dot product of v2 with v1
  K::FT length2 = v1.squared_length(); // Squared length of v1

  // Projection factor t
  K::FT t = dot_product / length2;

  // Compute the projection point using the parametric equation of the line
  return Point(p1.x() + t * v1.x(), p1.y() + t * v1.y());
}

// Fuction that returns the number of the best method for reducing obtuse
int min(int r1, int r2, int r3, int r4, int r5, int r6) {
  int min = r1;
  int min_num = 1;
  if (r2 < min) {
    min = r2;
    min_num = 2;
  }
  if (r3 < min) {
    min = r3;
    min_num = 3;
  }
  if (r4 < min) {
    min = r4;
    min_num = 4;
  }
  if (r5 < min) {
    min = r5;
    min_num = 5;
  }
  if (r6 < min) {
    min = r6;
    min_num = 6;
  }
  return min_num;
}

int min_local(int r1, int r2, int r3, int r4) {
  int min = r1;
  int min_num = 1;
  if (r2 < min) {
    min = r2;
    min_num = 2;
  }
  if (r3 < min) {
    min = r3;
    min_num = 3;
  }
  if (r4 < min) {
    min = r4;
    min_num = 4;
  }
  // if (r5 < min) {
  //   min = r5;
  //   min_num = 5;
  // }
  return min_num;
}

// Function to compute the midpoint (bisect)
Point midpoint_edge(const Point &p1, const Point &p2) {
  // Get the x coordinates
  K::FT x1 = p1.x();
  K::FT x2 = p2.x();

  // Get the y coordinates
  K::FT y1 = p1.y();
  K::FT y2 = p2.y();

  // Calculate the x of the midpoint
  K::FT midpoint_x = (x1 + x2) / 2.0;

  // Calculate the y of the midpoint
  K::FT midpoint_y = (y1 + y2) / 2.0;

  return Point(midpoint_x, midpoint_y);
}

// Simple fuction for putting steiner point in the first obtuse face
void insert_steiner_points_combined(CDT &cdt, int num_of_obtuse_before, Polygon &pol) {
  for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
       ++face) {
    CDT copy1, copy2, copy3, copy4, copy5, copy6;
    copy1 = cdt;
    copy2 = cdt;
    copy3 = cdt;
    copy4 = cdt;
    copy5 = cdt;

    copy6 = cdt;
    int testing = return_obtuse(copy1, pol);

    Point p1 = face->vertex(0)->point();
    Point p2 = face->vertex(1)->point();
    Point p3 = face->vertex(2)->point();

    // Find the obtuse edge
    int obtuse_edge = is_obtuse(p1, p2, p3);

    if ((pol.bounded_side(CGAL::centroid(p1, p2, p3)) !=
         CGAL::ON_UNBOUNDED_SIDE) &&
        (obtuse_edge != 0)) {
      Point steiner_point;
      Point midpoint;
      if (obtuse_edge == 1) {
        steiner_point = project_edge(p1, p2, p3); // Bisect edge (p2, p3)
        midpoint = midpoint_edge(p2, p3);
      } else if (obtuse_edge == 2) {
        steiner_point = project_edge(p2, p1, p3); // Bisect edge (p1, p3)
        midpoint = midpoint_edge(p1, p3);
      } else if (obtuse_edge == 3) {
        steiner_point = project_edge(p3, p1, p2); // Bisect edge (p1, p2)
        midpoint = midpoint_edge(p1, p2);
      }

      Face_handle neighbor = face->neighbor(obtuse_edge - 1);

      copy1.insert(CGAL::centroid(p1, p2, p3));
      copy2.insert_no_flip(steiner_point);
      if (pol.bounded_side(CGAL::circumcenter(p1, p2, p3)) !=
          CGAL::ON_UNBOUNDED_SIDE) {
        copy3.insert(CGAL::circumcenter(p1, p2, p3));
      }
      copy4.insert_no_flip(midpoint);
      copy5.insert(midpoint);
      copy6.insert(steiner_point);

      int ret1 = return_obtuse(copy1, pol);
      int ret2 = return_obtuse(copy2, pol);
      int ret3 = return_obtuse(copy3, pol);
      int ret4 = return_obtuse(copy4, pol);
      int ret5 = return_obtuse(copy5, pol);
      int ret6 = return_obtuse(copy6, pol);
      int min_ret = min(ret1, ret2, ret3, ret4, ret5, ret6);

      if (min_ret == 1) {
        if (ret1 < num_of_obtuse_before) {
          cdt.insert(CGAL::centroid(p1, p2, p3));

          // Save the centroid coordinates

          steiner_points_x_2.push_back((CGAL::centroid(p1, p2, p3).x()));
          steiner_points_y_2.push_back((CGAL::centroid(p1, p2, p3).y()));

          insert_steiner_points_combined(cdt, ret1, pol);
          return;
        }
      } else if (min_ret == 2) {
        if (ret2 < num_of_obtuse_before) {
          cdt.insert_no_flip(steiner_point);

          // Save the steiner point coordinates

          steiner_points_x_2.push_back(steiner_point.x());
          steiner_points_y_2.push_back(steiner_point.y());

          insert_steiner_points_combined(cdt, ret2, pol);
          return;
        }
      } else if (min_ret == 3) {
        if (ret3 < num_of_obtuse_before) {
          cdt.insert(CGAL::circumcenter(p1, p2, p3));

          // Save the circumcenter coordinates

          steiner_points_x_2.push_back((CGAL::circumcenter(p1, p2, p3).x()));
          steiner_points_y_2.push_back(CGAL::circumcenter(p1, p2, p3).y());

          insert_steiner_points_combined(cdt, ret3, pol);
          return;
        }
      } else if (min_ret == 4) {
        if (ret4 < num_of_obtuse_before) {
          cdt.insert_no_flip(midpoint);

          // Save the midpoint coordinates

          steiner_points_x_2.push_back(midpoint.x());
          steiner_points_y_2.push_back(midpoint.y());

          insert_steiner_points_combined(cdt, ret4, pol);
          return;
        }
      } else if (min_ret == 5) {
        if (ret5 < num_of_obtuse_before) {
          cdt.insert(midpoint);

          // Save the midpoint coordinates

          steiner_points_x_2.push_back(midpoint.x());
          steiner_points_y_2.push_back(midpoint.y());

          insert_steiner_points_combined(cdt, ret5, pol);
          return;
        }
      } else if (min_ret == 6) {
        if (ret6 < num_of_obtuse_before) {
          cdt.insert(steiner_point);

          // Save the midpoint coordinates

          steiner_points_x_2.push_back(steiner_point.x());
          steiner_points_y_2.push_back(steiner_point.y());

          insert_steiner_points_combined(cdt, ret6, pol);
          return;
        }
      }
    }
  }
  return;
}

// Check the face has an obtuse neighbor
bool find_obtuse_neighbor(CDT &cdt, Face_handle face, Face_handle &obtuse_neighbor, int &obtuse_neighbor_idx) {
  for (int i = 0; i < 3; ++i) {
    Face_handle neighbor = face->neighbor(i);
    if (!cdt.is_infinite(neighbor)) {
      Point p1 = neighbor->vertex(0)->point();
      Point p2 = neighbor->vertex(1)->point();
      Point p3 = neighbor->vertex(2)->point();

      // Check if the neighboring face is obtuse
      if (is_obtuse(p1, p2, p3) != 0 &&
          (boundary_polygon.bounded_side(CGAL::centroid(p1, p2, p3)) ==
           CGAL::ON_BOUNDED_SIDE)) {
        obtuse_neighbor = neighbor;
        obtuse_neighbor_idx = i;
        return true;
      }
    }
  }
  return false;
}

// Insert Steiner points between all neighboring obtuse faces
void insert_steiner_point_between_obtuse_neighbors(CDT &cdt, Polygon &pol) {
  bool found_pair = true;
  while (found_pair) {
    found_pair = false;

    // Iterate over all faces to find pairs of neighboring obtuse faces
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         ++face) {
      Point p1 = face->vertex(0)->point();
      Point p2 = face->vertex(1)->point();
      Point p3 = face->vertex(2)->point();

      // Check if the face is obtuse
      if (is_obtuse(p1, p2, p3) != 0 &&
          (boundary_polygon.bounded_side(CGAL::centroid(p1, p2, p3)) ==
           CGAL::ON_BOUNDED_SIDE)) {
        Face_handle obtuse_neighbor;
        int obtuse_neighbor_idx;

        // Check if there is an obtuse neighboring face
        if (find_obtuse_neighbor(cdt, face, obtuse_neighbor,
                                 obtuse_neighbor_idx)) {
          // Get the points of the neighboring obtuse face
          Point np1 = obtuse_neighbor->vertex(0)->point();
          Point np2 = obtuse_neighbor->vertex(1)->point();
          Point np3 = obtuse_neighbor->vertex(2)->point();

          // Compute the centroids of both obtuse faces
          Point centroid1 = CGAL::centroid(p1, p2, p3);
          Point centroid2 = CGAL::centroid(np1, np2, np3);

          // Compute the midpoint between the centroids
          Point steiner_point = midpoint_edge(centroid1, centroid2);

          // Insert the Steiner point in a CDT copy
          CDT temp_cdt = cdt;
          temp_cdt.insert(steiner_point);

          // Insert only if the obtuse count is reduced
          if (return_obtuse(temp_cdt, pol) < return_obtuse(cdt, pol) &&
              (boundary_polygon.bounded_side(steiner_point) ==
               CGAL::ON_BOUNDED_SIDE)) {
            cdt.insert(steiner_point);
            found_pair = true;
            break;
          }
        }
      }
    }
    // If there are no more obtose neigbor pairs then break
    if (!found_pair) {
      break;
    }
  }
}

bool point_not_in_vector(std::vector<Point> &points, Point p) {
  for (int i = 0; i < points.size(); i++) {
    if (points[i] == p) {
      return false;
    }
  }
  return true;
}

bool check_if_vertex_constrained(CDT &cdt, CDT::Edge_circulator &c1) {
  CDT::Edge startedge = *c1;
  CDT::Edge edge;
  do {
    edge = *c1;
    if (cdt.is_infinite(edge)) {
      // std::cout<<"Infinite edge"<<std::endl;
      c1++;
      continue;
    }
    if (cdt.is_constrained(edge)) {
      // std::cout<<"Constrained edge"<<std::endl;
      return true;

    } else {
      // std::cout<<"Not Constrained edge"<<std::endl;
    }
    c1++;

  } while (*c1 != startedge);
  return false;
}

void insert_circumcenter(CDT &cdt, Point p1, Point p2, Point p3, Polygon &bound) {
  Face_handle face;
  for (auto face2 = cdt.finite_faces_begin(); face2 != cdt.finite_faces_end();
       face2++) {
    Point facep1 = face2->vertex(0)->point();
    Point facep2 = face2->vertex(1)->point();
    Point facep3 = face2->vertex(2)->point();
    if ((facep1 == p1 && facep2 == p2 && facep3 == p3)) {
      face = face2;

      break;
    }
  }
  Point circcenter = CGAL::circumcenter(p1, p2, p3);
  for (int i = 0; i < 3; i++) {
    Face_handle neighbor = face->neighbor(i);
    if (cdt.is_infinite(neighbor) == true) {
      continue;
    }
    int constrained = 0;
    // for (int j=0; j<3; j++){
    //     if (neighbor->is_constrained(j)==true){
    //         constrained=1;
    //     }
    // }
    // if (constrained==1){
    //     continue;                                            //elegxos an o
    //     geitonas periexei constraints
    // mexri tora den dhmiourgei provlimata an den elegxoume otan kanoume insert
    // nea kai meta remove constraint an dhmiourgountai provlimata dokimaste me
    // uncomment
    // }
    int index_n;
    index_n = neighbor->index(face);
    Point np1 = neighbor->vertex((index_n + 2) % 3)->point();
    Point np2 = neighbor->vertex(index_n)->point();
    Point np3 = neighbor->vertex((index_n + 1) % 3)->point();
    Polygon npolygon;
    npolygon.push_back(np1); // DHMIOURGIA POLIGONOY GEITONON
    npolygon.push_back(np2);
    npolygon.push_back(np3);
    std::vector<Point> points;
    std::vector<CDT::Vertex_handle> vertices;

    if ((npolygon.bounded_side(circcenter) != CGAL::ON_UNBOUNDED_SIDE) &&
        (bound.bounded_side(CGAL::centroid(np1, np2, np3)) !=
         CGAL::ON_UNBOUNDED_SIDE)) {
      vertices.clear();
      points.clear();
      points.push_back(face->vertex((i + 2) % 3)->point());
      points.push_back(face->vertex(i)->point());
      points.push_back(face->vertex((i + 1) % 3)->point());
      points.push_back(np2);
      vertices.push_back(face->vertex((i + 2) % 3));
      vertices.push_back(face->vertex(i));
      vertices.push_back(face->vertex((i + 1) % 3));
      vertices.push_back(neighbor->vertex(index_n));

      for (int k = 0; k < vertices.size(); k++) {

        CDT::Vertex_handle vr = vertices[k];
        CDT::Edge_circulator circ = vr->incident_edges();
        if (check_if_vertex_constrained(cdt, circ) ==
            false) { // ELEGXOS VERTICES TOY POLIGOUNOU POU ANHKOYN SE
                     // CONSTRAINT. AN OXI MPOREI NA AFAIRETHEI. ALLIOS OXI
          cdt.remove_no_flip(vr);
        }
      }

      std::vector<CDT::Constraint_id> constraints;

      for (int k = 0; k < points.size() - 1; k++) {
        CDT::Constraint_id constraint_id =
            cdt.insert_constraint(points[k], points[k + 1]);
        constraints.push_back(
            constraint_id); // EISAGOSI SIMEION POLIGONOY ME CONSTRAINT
      }

      CDT::Constraint_id cid =
          cdt.insert_constraint(points[0], points[points.size() - 1]);

      constraints.push_back(cid);
      cdt.insert_no_flip(circcenter);

      int obtuse_count;
      int new_obtuse_count;
      do {
        obtuse_count = return_obtuse(cdt, bound);
        flip_edges(cdt, bound);
        new_obtuse_count = return_obtuse(
            cdt, bound); // EISAGOGI SIMEIOY KAI ELEGXOS ME FLIP AN EXEI MPOREI
                         // NA AFAIRETHEI TO EDGE POY DIAXORISEI TOYS GEITONES
                         // ME FLIP TOY. AN TO EDGE POY DIAXORISEI TA TRIGONA
                         // DEN MPOREI NA AFAIRETHEI EPEIDI KAI TA DIO VERTICES
                         // ANHKOYN SE KAPOIO CONSTRAINT

      } while (obtuse_count > new_obtuse_count);

      for (int k = 0; k < constraints.size(); k++) {
        cdt.remove_constraint(constraints[k]); // AFAIRESI CONSTRAINTS
      }

      return;
    }
  }
}

Point unify_triangles(CDT &cdt, Point p1, Point p2, Point p3, Polygon &bound) {
  std::vector<Point> points;
  std::vector<CDT::Vertex_handle> vertices;
  int obtuse_neighbors = 0;
  Point mean;
  // EISAGOSI SIMEIOY ME ENOSI TRIGONON. PAROMOIA ME THN CIRCUMCENTER ALLA
  // FTIAXNEI POLIGONO ME POLLAPLA TRIGONA.
  Face_handle face;

  for (auto face2 = cdt.finite_faces_begin(); face2 != cdt.finite_faces_end();
       face2++) {
    Point facep1 = face2->vertex(0)->point();
    Point facep2 = face2->vertex(1)->point();
    Point facep3 = face2->vertex(2)->point();
    if ((facep1 == p1 && facep2 == p2 && facep3 == p3)) {
      face = face2;

      break;
    }
  }

  for (int j = 0; j < 3; j++) {
    if (face->is_constrained(j) == true) {

      return mean;
    }
  }

  for (int i = 0; i < 3; i++) {

    int constrained = 0;
    Face_handle neighbor = face->neighbor(i);
    if (cdt.is_infinite(neighbor) == true) {
      if (point_not_in_vector(points, face->vertex((i + 1) % 3)->point()) ==
          true) {
        Point p = face->vertex((i + 1) % 3)->point();
        points.push_back(p);
        vertices.push_back(face->vertex((i + 1) % 3));
      }

      if (point_not_in_vector(points, face->vertex((i + 2) % 3)->point()) ==
          true) {

        points.push_back(face->vertex((i + 2) % 3)->point());
        vertices.push_back(face->vertex((i + 2) % 3));
      }
      continue;
    }
    for (int j = 0; j < 3; j++) {
      if (neighbor->is_constrained(j) == true) {

        constrained = 1;
      }
    }

    int index_n;
    index_n = neighbor->index(face);
    Point p1 = neighbor->vertex((index_n + 2) % 3)->point();
    Point p2 = neighbor->vertex(index_n)->point();
    Point p3 = neighbor->vertex((index_n + 1) % 3)->point();

    if ((constrained == 1) || (is_obtuse(p1, p2, p3) == 0) ||
        (bound.bounded_side(CGAL::centroid(p1, p2, p3)) ==
         CGAL::ON_UNBOUNDED_SIDE)) {
      if (point_not_in_vector(points, face->vertex((i + 1) % 3)->point()) ==
          true) {
        points.push_back(face->vertex((i + 1) % 3)->point());
        vertices.push_back(face->vertex((i + 1) % 3));
      }
      if (point_not_in_vector(points, face->vertex((i + 2) % 3)->point()) ==
          true) {

        points.push_back(face->vertex((i + 2) % 3)->point());
        vertices.push_back(face->vertex((i + 2) % 3));
      }
      continue;
    }
    obtuse_neighbors++;
    if (point_not_in_vector(points, p1) == true) {
      points.push_back(p1);
      vertices.push_back(neighbor->vertex((index_n + 2) % 3));
    }
    if (point_not_in_vector(points, p2) == true) {
      points.push_back(p2);
      vertices.push_back(neighbor->vertex((index_n) % 3));
    }
    if (point_not_in_vector(points, p3) == true) {

      points.push_back(p3);
      vertices.push_back(neighbor->vertex((index_n + 1) % 3));
    }
  }
  Polygon built(points.begin(), points.end());

  if (obtuse_neighbors == 0 ||
      CGAL::is_convex_2(points.begin(), points.end()) == false) {
    return mean;
  }
  K::FT sum_x = 0;
  K::FT sum_y = 0;
  K::FT size = points.size();

  for (int i = 0; i < points.size(); i++) {

    sum_x = sum_x + points[i].x();
    sum_y = sum_y + points[i].y();

    CDT::Vertex_handle vr = vertices[i];
    CDT::Edge_circulator circ = vr->incident_edges();
    if (check_if_vertex_constrained(cdt, circ) == false) {
      cdt.remove_no_flip(vr);
    }
  }

  mean = Point(sum_x / size, sum_y / size);

  std::vector<CDT::Constraint_id> constraints;

  for (int i = 0; i < points.size() - 1; i++) {
    CDT::Constraint_id constraint_id =
        cdt.insert_constraint(points[i], points[i + 1]);
    constraints.push_back(constraint_id);
  }

  CDT::Constraint_id cid =
      cdt.insert_constraint(points[0], points[points.size() - 1]);

  constraints.push_back(cid);
  cdt.insert_no_flip(mean);

  int obtuse_count;
  int new_obtuse_count;
  do {
    obtuse_count = return_obtuse(cdt, bound);
    flip_edges(cdt, bound);
    new_obtuse_count = return_obtuse(cdt, bound);

  } while (obtuse_count > new_obtuse_count);

  for (int i = 0; i < constraints.size(); i++) {
    cdt.remove_constraint(constraints[i]);
  }

  return mean;
}

void randomization_insert(CDT &cdt, Polygon &pol) {
    std::vector<Face_handle> obtuse_faces;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();
        if (is_obtuse(p1, p2, p3)) {
            obtuse_faces.push_back(face);
        }
    }

    if (obtuse_faces.empty()) {
        printf("No obtuse triangles found.\n");
        return;
    }

    //Calculate Gaussian standard deviation based on average triangle size
    double total_area = 0.0;
    int triangle_count = 0;

    for (const auto &face : obtuse_faces) {
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        //Calculate area of each triangle
        double area = std::abs(CGAL::to_double(
            CGAL::area(p1, p2, p3)
        ));
        total_area += area;
        triangle_count++;
    }

    //Compute Gaussian standard deviation
    double gaussian_stddev = std::sqrt(total_area / triangle_count) / 2.0;

    //Gaussian random number generator
    std::random_device rd;
    boost::random::mt19937 gen(rd());
    boost::random::normal_distribution<> d(0, gaussian_stddev);


    int initial_obtuse = return_obtuse(cdt, pol); //Count initial obtuse triangles
    int obtuse_reduced = 0;

    for (const auto &face : obtuse_faces) {
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        //Calculate the centroid of the triangle
        Point centroid = CGAL::centroid(p1, p2, p3);

        //Generate a random point near the centroid using Gaussian distribution
        Point random_point(
            centroid.x() + d(gen),
            centroid.y() + d(gen)
        );

        //Check if the random point is inside the polygon
        if (pol.bounded_side(random_point) == CGAL::ON_UNBOUNDED_SIDE) {
            continue;
        }

        // Create a temporary CDT to evaluate the effect of inserting the random point
        CDT temp_cdt = cdt;
        temp_cdt.insert(random_point);

        int new_obtuse = return_obtuse(temp_cdt, pol);

        //If the random point reduces obtuse triangles, insert it into the original CDT
        if (new_obtuse < initial_obtuse) {
            bool randomization = true;
            cdt.insert(random_point);
            initial_obtuse = new_obtuse;
            obtuse_reduced++;
            steiner_points_x_2.push_back(random_point.x());
            steiner_points_y_2.push_back(random_point.y());
        }
    }

    printf("Randomization completed. Obtuse triangles reduced by %d.\n", obtuse_reduced);
}



void local_method(CDT &cdt, Polygon &pol, int L) {
  int count = 0;
  auto start_time = std::chrono::steady_clock::now();

  while (count < L) {
    int num_of_obtuse_before = return_obtuse(cdt, pol);

    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end();
         ++face) {
      CDT copy1, copy2, copy3, copy4, copy5, copy6;
      copy1 = cdt;
      copy2 = cdt;
      copy3 = cdt;
      copy4 = cdt;
      copy5 = cdt;
      // printf("im stuck\n");

      Point p1 = face->vertex(0)->point();
      Point p2 = face->vertex(1)->point();
      Point p3 = face->vertex(2)->point();

      // Find the obtuse edge

      int obtuse_edge = is_obtuse(p1, p2, p3);

      if ((pol.bounded_side(CGAL::centroid(p1, p2, p3)) !=
           CGAL::ON_UNBOUNDED_SIDE) &&
          (obtuse_edge != 0)) {
        Point steiner_point;
        Point midpoint;
        if (obtuse_edge == 1) {
          steiner_point = project_edge(p1, p2, p3); // Bisect edge (p2, p3)
          midpoint = midpoint_edge(p2, p3);
        } else if (obtuse_edge == 2) {
          steiner_point = project_edge(p2, p1, p3); // Bisect edge (p1, p3)
          midpoint = midpoint_edge(p1, p3);
        } else if (obtuse_edge == 3) {
          steiner_point = project_edge(p3, p1, p2); // Bisect edge (p1, p2)
          midpoint = midpoint_edge(p1, p2);
        }

        Face_handle neighbor = face->neighbor(obtuse_edge - 1);

        copy1.insert(CGAL::centroid(p1, p2, p3));
        copy2.insert_no_flip(steiner_point);
        if (pol.bounded_side(CGAL::circumcenter(p1, p2, p3)) !=
            CGAL::ON_UNBOUNDED_SIDE) {
          // insert_circumcenter(copy3, p1, p2, p3, pol);
          copy3.insert(CGAL::circumcenter(p1,p2,p3));
        }
        copy4.insert_no_flip(midpoint);
        // Point mean = unify_triangles(copy5, p1, p2, p3, pol);

        int ret1 = return_obtuse(copy1, pol);
        int ret2 = return_obtuse(copy2, pol);
        int ret3 = return_obtuse(copy3, pol);
        int ret4 = return_obtuse(copy4, pol);
        // int ret5 = return_obtuse(copy5, pol);

        int min_ret = min_local(ret1, ret2, ret3, ret4);

        if (min_ret == 1) {
          if (ret1 < num_of_obtuse_before) {
            // cdt.insert(CGAL::centroid(p1,p2,p3));
            
            cdt = copy1;

            // Save the centroid coordinates

            steiner_points_x_2.push_back((CGAL::centroid(p1, p2, p3).x()));
            steiner_points_y_2.push_back((CGAL::centroid(p1, p2, p3).y()));

            break;
          }
        } else if (min_ret == 2) {
          if (ret2 < num_of_obtuse_before) {
            // cdt.insert_no_flip(steiner_point);
            
            cdt = copy2;

            // Save the steiner point coordinates

            steiner_points_x_2.push_back(steiner_point.x());
            steiner_points_y_2.push_back(steiner_point.y());

            break;
          }
        } else if (min_ret == 3) {
          if (ret3 < num_of_obtuse_before) {
            // cdt.insert(CGAL::circumcenter(p1,p2,p3));

            
            cdt = copy3;

            // Save the circumcenter coordinates

            steiner_points_x_2.push_back((CGAL::circumcenter(p1, p2, p3).x()));
            steiner_points_y_2.push_back(CGAL::circumcenter(p1, p2, p3).y());

            break;
          }
        } else if (min_ret == 4) {
          if (ret4 < num_of_obtuse_before) {

            
            cdt = copy4;

            // Save the midpoint coordinates

            steiner_points_x_2.push_back(midpoint.x());
            steiner_points_y_2.push_back(midpoint.y());

            break;
          }
        // } else if (min_ret == 5) {

        //   if (ret5 < num_of_obtuse_before) {

        //     cdt.clear();
        //     cdt = copy5;

        //     // Save the midpoint coordinates

        //     steiner_points_x_2.push_back(mean.x());
        //     steiner_points_y_2.push_back(mean.y());

        //     return;
        //   }
        }
      }
    }


    count++;
    //Check elapsed time
      auto current_time = std::chrono::steady_clock::now();
      auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

      //If 10-15 seconds have passed with no improvement calling randomization
      if (elapsed_seconds >= 60) {
          printf("Deadlock detected after %ld seconds. Calling randomization...\n", elapsed_seconds);
          randomization_insert(cdt, pol);
          start_time = std::chrono::steady_clock::now(); 
      }
  }

  return;
}

void sa_method(CDT &cdt, Polygon &pol, double alpha, double beta, int L) {
  int initial_obtuse = return_obtuse(cdt, pol);
  int initial_steiner_points = steiner_points_x_2.size();
  double energy = alpha * initial_obtuse + beta * initial_steiner_points;
  srand(time(0));

  double T = 1.0;
  bool improvement;

  // Start timing
  auto start_time = std::chrono::steady_clock::now();
  while (T > 0.0) {
    improvement = false;

    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
      Point p1 = face->vertex(0)->point();
      Point p2 = face->vertex(1)->point();
      Point p3 = face->vertex(2)->point();

      int obtuse_edge = is_obtuse(p1, p2, p3);
      if (obtuse_edge == 0)
        continue;

      Point project_e;
      Point midpoint;
      if (obtuse_edge == 1) {
        project_e = project_edge(p1, p2, p3); // Bisect edge (p2, p3)
        midpoint = midpoint_edge(p2, p3);
      } else if (obtuse_edge == 2) {
        project_e = project_edge(p2, p1, p3); // Bisect edge (p1, p3)
        midpoint = midpoint_edge(p1, p3);
      } else if (obtuse_edge == 3) {
        project_e = project_edge(p3, p1, p2); // Bisect edge (p1, p2)
        midpoint = midpoint_edge(p1, p2);
      }

      std::vector<Point>
          steiner_points; // Vector of points for the 5 possible steiner points

      Face_handle neighbor;
      int obtuse_neighbor_idx;
      Point neigbr;
      // if (find_obtuse_neighbor(cdt, face, neighbor, obtuse_neighbor_idx)) {
      //     Point np1 = neighbor->vertex(0)->point();
      //     Point np2 = neighbor->vertex(1)->point();
      //     Point np3 = neighbor->vertex(2)->point();

      //     if (is_obtuse(np1, np2, np3) != 0) {
      //         neigbr = midpoint_edge(CGAL::centroid(p1, p2, p3),
      //         CGAL::centroid(np1, np2, np3));
      //         steiner_points.push_back(neigbr);
      //     }
      // }

      steiner_points.push_back(project_e);
      steiner_points.push_back(midpoint);
      steiner_points.push_back(CGAL::circumcenter(p1, p2, p3));
      steiner_points.push_back(CGAL::centroid(p1, p2, p3));

      Point selected_point = steiner_points[rand() % steiner_points.size()];

      if (pol.bounded_side(selected_point) == CGAL::ON_UNBOUNDED_SIDE)
        continue;

      CDT temp_cdt = cdt;
      if (selected_point == CGAL::centroid(p1, p2, p3) ||
          selected_point == CGAL::circumcenter(p1, p2, p3)) {
        temp_cdt.insert(selected_point);
      } else {
        temp_cdt.insert_no_flip(selected_point);
      }

      int new_obtuse = return_obtuse(temp_cdt, pol);
      int new_steiner_points = steiner_points_x_2.size() + 1;

      double new_energy = alpha * new_obtuse + beta * new_steiner_points;

      double d_energy = new_energy - energy;

      if (d_energy < 0 ||
          (exp((-1) * d_energy / T) >=
           ((double)rand() / RAND_MAX))) { // Propability between 0.0 and 1.0
        if (selected_point == CGAL::centroid(p1, p2, p3) ||
            selected_point == CGAL::circumcenter(p1, p2, p3)) {
          cdt.insert(selected_point);

        } else {
          cdt.insert_no_flip(selected_point);
        }

        energy = new_energy;
        improvement = true;
        steiner_points_x_2.push_back((selected_point.x()));
        steiner_points_y_2.push_back((selected_point.y()));

        break;
      }
    }

    T -= 1.0 / L;

      //Check elapsed time
      auto current_time = std::chrono::steady_clock::now();
      auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();

      //If 10-15 seconds have passed with no improvement calling randomization
      if (elapsed_seconds >= 60 && !improvement) {
          printf("Deadlock detected after %ld seconds. Calling randomization...\n", elapsed_seconds);
          randomization_insert(cdt, pol);
          start_time = std::chrono::steady_clock::now(); 
      }
  }
}

double findcircumradius(Face_handle &face) {
  // Find the circumradius of the face
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();

  Point circenter = CGAL::circumcenter(p1, p2, p3);
  K::FT sqdistance = CGAL::squared_distance(circenter, p1);
  double sqdistanced = CGAL::to_double(sqdistance);
  double radius = double(std::sqrt(sqdistanced));

  return radius;
}

double findheight(Face_handle face) {
  // sinartiseis gia circumradius kai ipsos gia na vrethei to ρ
  Point p1 = face->vertex(0)->point();
  Point p2 = face->vertex(1)->point();
  Point p3 = face->vertex(2)->point();
  int obt_angle = is_obtuse(p1, p2, p3);
  Point projection;
  K::FT sqdistance;
  if (obt_angle == 1) {
    projection = project_edge(p1, p2, p3);
    sqdistance = CGAL::squared_distance(p1, projection);

  } else if (obt_angle == 2) {
    projection = project_edge(p2, p1, p3);
    sqdistance = CGAL::squared_distance(p2, projection);
  } else {
    projection = project_edge(p3, p1, p2);
    sqdistance = CGAL::squared_distance(p3, projection);
  }

  double sqdistanced = CGAL::to_double(sqdistance);
  double distance = std::sqrt(sqdistanced);

  return distance;
}

bool equal_faces_points(Point f1p1, Point f1p2, Point f1p3, Point f2p1, Point f2p2, Point f2p3) {
  int count = 0;

  if (f1p1 == f2p1 || f1p1 == f2p2 || f1p1 == f2p3) {
    count++;
  }
  if (f1p2 == f2p1 || f1p2 == f2p2 || f1p2 == f2p3) {
    count++;
  }
  if (f1p2 == f2p1 || f1p2 == f2p2 || f1p2 == f2p3) {
    count++;
  }
  if (count == 3) {
    return true;
  }
  return false;
}

void track_changed_triangles_and_conflicts(CDT &originalcdt, CDT &newCdt, std::vector<std::vector<Point>> &conflicts, std::vector<int> &array_of_obtuse_am, int num_of_obtuse) {
  std::vector<Point> changed_faces;
  for (auto face = originalcdt.finite_faces_begin();
       face != originalcdt.finite_faces_end(); face++) {
    Face_handle facedup;

    facedup = face;
    int check = 0;
    for (auto face2 = newCdt.finite_faces_begin();
         face2 != newCdt.finite_faces_end(); face2++) {
      Face_handle facedup2;
      facedup2 = face2;

      if (equal_faces_points(
              facedup->vertex(0)->point(), facedup->vertex(1)->point(),
              facedup->vertex(2)->point(), facedup2->vertex(0)->point(),
              facedup2->vertex(1)->point(),
              facedup2->vertex(2)->point()) == true) {
        check = 1;
        break;
      }
    }
    if (check == 1) {
      continue;
    } else {
      changed_faces.push_back(face->vertex(0)->point());
      changed_faces.push_back(face->vertex(1)->point());
      changed_faces.push_back(face->vertex(2)->point());
    }

    // Face_handle new_face=newCdt.locate(face->vertex(0)->point());
    // if (new_face!=newCdt.faces_end() &&!check_faces_points(face,new_face)){
    //   changed_faces.push_back(face);
    // }
  }
  bool conflict = false;
  int conflict_index = -1;
  if (changed_faces.size() > 0) {
    for (int i = 0; i < conflicts.size(); i++) {
      std::vector<Point> temp_vector;
      temp_vector = conflicts[i];
      for (int j = 0; j < temp_vector.size(); j = j + 3) {
        for (int k = 0; i < changed_faces.size(); k = k + 3) {
          if (equal_faces_points(temp_vector[i], temp_vector[j + 1],
                                 temp_vector[j + 2], changed_faces[k],
                                 changed_faces[k + 1],
                                 changed_faces[k + 2]) == true) {
            conflict = true;
            conflict_index = i;
            break;
          }
          if (conflict == true) {
            break;
          }
        }
        if (conflict == true) {
          break;
        }
      }
      if (conflict == true) {
        if (num_of_obtuse < array_of_obtuse_am[conflict_index]) {
          array_of_obtuse_am.erase(array_of_obtuse_am.begin() + conflict_index);
          conflicts.erase(conflicts.begin() + conflict_index);
          // phertracker.erase()
          //...
          conflicts.push_back(changed_faces);
          array_of_obtuse_am.push_back(num_of_obtuse);
          // phertracker.push_back()
          //
        }

      } else {
        conflicts.push_back(changed_faces);
        array_of_obtuse_am.push_back(num_of_obtuse);
        // phertracker.push_back()
        //
      }
    }
  }
}




void track_changed_triangles_and_conflicts(std::vector<std::vector<Point>> &conflicts,std::vector<int> &array_of_obtuse_am,std::vector<Point> &points_inserted,
    Point point_coord, std::vector<double> &pher_per_point, double pheromone_val,  int num_of_obtuse,
    std::vector<int> &method_per_point,int method_used, std::vector<Point> changed_faces) {
  
  // for (auto face = originalcdt.finite_faces_begin();
  //      face != originalcdt.finite_faces_end(); face++) {
  //   Face_handle facedup;

  //   facedup = face;
  //   int check = 0;
  //   for (auto face2 = newCdt.finite_faces_begin();
  //        face2 != newCdt.finite_faces_end(); face2++) {
  //     Face_handle facedup2;
  //     facedup2 = face2;

  //     if (equal_faces_points(
  //             facedup->vertex(0)->point(), facedup->vertex(1)->point(),
  //             facedup->vertex(2)->point(), facedup2->vertex(0)->point(),
  //             facedup2->vertex(1)->point(),
  //             facedup2->vertex(2)->point()) == true) {
  //       check = 1;
  //       break;
  //     }
  //   }
  //   if (check == 1) {
  //     continue;
  //   } else {
  //     changed_faces.push_back(face->vertex(0)->point());
  //     changed_faces.push_back(face->vertex(1)->point());
  //     changed_faces.push_back(face->vertex(2)->point());
  //   }

  //   // Face_handle new_face=newCdt.locate(face->vertex(0)->point());
  //   // if (new_face!=newCdt.faces_end() &&!check_faces_points(face,new_face)){
  //   //   changed_faces.push_back(face);
  //   // }
  // }
  bool conflict = false;
  int conflict_index = -1;
  if (changed_faces.size() > 0) {
    for (int i = 0; i < conflicts.size(); i++) {
      std::vector<Point> temp_vector=conflicts[i];
      
      for (int j = 0; j < temp_vector.size(); j = j + 3) {
        for (int k = 0; i < changed_faces.size(); k = k + 3) {
          if (equal_faces_points(temp_vector[i], temp_vector[j + 1],
                                 temp_vector[j + 2], changed_faces[k],
                                 changed_faces[k + 1],
                                 changed_faces[k + 2]) == true) {
            conflict = true;
            conflict_index = i;
            break;
          }
          if (conflict == true) {
            break;
          }
        }
        if (conflict == true) {
          break;
        }
      }
      if (conflict == true) {
        if (num_of_obtuse < array_of_obtuse_am[conflict_index]) {
          array_of_obtuse_am.erase(array_of_obtuse_am.begin() + conflict_index);
          conflicts.erase(conflicts.begin() + conflict_index);
          pher_per_point.erase(pher_per_point.begin() + conflict_index);
          method_per_point.erase(method_per_point.begin() + conflict_index);
          points_inserted.erase(points_inserted.begin() + conflict_index);
          // phertracker.erase()
          //...
          conflicts.push_back(changed_faces);
          array_of_obtuse_am.push_back(num_of_obtuse);
        }
          

      } 
      else {
        conflicts.push_back(changed_faces);
        array_of_obtuse_am.push_back(num_of_obtuse);
        pher_per_point.push_back(pheromone_val);
        method_per_point.push_back(method_used);
        points_inserted.push_back(point_coord);
        
      }
    }
  }
  
}

void insert_circumcenter_for_ant(CDT &cdt, Point p1, Point p2, Point p3,
                         Polygon &bound,std::vector<Point> &faces_changed) {
  Face_handle face;
  for (auto face2 = cdt.finite_faces_begin(); face2 != cdt.finite_faces_end();
       face2++) {
    Point facep1 = face2->vertex(0)->point();
    Point facep2 = face2->vertex(1)->point();
    Point facep3 = face2->vertex(2)->point();
    if ((facep1 == p1 && facep2 == p2 && facep3 == p3)) {
      face = face2;

      break;
    }
  }
  Point circcenter = CGAL::circumcenter(p1, p2, p3);
  for (int i = 0; i < 3; i++) {
    Face_handle neighbor = face->neighbor(i);
    if (cdt.is_infinite(neighbor) == true) {
      continue;
    }
    int constrained = 0;
    // for (int j=0; j<3; j++){
    //     if (neighbor->is_constrained(j)==true){
    //         constrained=1;
    //     }
    // }
    // if (constrained==1){
    //     continue;                                            //elegxos an o
    //     geitonas periexei constraints
    // mexri tora den dhmiourgei provlimata an den elegxoume otan kanoume insert
    // nea kai meta remove constraint an dhmiourgountai provlimata dokimaste me
    // uncomment
    // }
    int index_n;
    index_n = neighbor->index(face);
    Point np1 = neighbor->vertex((index_n + 2) % 3)->point();
    Point np2 = neighbor->vertex(index_n)->point();
    Point np3 = neighbor->vertex((index_n + 1) % 3)->point();
    Polygon npolygon;
    npolygon.push_back(np1); // DHMIOURGIA POLIGONOY GEITONON
    npolygon.push_back(np2);
    npolygon.push_back(np3);
    std::vector<Point> points;
    std::vector<CDT::Vertex_handle> vertices;

    if ((npolygon.bounded_side(circcenter) != CGAL::ON_UNBOUNDED_SIDE) &&
        (bound.bounded_side(CGAL::centroid(np1, np2, np3)) !=
         CGAL::ON_UNBOUNDED_SIDE)) {
      vertices.clear();
      points.clear();
      points.push_back(face->vertex((i + 2) % 3)->point());
      points.push_back(face->vertex(i)->point());
      points.push_back(face->vertex((i + 1) % 3)->point());
      points.push_back(np2);
      vertices.push_back(face->vertex((i + 2) % 3));
      vertices.push_back(face->vertex(i));
      vertices.push_back(face->vertex((i + 1) % 3));
      vertices.push_back(neighbor->vertex(index_n));

      for (int k = 0; k < vertices.size(); k++) {

        CDT::Vertex_handle vr = vertices[k];
        CDT::Edge_circulator circ = vr->incident_edges();
        if (check_if_vertex_constrained(cdt, circ) ==
            false) { // ELEGXOS VERTICES TOY POLIGOUNOU POU ANHKOYN SE
                     // CONSTRAINT. AN OXI MPOREI NA AFAIRETHEI. ALLIOS OXI
          cdt.remove_no_flip(vr);
        }
      }

      std::vector<CDT::Constraint_id> constraints;

      for (int k = 0; k < points.size() - 1; k++) {
        CDT::Constraint_id constraint_id =
            cdt.insert_constraint(points[k], points[k + 1]);
        constraints.push_back(
            constraint_id); // EISAGOSI SIMEION POLIGONOY ME CONSTRAINT
      }

      CDT::Constraint_id cid =
          cdt.insert_constraint(points[0], points[points.size() - 1]);

      constraints.push_back(cid);
      cdt.insert_no_flip(circcenter);

      int obtuse_count;
      int new_obtuse_count;
      do {
        obtuse_count = return_obtuse(cdt, bound);
        flip_edges(cdt, bound);
        new_obtuse_count = return_obtuse(
            cdt, bound); // EISAGOGI SIMEIOY KAI ELEGXOS ME FLIP AN EXEI MPOREI
                         // NA AFAIRETHEI TO EDGE POY DIAXORISEI TOYS GEITONES
                         // ME FLIP TOY. AN TO EDGE POY DIAXORISEI TA TRIGONA
                         // DEN MPOREI NA AFAIRETHEI EPEIDI KAI TA DIO VERTICES
                         // ANHKOYN SE KAPOIO CONSTRAINT

      } while (obtuse_count > new_obtuse_count);

      for (int k = 0; k < constraints.size(); k++) {
        cdt.remove_constraint(constraints[k]); // AFAIRESI CONSTRAINTS
      }
      faces_changed.push_back(face->vertex(0)->point()); 
      faces_changed.push_back(face->vertex(1)->point());
      faces_changed.push_back(face->vertex(2)->point());
      faces_changed.push_back(np1); 
      faces_changed.push_back(np2);
      faces_changed.push_back(np3);

      return;
    }
  }
}
    
void antwork(CDT &cdt,Polygon pol,std::vector<double> phertracker,std::vector<Point> points_inserted,
std::vector<std::vector<Point>> &faces_changed_per_ant,std::vector<int> &obtuse_am_per_point,std::vector<int> &method_per_point,int xi, int psi,double alpha,double beta ){
  std::vector<Point> faces_changed;
    srand(time(0));
    int amountofobtuse=return_obtuse(cdt,pol);
    int random_triangle=rand()%amountofobtuse;
    int count=0;
    for (auto face=cdt.finite_faces_begin();
    face!=cdt.finite_faces_end();face++){
        Point p1=face->vertex(0)->point();
        Point p2=face->vertex(1)->point();
        Point p3=face->vertex(2)->point();
        int obtuse_angle=is_obtuse(p1,p2,p3);
        if (obtuse_angle !=0 && (pol.bounded_side(CGAL::centroid(p1,p2,p3))
        !=CGAL::ON_UNBOUNDED_SIDE)){

        if (count==random_triangle){
        Face_handle curface=face;
        double facecr=findcircumradius(curface)/findheight(curface);

        double hmidpoint=3-2*facecr/facecr;
        if (hmidpoint<0){
            hmidpoint=0;
        }
        double hproj=(facecr-1)/facecr;
        if (hproj<0){
            hproj=0;
        }
        int hcirc=facecr/(2+facecr);
        double sumofmethods=pow(phertracker[1],xi)*pow(hmidpoint,psi) +
        pow(phertracker[0],xi)*pow(hproj,psi) + pow(phertracker[2],xi)+pow(hcirc,psi);
        double thmid=pow(phertracker[1],xi)*pow(hmidpoint,psi);
        double thproj=pow(phertracker[0],xi)*pow(hproj,psi);
        double thcirc=pow(phertracker[2],xi)+pow(hcirc,psi);
        double probmid=thmid/sumofmethods;
        double probproj=probmid+thproj/sumofmethods;
        double probcirc=probproj+thcirc/sumofmethods;
        double random_num=rand()/(RAND_MAX);
        Point midpoint;
        Point projection;
        cdt.clear();

        if (obtuse_angle == 1)
                {
                    projection = project_edge(p1, p2, p3); // Bisect edge
                    midpoint = midpoint_edge(p2, p3);
                }
                else if (obtuse_angle == 2)
                {
                    projection = project_edge(p2, p3, p1); // Bisect edge
                    midpoint = midpoint_edge(p3, p1);
                }
                else if (obtuse_angle == 3)
                {
                    projection = project_edge(p3, p1, p2); // Bisect edge
                    midpoint = midpoint_edge(p1, p2);
                }

        if (random_num<=probmid){
          faces_changed.push_back(face->vertex(0)->point());
          faces_changed.push_back(face->vertex(1)->point());
          faces_changed.push_back(face->vertex(2)->point());
          Face_handle neighbor=face->neighbor(obtuse_angle-1);
          if (cdt.is_infinite(neighbor) == false){
            faces_changed.push_back(neighbor->vertex(0)->point());
            faces_changed.push_back(neighbor->vertex(1)->point());
            faces_changed.push_back(neighbor->vertex(2)->point());

          }
          cdt.insert_no_flip (midpoint);
          int num_of_obtuse=return_obtuse(cdt,pol);
          double pher_val=1/(alpha*num_of_obtuse + beta* (steiner_points_x_2.size()+1));
          track_changed_triangles_and_conflicts(faces_changed_per_ant,obtuse_am_per_point,points_inserted,midpoint,phertracker,pher_val,num_of_obtuse,method_per_point,1,faces_changed);


        }
        else if(random_num<=probproj){
          faces_changed.push_back(face->vertex(0)->point());
          faces_changed.push_back(face->vertex(1)->point());
          faces_changed.push_back(face->vertex(2)->point());
          Face_handle neighbor=face->neighbor(obtuse_angle-1);
          if (cdt.is_infinite(neighbor) == false){
            faces_changed.push_back(neighbor->vertex(0)->point());
            faces_changed.push_back(neighbor->vertex(1)->point());
            faces_changed.push_back(neighbor->vertex(2)->point());
            int num_of_obtuse=return_obtuse(cdt,pol);
            double pher_val=1/(alpha*num_of_obtuse + beta* (steiner_points_x_2.size()+1));
            track_changed_triangles_and_conflicts(faces_changed_per_ant,obtuse_am_per_point,points_inserted,projection,phertracker,pher_val,num_of_obtuse,method_per_point,0,faces_changed);
          }
            cdt.insert_no_flip (projection);

        }
        else if (random_num<=probcirc){
          Point circec=CGAL::circumcenter(p1,p2,p3);
          std::vector<Point> changed_faces;
          insert_circumcenter_for_ant(cdt,p1,p2,p3,pol,changed_faces);
          
          if (changed_faces.size()>0){

          int num_of_obtuse=return_obtuse(cdt,pol);
          double pher_val=1/(alpha* num_of_obtuse+ beta* (steiner_points_x_2.size()+1));
          track_changed_triangles_and_conflicts(faces_changed_per_ant,obtuse_am_per_point,points_inserted,circec,phertracker,pher_val,num_of_obtuse,method_per_point,2,faces_changed);
          }
        }
        }
        count++;
        }
        int num_obtuse;
        

    }
}

void antmethod(CDT &cdt, Polygon &pol, double alpha, double beta, double xi,
               double psi, double lambda, int kappa, int L) {
  int cycleamount = 50;
  int bestobtcount = return_obtuse(cdt, pol);
  double evap=0.5;
  std::vector<double> phertracker;
  std::vector<int> point_method;
  std::vector<double> pher_per_ant;
  std::vector<std::vector<Point>> conflicts;
  std::vector<Point> points_inserted;
  std::vector<int> obtuse_per_ant;
  int num_of_obtuse=return_obtuse(cdt,pol);
  phertracker.push_back(0.5);                     //projection pher
  phertracker.push_back(0.5);                     //midpoint pher
  phertracker.push_back(0.5);                     //circumcenter pher
  CDT besttriang;
  for (int i = 0; i < CYCLE_MAX; i++) {
    pher_per_ant.clear();
    points_inserted.clear();
    conflicts.clear();
    point_method.clear();
    obtuse_per_ant.clear();
    CDT copy = cdt;

    for (int j = 0; j < 50; i++) {

      
      antwork(cdt,pol,pher_per_ant,points_inserted,conflicts,obtuse_per_ant,point_method,xi,psi,alpha,beta);
    }
    for (int k=0; k<phertracker.size(); k++){
      phertracker[k]=evap*phertracker[k];         //evaporation
    }
    for (int k=0; k<pher_per_ant.size(); k++){
      phertracker[k] = phertracker[k] + pher_per_ant[k];      //addition of pheromone to each method
    }
    for (int k=0; k<points_inserted.size();k++){
      if (point_method[k]==0 || point_method[k]==1){
        cdt.insert_no_flip(points_inserted[k]);

      }
      else{
        copy.insert(points_inserted[k]);
      }
      
      
    }
    int new_num_of_obtuse= return_obtuse(copy,pol);
    if (new_num_of_obtuse<num_of_obtuse){
      cdt=copy;
      num_of_obtuse=new_num_of_obtuse;
      for (int k=0; k<points_inserted.size(); k++){
        steiner_points_x_2.push_back(points_inserted[k].x());
        steiner_points_y_2.push_back(points_inserted[k].y());

      }

    }

  }
  // check energy
  // update pheromone
}


bool check_for_closed_constraints(json::array constraints, std::vector<int> region_bound) {
    std::queue<int> fifo_queue;
    std::queue<int> prev_queue;
    std::unordered_set<int> visited; 

    //Iterate through all constraints
    for (const auto &constraint : constraints) {
        int first_index = constraint.as_array()[0].as_int64();
        int second_index = constraint.as_array()[1].as_int64();
        bool is_visited = visited.find(first_index) != visited.end();
        int boundary_cross_count = 0;

        //If the point has not been visited, start BFS from it
        if (!is_visited) {
            fifo_queue.push(first_index);
            prev_queue.push(-1);

            //Check if the first index is on the region boundary
            if (std::find(region_bound.begin(), region_bound.end(), first_index) != region_bound.end()) {
                boundary_cross_count = 1;
            }
        }

        //BFS to explore constraints
        while (!fifo_queue.empty()) {
            int index = fifo_queue.front();
            fifo_queue.pop();
            int prev = prev_queue.front();
            prev_queue.pop();

            for (const auto &constraint2 : constraints) {
                int next_point = -1;
                bool found = false;

                //Determine the next point to explore
                if (constraint2.as_array()[0].as_int64() == index) {
                    next_point = constraint2.as_array()[1].as_int64();
                    found = true;
                } else if (constraint2.as_array()[1].as_int64() == index) {
                    next_point = constraint2.as_array()[0].as_int64();
                    found = true;
                }

                if (found && next_point != prev) {
                    //If the point is already visited, a closed constraint is found
                    if (visited.find(next_point) != visited.end()) {
                        return true;
                    }

                    //Check if the next point lies on the boundary
                    int edge_count_1 = std::count(region_bound.begin(), region_bound.end(), next_point);
                    int edge_count_2 = std::count(region_bound.begin(), region_bound.end(), index);

                    if (edge_count_1 == 1 && edge_count_2 == 0) {
                        boundary_cross_count++;
                    }

                    //If boundary cross count exceeds 1, a closed constraint is found
                    if (boundary_cross_count > 1) {
                        return true;
                    }

                    //Add the next point to the BFS queue
                    fifo_queue.push(next_point);
                    prev_queue.push(index);
                }
            }
        }
    }
    return false;
}



int main(int argc, char *argv[]) {
  // Check for correct arguments
  bool preselected_params = false;
  bool delaunay = true;
  boost::json::object parameters = {};
  std::string method, best_method; 

  if (argc < 5) {
    std::cout << "Error" << std::endl;
    return 1;
  }
  std::string input_file, output_file;

  // Parse arguments
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "-i") == 0 && i + 1 < argc) {
      input_file = argv[++i];
    } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
      output_file = argv[++i];
    }
    if(std::strcmp(argv[i], "–preselected_params") == 0){
      preselected_params = true;
    }
  }

  // JSON file
  std::ifstream input_stream(input_file);
  boost::json::value data = boost::json::parse(input_stream);
  input_stream.close();

  // Take all the data from JSON file//

  // Instance_uid
  std::string instance_uid = data.at("instance_uid").as_string().c_str();

  // Number of points
  int num_points = data.at("num_points").as_int64();

  // Points_x and points_y from JSON file
  std::vector<int> points_x;
  for (const auto &val : data.at("points_x").as_array())
    points_x.push_back(val.as_int64());

  std::vector<int> points_y;
  for (const auto &val : data.at("points_y").as_array())
    points_y.push_back(val.as_int64());

  // The region boundary
  const auto &reg = data.at("region_boundary").as_array();
  std::vector<int> region_boundary;

  for (const auto &index : reg)
    region_boundary.push_back(index.as_int64());

  // Create the region boundary polygon
  std::vector<Point> boundary;
  for (int p : region_boundary) {
    double x = points_x[p];
    double y = points_y[p];
    boundary.push_back(Point(x, y));
  }
  boundary_polygon = Polygon(boundary.begin(), boundary.end());

  // Number of constraints
  int num_constraints = data.at("num_constraints").as_int64();

  // Additional constraints from JSON file
  json::array constraints = data.at("additional_constraints").as_array();

  // Initialize the Constrained Delanay Triangulation (CDT)
  CDT cdt;

  // Insert the Points (x,y cordinates)
  std::vector<Point> points;
  for (size_t i = 0; i < points_x.size(); ++i) {
    points.push_back(Point(points_x[i], points_y[i]));
    cdt.insert(points.back());
  }

  // Insert constrained edges based on the provided indices
  for (const Point &p : points)
    cdt.insert(p);

  // Insert constraints
  for (const auto &constraint : constraints) {
    int first_index = constraint.as_array()[0].as_int64();
    int second_index = constraint.as_array()[1].as_int64();
    cdt.insert_constraint(points[first_index], points[second_index]);
  }

  // Check if the Triangulation is valid
  assert(cdt.is_valid());

  for (int i = region_boundary.size() - 1; i > 0; i--) {
    cdt.insert_constraint(
        Point(points_x[region_boundary[i]], points_y[region_boundary[i]]),
        Point(points_x[region_boundary[i - 1]],
              points_y[region_boundary[i - 1]]));
  }
  cdt.insert_constraint(
      Point(points_x[region_boundary[0]], points_y[region_boundary[0]]),
      Point(points_x[region_boundary[region_boundary.size() - 1]],
            points_y[region_boundary[region_boundary.size() - 1]]));
  int category = 5; // Default to Category 5
  bool res = check_for_closed_constraints(constraints, region_boundary);

  if (boundary_polygon.is_convex()) {
    if (num_constraints == 0) {
        category = 1;        //Convex boundary without constraints
    } else if (res == true) {
        category = 3;     //Convex boundary with closed constraints
    } else {
        category = 2;    //Convex boundary with open constraints
    }
  } else {
    //Check if all boundary segments are parallel to x or y axes
    bool paralleltox_or_y = true;

    for (int i = 1; i < region_boundary.size(); i++) {
        int pointx = points_x[region_boundary[i]] - points_x[region_boundary[i - 1]];
        int pointy = points_y[region_boundary[i]] - points_y[region_boundary[i - 1]];
        if (pointx != 0 && pointy != 0) {
            paralleltox_or_y = false;
            break;
        }
    }
    int regsize = region_boundary.size();
    int pointx = points_x[region_boundary[0]] - points_x[region_boundary[regsize - 1]];
    int pointy = points_y[region_boundary[0]] - points_y[region_boundary[regsize - 1]];

    if (pointx != 0 && pointy != 0) {
        paralleltox_or_y = false;
    }
    if (paralleltox_or_y == true && num_constraints==0) {
        category = 4;               //Non-convex boundary with axis-parallel segments
    } else {
        category = 5;               //Irregular, non-convex boundary
    }
  }

  if(!preselected_params) {
    // Delaunay, method and parameter variebles
    delaunay = data.at("delaunay").as_bool();
    method = data.at("method").as_string().c_str();
    parameters = data.at("parameters").as_object();
  }
  
  // After classifying the input and determining the category
  if (category == 1) {
    printf("Case A: Convex boundary without restrictions.\n");
    best_method = "local";
    printf("Choosen method: Local search.\n");

  } else if (category == 2) {
    printf("Case B: Convex boundary with open constraints.\n");
    best_method = "sa";
    printf("Choosen method: Simulated annealing.\n");

  } else if (category == 3) {
    printf("Case C: Convex boundary with closed constraints.\n");
    best_method = "sa";
    printf("Choosen method: Simulated annealing.\n");

  } else if (category == 4) {
    printf("Case D: Non-convex boundary with straight segments parallel to the axes.\n");
    best_method = "sa";
    printf("Choosen method: Simulated annealing.\n");

  } else if (category == 5) {
    printf("Case E: Non-convex boundary, irregular, not included in categories A-D.\n");
    best_method = "local";
    printf("Choosen method: Local search.\n");

  } else {
    printf("Unknown category. Choosed  Simulated annealing by default.\n");
    best_method = "sa";
  }

  int obtuse_count;

  if (delaunay == true) {
    flip_edges(cdt, boundary_polygon);
    //insert_steiner_point_between_obtuse_neighbors(cdt, boundary_polygon);
    insert_steiner_points_combined(cdt, return_obtuse(cdt, boundary_polygon), boundary_polygon);
  }

  // If delanay parameter is YES the we use the delanay_original cdt, if is NO
  // we use the cdt that is worked in first part of project Method
  if (best_method == "sa") {
    double alpha = 5.0;
    double beta = 0.8;
    int L = 200;
  
    if(preselected_params==true) {
      if(method == "sa" || method == "ant") {            //An exei dothei sto input.json ws method ontws to sa h to ant tote pernoume ta argouments tou input
         alpha = parameters.at("alpha").as_double();
         beta = parameters.at("beta").as_double();
         L = parameters.at("L").as_int64();
      }
    } 
    else{
      parameters.emplace("alpha", 5.0);      
      parameters.emplace("beta", 0.8);       
      parameters.emplace("L", 200);
    }
    
    printf("Alpha: %.1f\n",  alpha);
    printf("Beta: %.1f\n",   beta);
    printf("L: %d\n", L);
    sa_method(cdt, boundary_polygon, alpha, beta, L);
    obtuse_count = return_obtuse(cdt, boundary_polygon);
  } 
  else if (best_method == "ant") {
    double alpha = 4.0;
    double beta = 0.8;
    double xi = 1.0;
    double psi = 3.0;
    double lambda = 0.5;
    int kappa = 10;
    int L = 200;

    if(!preselected_params) {
      if(method == "ant") {
         alpha = parameters.at("alpha").as_double();
         beta = parameters.at("beta").as_double();
         xi = parameters.at("xi").as_double();
         psi = parameters.at("psi").as_double();
         lambda = parameters.at("lambda").as_double();
         kappa = parameters.at("kappa").as_int64();
         L = parameters.at("L").as_int64();
      }
    } 
    else{
      parameters.emplace("alpha", 4.0);      
      parameters.emplace("beta", 0.8);    
      parameters.emplace("xi", 1.0);      
      parameters.emplace("psi", 3.0);      
      parameters.emplace("lambda", 0.5);
      parameters.emplace("kappa", 10);
      parameters.emplace("L", 200);
    }
   
    std::cout << std::fixed << std::setprecision(1);
    printf("Alpha: %.1f\n", alpha);
    printf("Beta: %.1f\n", beta);
    printf("Xi: %.1f\n", xi);
    printf("Psi: %.1f\n", psi);
    printf("Lambda: %.1f\n", lambda);
    printf("kappa: %d\n", kappa);
    printf("L: %d\n", L);
    //ant_method(cdt, boundary_polygon, best_alpha, best_beta, best_xi, best_psi, best_lambda, best_kappa, best_L);
    obtuse_count = return_obtuse(cdt, boundary_polygon);

  } 
  else if (best_method == "local") {
    int L = 1000;

    if(!preselected_params){
       L = parameters.at("L").as_int64();
    }
    else{
      parameters.emplace("L", 1000);
    }
    
    printf("L: %d\n", L);
    local_method(cdt, boundary_polygon, L);
    obtuse_count = return_obtuse(cdt, boundary_polygon);
  }

  for (int i = 0; i < steiner_points_x_2.size(); i++) {
    steiner_points_x.push_back(return_rational(steiner_points_x_2[i]));
    steiner_points_y.push_back(return_rational(steiner_points_y_2[i]));
  }

  std::cout << "Number of obtuse: " << obtuse_count << std::endl;
  std::cout << "Number of steiner points: " << steiner_points_x.size() << std::endl;

  // Create JSON output
  json::object output;
  output["content_type"] = "CG_SHOP_2025_Solution";
  output["instance_uid"] = instance_uid;
  output["steiner_points_x"] = json::array(steiner_points_x.begin(), steiner_points_x.end());
  output["steiner_points_y"] = json::array(steiner_points_y.begin(), steiner_points_y.end());

  json::array edges;
  int count = 0;

  for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
    CDT::Face_handle face = edge->first;
    int idx = edge->second;
    count++;

    Point p1 = face->vertex((idx + 1) % 3)->point();
    Point p2 = face->vertex((idx + 2) % 3)->point();
    int v1, v2;
    Point midpoint((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2);
    if (boundary_polygon.bounded_side(midpoint) != CGAL::ON_UNBOUNDED_SIDE) {
      for (int i = 0; i < points.size(); ++i) {
        if (points[i] == p1) {
          v1 = i;
          break;
        }
        v1 = -1;
      }
      if (v1 == -1) {
        for (int i = 0; i < steiner_points_x_2.size(); ++i) {
          Point steiner(steiner_points_x_2[i], steiner_points_y_2[i]);
          if (steiner == p1) {
            v1 = num_points + i;
            break;
          }
          v1 = -1;
        }
      }
      for (int i = 0; i < points.size(); ++i) {
        if (points[i] == p2) {
          v2 = i;
          break;
        }
        v2 = -1;
      }
      if (v2 == -1) {
        for (int i = 0; i < steiner_points_x_2.size(); ++i) {
          Point steiner(steiner_points_x_2[i], steiner_points_y_2[i]);
          if (steiner == p2) {
            v2 = num_points + i;
            break;
          }
          v2 = -1;
        }
      }
      if (v1 != -1 && v2 != -1) {
        edges.push_back({v1, v2});
      }
    }
  }

  output["edges"] = edges;
  output["obtuse_count"] = obtuse_count;
  output["method"] = best_method;
  output["parameters"] = parameters;
  output["randomization"] = randomization;

  std::ofstream output_stream(output_file);
  output_stream.is_open();
  output_stream << json::serialize(output);
  output_stream.close();

  // Draw the triangulation using CGAL's draw fuction
  CGAL::draw(cdt);

  return 0;
}