#ifndef CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H

#define CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H



#include <CGAL/Constrained_Delaunay_triangulation_2.h>



template <class Gt, class Tds = CGAL::Default, class Itag = CGAL::Default>

class Custom_Constrained_Delaunay_triangulation_2

    : public CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag> {

public:

    using Base = CGAL::Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>;

    using typename Base::Face_handle;

    using typename Base::Point;

    using typename Base::Vertex_handle;

    using typename Base::Locate_type;



    // Constructors

    Custom_Constrained_Delaunay_triangulation_2(const Gt& gt = Gt())

        : Base(gt) {}



    Custom_Constrained_Delaunay_triangulation_2(typename Base::List_constraints& lc, const Gt& gt = Gt())

        : Base(lc, gt) {}



    template <class InputIterator>

    Custom_Constrained_Delaunay_triangulation_2(InputIterator it, InputIterator last, const Gt& gt = Gt())

        : Base(it, last, gt) {}



    // New insert method without flips

    Vertex_handle insert_no_flip(const Point& a, Face_handle start = Face_handle()) {

        // Call Ctr::insert without flip_around

        Vertex_handle va = this->Base::Ctr::insert(a, start); // Directly call Ctr::insert from the base

        return va;

    }



    // Another insert method with known location

    Vertex_handle insert_no_flip(const Point& a, Locate_type lt, Face_handle loc, int li) {

        Vertex_handle va = this->Base::Ctr::insert(a, lt, loc, li); // Directly call Ctr::insert from the base

        return va;

    }

        void remove_no_flip(Vertex_handle v) {

        this->Base::Ctr::remove(v); // Directly call Ctr::remove from the base

        

    }

    void insert_constraint_no_flip(const Point& a, Point& b) {

        this->Base::Ctr::insert_constraint(a,b); // Directly call Ctr::insert_constraint from the base

        

    }

    template <class OutputIterator>
    void remove_constrained_edge(Face_handle loc, int li,OutputIterator out){
        this->Base::Ctr::remove_constrained_edge(loc, li,out);  //Directly call Ctr::remove_constrained_edge from the base
    }

};



#endif // CGAL_CUSTOM_CONSTRAINED_DELAUNAY_TRIANGULATION_2_H