#pragma once

#include <string>
#include <vector>
#include <iostream>
#include "split.hpp"

namespace namu {
 
    class Tree;
    class TreeManip;
    class TreeIO;
    //class Likelihood;
    //class Updater;

    class Node {
        friend class Tree;
        friend class TreeManip;
        friend class TreeIO;
        friend class NodeTester;
        //friend class Likelihood;
        //friend class Updater;

        public:
                                    Node();
                                    ~Node();

            Node *                  get_parent() const      {return this->_parent;}
            Node *                  get_left_child() const  {return this->_left_child;}
            Node *                  get_right_sib() const   {return this->_right_sib;}
            Node *                  get_rightmost_child() const; 
            int                     get_number() const      {return this->_number;}
            std::string             get_name() const        {return this->_name;}
            Split                   get_spit() const        {return this->_split;}
            unsigned                get_number_of_children() const;

            double                  get_height() const      {return this->_height;}
            void                    set_height(double height);

            bool                    has_parent() const;
            bool                    has_child() const;
            bool                    is_root() const;
            bool                    is_leaf() const;
            bool                    is_polytomy() const;
            double                  get_edge_length() const;

            static const double _smallest_edge_length;

            typedef std::vector<Node>       Vector;
            typedef std::vector<Node *>     PtrVector;

        private:
 
            void                clear();
            Node *              _left_child;
            Node *              _right_sib;
            Node *              _parent;
            int                 _number;
            std::string         _name;
            double              _height;
            Split               _split;
    };

    inline Node::Node() {
        // std::cout << "Creating Node object" << std::endl;
        this->clear();
    }

    inline Node::~Node() {
        // std::cout << "Destroying Node object" << std::endl;
    }

    inline void Node::clear() {
        _left_child = 0;
        _right_sib = 0;
        _parent = 0;
        _number = -1;
        _name = "";
        _height = 0.0;
    }

    inline Node * Node::get_rightmost_child() const {
        if (this->_left_child == 0) {
            return this->_left_child;
        }
        Node * nd = this->get_left_child();
        Node * sib = nd->get_right_sib();
        while(sib) {
            nd = sib;
            sib = nd->get_right_sib();
        }
        return nd;
    }

    inline unsigned Node::get_number_of_children() const {
        unsigned n_children = 0;
        Node * nd = this->get_left_child();
        if (nd) {
            ++n_children;
            Node * sib = nd->get_right_sib();
            while(sib) {
                nd = sib;
                sib = nd->get_right_sib();
                ++n_children;
            }
        }
        return n_children;
    }

    inline void Node::set_height(double height) {
        this->_height = height;
    }

    inline bool Node::has_parent() const {
        return this->_parent == 0 ? false : true;
    }

    inline bool Node::has_child() const {
        return this->_left_child == 0 ? false : true;
    }

    inline bool Node::is_root() const {
        return this->_parent == 0 ? true : false;
    }

    inline bool Node::is_leaf() const {
        return this->_left_child == 0 ? true : false;
    }

    inline bool Node::is_polytomy() const {
        return this->get_number_of_children() > 2 ? true : false;
    }

    inline double Node::get_edge_length() const {
        if (this->has_parent()) {
            return this->_parent->get_height() - this->get_height();
        }
        return 0.0;
    }
}
