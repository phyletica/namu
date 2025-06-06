#pragma once

#include <memory>
#include <iostream>
#include "node.hpp"

namespace phycoeval {

    class TreeManip;
    class TreeIO;
    //class Likelihood;
    //class Updater;

    class Tree {
        friend class TreeManip;
        friend class TreeIO;
        //friend class Likelihood;
        //friend class Updater;

        public:

                                        Tree();
                                        ~Tree();

            unsigned                    num_leaves() const;
            unsigned                    num_internals() const;
            unsigned                    num_nodes() const;

        private:

            void                        clear();

            Node *                      _root;
            unsigned                    _nleaves;
            unsigned                    _ninternals;
            Node::PtrVector             _preorder;
            Node::PtrVector             _levelorder;
            Node::Vector                _nodes;

        public:

            typedef std::shared_ptr< Tree > SharedPtr;
    };

    inline Tree::Tree() {
        std::cout << "Constructing a Tree" << std::endl;
        this->clear();
    }

    inline Tree::~Tree() {
        std::cout << "Destroying a Tree" << std::endl;
    }

    inline void Tree::clear() {
        this->_root = 0;
        this->_nodes.clear();
        this->_preorder.clear();
        this->_levelorder.clear();
    }

    inline unsigned Tree::num_leaves() const {
        return _nleaves;
    }

    inline unsigned Tree::num_internals() const {
        return _ninternals;
    }

    inline unsigned Tree::num_nodes() const {
        return (unsigned)this->_nodes.size();
    }
}
