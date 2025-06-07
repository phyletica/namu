#pragma once

#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#include <regex>
#include <limits>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/format.hpp>
#include <ncl/nxsmultiformat.h>
#include "namu/tree.hpp"
#include "namu/error.hpp"

namespace namu {

    class TreeIO;

    class TreeManip {

        friend class TreeIO;

        public:
                                            TreeManip();
                                            TreeManip(Tree::SharedPtr t);
                                            ~TreeManip();

            void                            set_tree(Tree::SharedPtr t);
            Tree::SharedPtr                 get_tree() const;

            double                          calc_tree_length() const;
            unsigned                        count_edges() const;
            unsigned                        count_nodes() const;
            void                            scale_all_internal_heights(double scaler);

            std::string                     make_newick(unsigned precision, bool use_names = false) const;

            void                            store_splits(std::set<Split> & split_set);

            static bool                     compare_node_vecs(
                                                const Node::PtrVector & a,
                                                const Node::PtrVector & b);

            void                            clear();

        private:

            void                            refresh_preorder();
            void                            refresh_levelorder();
            void                            refresh_div_event_order();
            Node *                          find_next_preorder(Node * nd) const;
            void                            renumber_internals();
            bool                            no_gaps_in_div_events() const;

            Tree::SharedPtr                 _tree;
            std::vector< Node::PtrVector >  _div_events;

        public:

            typedef std::shared_ptr< TreeManip >        SharedPtr;
            typedef std::vector< TreeManip::SharedPtr > SharedPtrVector;
    };

    inline TreeManip::TreeManip() {
        std::cout << "Constructing a TreeManip" << std::endl;
        this->clear();
    }

    inline TreeManip::TreeManip(Tree::SharedPtr t) {
        std::cout << "Constructing a TreeManip with a supplied tree" << std::endl;
        this->clear();
        this->set_tree(t);
    }

    inline TreeManip::~TreeManip() {
        std::cout << "Destroying a TreeManip" << std::endl;
    }

    inline void TreeManip::clear() {
        this->_tree.reset();
        this->_div_events.clear();
    }

    inline void TreeManip::set_tree(Tree::SharedPtr t) {
        assert(t);
        this->_tree = t;
    }

    inline Tree::SharedPtr TreeManip::get_tree() const {
        return this->_tree;
    }

    inline double TreeManip::calc_tree_length() const {
        double length = 0.0;
        for (auto node : this->_tree->_preorder) {
            length += node->get_edge_length();
        }
        return length;
    }

    inline unsigned TreeManip::count_edges() const {
        unsigned num_nodes = this->count_nodes();
        if (num_nodes < 2) {
            return 0;
        }
        return num_nodes - 1;
    }

    inline unsigned TreeManip::count_nodes() const {
        return (unsigned)this->_tree->_preorder.size();
    }

    inline void TreeManip::scale_all_internal_heights(double scaler) {
        for (auto node : this->_tree->_preorder) {
            if (node->has_child()) {
                node->set_height(scaler * node->get_height());
            }
        }
    }

    inline bool TreeManip::compare_node_vecs(
            const Node::PtrVector & a,
            const Node::PtrVector & b) {
        if (a.empty() && (! b.empty())) {
            return false;
        }
        if (b.empty() && (! a.empty())) {
            return true;
        }
        if (b.empty() && a.empty()) {
            return false;
        }
        return a.at(0)->_height < b.at(0)->_height;
    }

    inline void TreeManip::refresh_div_event_order() {
        std::sort(
                this->_div_events.begin(),
                this->_div_events.end(),
                TreeManip::compare_node_vecs);
    }

    inline bool TreeManip::no_gaps_in_div_events() const {
        bool found_end = false;
        for (auto node_vec : this->_div_events) {
            if (found_end && node_vec.size() > 0) {
                return false;
            }
            if (node_vec.size() == 0) {
                found_end = true;
            }
        }
        return true;
    }

    inline Node * TreeManip::find_next_preorder(Node * nd) const {
        assert(nd);
        Node * next = 0;
        if (!nd->_left_child && !nd->_right_sib) {
            // node has no children or siblings, so next preorder is the right
            // sibling of the first ancestral node that has a right sibling
            Node * anc = nd->_parent;
            while (anc && !anc->_right_sib) {
                anc = anc->_parent;
            }
            if (anc) {
                // We found an ancestor with a right sibling
                next = anc->_right_sib;
            }
            else {
                // node is last preorder node in the tree
                next = 0;
            }
        }
        else if (nd->_right_sib && !nd->_left_child) {
            // node has no children (it's a tip), but does have a right sibling
            next = nd->_right_sib;
        }
        else if (nd->_left_child && !nd->_right_sib) {
            // node has children but no right siblings
            next = nd->_left_child;
        }
        else {
            // node has both children and right siblings
            next = nd->_left_child;
        }
        return next;
    }

    inline void TreeManip::refresh_preorder() {
        this->_tree->_preorder.clear();
        this->_tree->_preorder.reserve(this->_tree->_nodes.size());

        if (! this->_tree->_root) {
            return;
        }

        Node * nd = this->_tree->_root;
        this->_tree->_preorder.push_back(nd);

        while (true) {
            nd = this->find_next_preorder(nd);
            if (nd)
                this->_tree->_preorder.push_back(nd);
            else
                break;
        }
    }

    inline void TreeManip::refresh_levelorder() {
        if (! this->_tree->_root) {
            return;
        }

        // q is the buffer queue
        std::queue<Node *> q;

        // _tree->_levelorder is the stack vector
        this->_tree->_levelorder.clear();
        this->_tree->_levelorder.reserve(this->_tree->_nodes.size());

        Node * nd = this->_tree->_root;

        q.push(nd);

        while (! q.empty()) {
            // pop nd off front of queue
            nd = q.front(); q.pop();

            // and push it onto the stack
            this->_tree->_levelorder.push_back(nd);

            // add all children of nd to back of the queue
            Node * child = nd->_left_child;
            if (child) {
                q.push(child);
                child = child->_right_sib;
                while(child) {
                    q.push(child);
                    child = child->_right_sib;
                }
            }
        }
    }

    inline void TreeManip::renumber_internals() {
        assert(this->_tree->_preorder.size() > 0);

        // Renumber internals in postorder sequence
        unsigned curr_internal = this->_tree->_nleaves;
        unsigned internal_count = 0;
        for (auto nd : boost::adaptors::reverse(this->_tree->_preorder)) {
            if (nd->_left_child) {
                // nd is an internal node
                nd->_number = curr_internal++;
                ++internal_count;
            }
        }

        assert(this->_tree->_root->_number > 0);
        assert(this->_tree->_root->_number == (int)(this->_tree->_preorder.size() - 1));

        this->_tree->_ninternals = internal_count;

        // If the tree has polytomies, then there are Node objects stord in
        // this->_tree->_nodes vector that have not yet been numbered. These
        // can be identified, because their _number is currently equal to -1.
        for (auto & nd : this->_tree->_nodes) {
            if (nd._number == -1)
                nd._number = curr_internal++;
        }
    }

    inline void TreeManip::store_splits(std::set<Split> & split_set) {    
        // Start by clearing and resizing all splits
        for (auto & nd : _tree->_nodes) {
            nd._split.resize(_tree->_nleaves);
        }

        // Now do a postorder traversal and add the bit corresponding
        // to the current node in its parent node's split
        for (auto nd : boost::adaptors::reverse(_tree->_preorder)) {
            if (nd->_left_child) {
                // add this internal node's split to splitset
                split_set.insert(nd->_split);
            }
            else {
                // set bit corresponding to this leaf node's number
                nd->_split.set_bit_at(nd->_number);
            }

            if (nd->_parent) {
                // parent's bits are the union of the bits set in all its children
                nd->_parent->_split.add_split(nd->_split);
            }
        }
    }

    inline std::string TreeManip::make_newick(unsigned precision, bool use_names) const {
        std::string newick;
        const boost::format tip_node_name_format( boost::str(boost::format("%%s:%%.%df") % precision) );
        const boost::format tip_node_number_format( boost::str(boost::format("%%d:%%.%df") % precision) );
        const boost::format internal_node_format( boost::str(boost::format("):%%.%df") % precision) );
        std::stack<Node *> node_stack;

        for (auto nd : this->_tree->_preorder) {
            if (nd->_left_child) {
                newick += "(";
                node_stack.push(nd);
            }
            else {
                if (use_names) {
                    std::string name = nd->_name;
                    if (! name.empty()) {
                        name = '\'' + name + '\'';
                    }
                    newick += boost::str(boost::format(tip_node_name_format)
                            % name
                            % nd->get_edge_length());
                }
                else {
                    newick += boost::str(boost::format(tip_node_number_format)
                            % (nd->_number + 1)
                            % nd->get_edge_length());
                }
                if (nd->_right_sib) {
                    newick += ",";
                }
                else {
                    Node * popped = (node_stack.empty() ? 0 : node_stack.top());
                    while (popped && !popped->_right_sib) {
                        node_stack.pop();
                        if (node_stack.empty()) {
                            newick += ")";
                            popped = 0;
                        }
                        else {
                            newick += boost::str(boost::format(internal_node_format) % popped->get_edge_length());
                            popped = node_stack.top();
                        }
                    }
                    if (popped && popped->_right_sib) {
                        node_stack.pop();
                        newick += boost::str(boost::format(internal_node_format) % popped->get_edge_length());
                        newick += ",";
                    }
                }
            }
        }

        return newick;
    }

}
