#pragma once

#include <cassert>
#include <memory>
#include <stack>
#include <queue>
#include <set>
#include <regex>
#include <limits>
#include <unordered_map>
#include <fstream>
#include <cmath>
#include <typeinfo>
#include <boost/range/adaptor/reversed.hpp>
#include <boost/format.hpp>
#include <ncl/nxsmultiformat.h>
#include "namu/tree.hpp"
#include "namu/tree_manip.hpp"
#include "namu/error.hpp"
#include "namu/string_util.hpp"

namespace namu {

    class TreeIO {

        public:
                                        TreeIO() {};
                                        ~TreeIO() {};

            TreeManip::SharedPtrVector  parse_from_newick(
                                            const std::string & newick,
                                            const std::string & ncl_file_format = "relaxedphyliptree",
                                            const double relative_ultrametricity_tolerance = 1e-6) const;
            TreeManip::SharedPtrVector  parse_from_path(
                                            const std::string & path,
                                            const std::string & ncl_file_format = "relaxedphyliptree",
                                            const double relative_ultrametricity_tolerance = 1e-6) const;
            TreeManip::SharedPtrVector  parse_from_stream(
                                            std::istream & tree_stream,
                                            const std::string & ncl_file_format = "relaxedphyliptree",
                                            const double relative_ultrametricity_tolerance = 1e-6) const;

        private:

            TreeManip::SharedPtr        get_from_ncl_tree_description(
                                            const NxsFullTreeDescription & tree_description,
                                            const NxsTaxaBlock * taxa_block,
                                            const double relative_ultrametricity_tolerance = 1e-6) const;
            double                      get_dist_from_root(const NxsSimpleNode * ncl_node) const;
            double                      get_height_from_edge_lengths(
                                            const NxsSimpleNode * ncl_node,
                                            const double max_root_to_tip_dist) const ;
            double                      get_max_root_to_tip_dist(
                                            NxsSimpleTree & ncl_simple_tree) const;
            double                      get_min_leaf_height(
                                            NxsSimpleTree & ncl_simple_tree) const;
            bool                        parsed_tree_is_ultrametric(
                                            const NxsSimpleTree & ncl_simple_tree,
                                            const double abs_tolerance) const;
            void                        parse_node_comments(
                                            const NxsSimpleNode * ncl_node,
                                            std::map<std::string, std::string> & comment_map) const;
            void                        process_child_nodes_from_ncl(
                                            TreeManip::SharedPtr tm,
                                            Node * parent_node,
                                            const NxsSimpleNode * parent_ncl_node,
                                            const NxsTaxaBlock * taxa_block,
                                            const std::unordered_map<const NxsSimpleNode *, int> node_to_preorder_index,
                                            const unsigned num_leaves,
                                            const double max_root_to_tip_dist,
                                            const double min_leaf_height,
                                            const double abs_tolerance,
                                            const bool using_height_comments) const;
            void                        process_leaf_node_from_ncl(
                                            Node * new_node,
                                            const NxsSimpleNode * ncl_node,
                                            const NxsTaxaBlock * taxa_block,
                                            const double max_root_to_tip_dist,
                                            const double min_leaf_height,
                                            const double abs_tolerance,
                                            const bool using_height_comments) const;
            void                        process_internal_node_from_ncl(
                                            TreeManip::SharedPtr tm,
                                            Node * new_node,
                                            const NxsSimpleNode * ncl_node,
                                            const double max_root_to_tip_dist,
                                            const bool using_height_comments) const;

    };

    inline TreeManip::SharedPtrVector TreeIO::parse_from_newick(
            const std::string & newick,
            const std::string & ncl_file_format,
            const double relative_ultrametricity_tolerance) const {
        std::istringstream in_stream(newick);
        try {
            return this->parse_from_stream(in_stream, ncl_file_format,
                    relative_ultrametricity_tolerance);
        }
        catch(const std::exception & e) {
            std::ostringstream msg;
            msg << "Problem parsing newick string: " << newick << std::endl;
            msg << "Exception caught: " << std::endl << typeid(e).name() << std::endl;
            msg << "Exception message: " << std::endl << e.what() << std::endl;
            throw NamuX(msg.str());
        }
    }

    inline TreeManip::SharedPtrVector TreeIO::parse_from_path(
            const std::string & path,
            const std::string & ncl_file_format,
            const double relative_ultrametricity_tolerance) const {
        std::ifstream in_stream;
        in_stream.open(path);
        if (! in_stream.is_open()) {
            std::ostringstream msg;
            msg << "Could not open file: " << path;
            throw NamuX(msg.str());
        }
        try {
            return this->parse_from_stream(in_stream, ncl_file_format,
                    relative_ultrametricity_tolerance);
        }
        catch(const std::exception & e) {
            std::ostringstream msg;
            msg << "Problem parsing file: " << path << std::endl;
            msg << "Exception caught: " << std::endl << typeid(e).name() << std::endl;
            msg << "Exception message: " << std::endl << e.what() << std::endl;
            throw NamuX(msg.str());
        }
    }

    inline TreeManip::SharedPtrVector TreeIO::parse_from_stream(
            std::istream & tree_stream,
            const std::string & ncl_file_format,
            const double relative_ultrametricity_tolerance) const {
        MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
        try {
            nexus_reader.ReadStream(tree_stream, ncl_file_format.c_str());
        }
        catch(...) {
            nexus_reader.DeleteBlocksFromFactories();
            throw;
        }
        unsigned int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
        if (num_taxa_blocks != 1) {
            std::ostringstream msg;
            msg << "Expecting 1 taxa block, but found " << num_taxa_blocks;
            throw NamuX(msg.str());
        }

        NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(0);

        unsigned int num_tree_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
        if (num_tree_blocks != 1) {
            std::ostringstream msg;
            msg << "Expecting 1 tree block, but found " << num_tree_blocks;
            throw NamuX(msg.str());
        }

        NxsTreesBlock * tree_block = nexus_reader.GetTreesBlock(taxa_block, 0);
        tree_block->ProcessAllTrees();
        unsigned int num_trees = tree_block->GetNumTrees();
        if (num_trees == 0) {
            std::ostringstream msg;
            msg << "Found no trees in tree block";
            throw NamuX(msg.str());
        }

        TreeManip::SharedPtrVector tree_manips;
        tree_manips.reserve(num_trees);
        for (unsigned i = 0; i < num_trees; ++i) {
            tree_manips.push_back(
                this->get_from_ncl_tree_description(
                    tree_block->GetFullTreeDescription(i),
                    taxa_block,
                    relative_ultrametricity_tolerance));
        }

        nexus_reader.DeleteBlocksFromFactories();
        return tree_manips;
    }

    inline TreeManip::SharedPtr TreeIO::get_from_ncl_tree_description(
            const NxsFullTreeDescription & tree_description,
            const NxsTaxaBlock * taxa_block,
            double relative_ultrametricity_tolerance) const {
        // Tree must be processed to create the NxsSimpleTree
        if (! tree_description.IsProcessed()) {
            throw NamuX("Input tree was not processed by NCL");
        }

        int default_int_edge_length = 0;
        double default_double_edge_length = 0.0;
        NxsSimpleTree ncl_simple_tree(tree_description,
                default_int_edge_length,
                default_double_edge_length);

        double max_root_to_tip_dist = this->get_max_root_to_tip_dist(ncl_simple_tree);
        double abs_tolerance = max_root_to_tip_dist * relative_ultrametricity_tolerance;
        // bool tree_is_ultrametric = this->parsed_tree_is_ultrametric(
        //         ncl_simple_tree,
        //         abs_tolerance);
        // if (! tree_is_ultrametric) {
        //     throw NamuX("Input tree not ultrametric");
        // }
        // const std::vector<NxsSimpleNode *> & leaves = ncl_simple_tree.GetLeavesRef();
        unsigned num_leaves = ncl_simple_tree.GetLeavesRef().size();
        // int num_leaves_int = leaves.size();
        // std::vector<std::string> leaf_labels;
        // leaf_labels.reserve(num_leaves);
        // for (auto leaf_node : leaves) {
        //     NxsString label = taxa_block->GetTaxonLabel(leaf_node->GetTaxonIndex());
        //     leaf_labels.push_back(label);
        // }
        // // Sort labels to ensure we always give the same label the same
        // // leaf node number 
        // std::sort(leaf_labels.begin(), leaf_labels.end());
        // // Create map of labels to node numbers
        // std::unordered_map<std::string, int> leaf_label_to_num_map;
        // leaf_label_to_num_map.reserve(num_leaves);
        // for (int i = 0; i < num_leaves_int; ++i) {
        //     assert(leaf_label_to_num_map.count(leaf_labels.at(i)) == 0);
        //     leaf_label_to_num_map[leaf_labels.at(i)] = i;
        // }

        // Create map of nodes to preorder indices
        // Making sure all persistant nodes (leaves and root) get the first
        // 0--num_leaves indices
        const std::vector<const NxsSimpleNode *> preorder_nodes = ncl_simple_tree.GetPreorderTraversal();
        unsigned num_nodes = preorder_nodes.size();
        std::unordered_map<const NxsSimpleNode *, int> node_to_preorder_index;
        node_to_preorder_index.reserve(preorder_nodes.size());
        unsigned next_leaf_index = 0;
        unsigned next_internal_index = num_leaves;
        for (unsigned i = 0; i < preorder_nodes.size(); ++i) {
            if (preorder_nodes.at(i)->IsTip()) {
                auto r = node_to_preorder_index.insert({preorder_nodes.at(i), next_leaf_index});
                assert(r.second);
                ++next_leaf_index;
            } else {
                auto r = node_to_preorder_index.insert({preorder_nodes.at(i), next_internal_index});
                assert(r.second);
                ++next_internal_index;
            }
        }
        assert(next_leaf_index == num_leaves);
        assert(next_internal_index == num_nodes);

        const NxsSimpleNode * ncl_root = ncl_simple_tree.GetRootConst();
        assert(node_to_preorder_index.at(ncl_root) == (int)num_leaves);
        NxsSimpleEdge root_edge = ncl_root->GetEdgeToParent();
        std::map<std::string, std::string> root_info;
        for (auto nxs_comment : root_edge.GetUnprocessedComments()) {
            std::string comment = str_util::strip(nxs_comment.GetText());
            if (str_util::startswith(comment, "&")) {
                std::string root_info_str = comment.substr(1);
                str_util::parse_map(
                        root_info_str,
                        root_info,
                        ',',
                        '=');
            }
        }
        // Preferably, we want to get the node heights and node height
        // indices from the metadata in the comments.
        // We can check the root if this information exists.
        // If so, we shoud throw an error if any other nodes lack these
        // data.
        // If the root does not have these data, we could still parse the
        // tree and calculate heights from the edge lengths; without height
        // indices, the resulting tree would have no shared node heights.
        bool using_height_comments = false;
        if (root_info.count("height_index") > 0) {
            using_height_comments = true;
        }

        double min_leaf_height = 0.0;
        if (using_height_comments) {
            min_leaf_height = this->get_min_leaf_height(ncl_simple_tree);
        }

        // Prepare data members
        std::shared_ptr<TreeManip> tm = std::make_shared<TreeManip>();
        tm->_tree.reset(new Tree());
        tm->_div_events.clear();
        // _tree->_ninternals is updated by call to renumber_internals
        tm->_tree->_nleaves = num_leaves;

        unsigned max_num_nodes = ((2 * num_leaves) - 1);
        tm->_tree->_nodes.resize(max_num_nodes);
        tm->_div_events.resize(max_num_nodes - num_leaves);
        for (auto & nd : tm->_tree->_nodes) {
            nd._number = -1;
        }

        try {
            // Make the root node
            int node_index = node_to_preorder_index.at(ncl_root);
            assert(node_index == (int)num_leaves);
            Node * root = &tm->_tree->_nodes.at(node_index);
            tm->_tree->_root = root;
            this->process_internal_node_from_ncl(
                    tm,
                    root,
                    ncl_root,
                    max_root_to_tip_dist,
                    using_height_comments);
            if (! using_height_comments) {
                tm->_div_events.at(node_index - num_leaves).push_back(root);
            }

            // Recursively add nodes down the tree toward the leaves
            this->process_child_nodes_from_ncl(
                    tm,
                    root,
                    ncl_root,
                    taxa_block,
                    node_to_preorder_index,
                    num_leaves,
                    max_root_to_tip_dist,
                    min_leaf_height,
                    abs_tolerance,
                    using_height_comments);
        }
        catch(const std::exception & e) {
            tm->clear();
            throw e;
        }

        bool gap_in_divs = (! tm->no_gaps_in_div_events());
        if (gap_in_divs) {
            std::ostringstream msg;
            msg << "ERROR: Found gap in height indices in newick tree";
            throw NamuX(msg.str());
        }
        tm->refresh_preorder();
        tm->refresh_levelorder();
        tm->renumber_internals();
        return tm;
    }

    inline double TreeIO::get_dist_from_root(const NxsSimpleNode * ncl_node) const {
        double dist = 0.0;
        const NxsSimpleNode * parent = ncl_node->GetEdgeToParent().GetParent();
        while (parent) {
            dist += ncl_node->GetEdgeToParent().GetDblEdgeLen();
            ncl_node = parent;
            parent = ncl_node->GetEdgeToParent().GetParent();
        }
        return dist;
    }

    inline double TreeIO::get_height_from_edge_lengths(
            const NxsSimpleNode * ncl_node,
            const double max_root_to_tip_dist) const {
        double dist_from_root = this->get_dist_from_root(ncl_node);
        assert(max_root_to_tip_dist >= dist_from_root);
        double height = max_root_to_tip_dist - dist_from_root;
        return height;
    }

    inline double TreeIO::get_max_root_to_tip_dist(NxsSimpleTree & ncl_simple_tree) const {
        const std::vector<NxsSimpleNode *> & leaves = ncl_simple_tree.GetLeavesRef();
        double max_dist = -1.0;
        for (auto leaf_node : leaves) {
            double dist = this->get_dist_from_root(leaf_node);
            if (dist > max_dist) {
                max_dist = dist;
            }
        }
        return max_dist;
    }

    inline double TreeIO::get_min_leaf_height(NxsSimpleTree & ncl_simple_tree) const {
        const std::vector<NxsSimpleNode *> & leaves = ncl_simple_tree.GetLeavesRef();
        double min_height = std::numeric_limits<double>::infinity();
        for (auto leaf_node : leaves) {
            std::map<std::string, std::string> comment_map;
            this->parse_node_comments(leaf_node, comment_map);
            if (comment_map.count("height") < 1) {
                throw NamuX(
                    "TreeIO::get_min_leaf_height: A leaf didn't have a height comment"
                );
            }
            double height;
            std::stringstream h_converter(comment_map["height"]);
            if (! (h_converter >> height)) {
                throw NamuX("could not convert node height \'" +
                        h_converter.str() + "\'");
            }
            if (height < min_height) {
                min_height = height;
            }
        }
        return min_height;
    }

    inline bool TreeIO::parsed_tree_is_ultrametric(
            const NxsSimpleTree & ncl_simple_tree,
            const double abs_tolerance) const {
        std::vector< std::vector<double> > dists_to_mrca = ncl_simple_tree.GetDblPathDistances(true);
        unsigned int num_leaves = dists_to_mrca.size();
        assert(dists_to_mrca.at(0).size() == num_leaves);

        // Check if distance from leaf i to the MRCA of i and j is (almost)
        // equal to the distance from leaf j to the MRCA of i and j. If
        // this is true for all pairs of tips, then the tree is
        // ultrametric.
        for (unsigned int i = 0; i < num_leaves - 1; ++i) {
            for (unsigned int j = i + 1; j < num_leaves; ++j) {
                double height_diff = fabs(dists_to_mrca.at(i).at(j) - dists_to_mrca.at(j).at(i));
                if (height_diff > abs_tolerance) {
                    return false;
                }
            }
        }
        return true;
    }

    inline void TreeIO::parse_node_comments(
            const NxsSimpleNode * ncl_node,
            std::map<std::string, std::string> & comment_map) const {
        NxsSimpleEdge ncl_edge = ncl_node->GetEdgeToParent();
        for (auto nxs_comment : ncl_edge.GetUnprocessedComments()) {
            std::string raw_comment = str_util::strip(nxs_comment.GetText());
            if (str_util::startswith(raw_comment, "&")) {
                std::string comment = raw_comment.substr(1);
                str_util::parse_map(
                        comment,
                        comment_map,
                        ',',
                        '=');
            }
        }
    }

    inline void TreeIO::process_child_nodes_from_ncl(
            TreeManip::SharedPtr tm,
            Node * parent_node,
            const NxsSimpleNode * parent_ncl_node,
            const NxsTaxaBlock * taxa_block,
            const std::unordered_map<const NxsSimpleNode *, int> node_to_preorder_index,
            const unsigned num_leaves,
            const double max_root_to_tip_dist,
            const double min_leaf_height,
            const double abs_tolerance,
            const bool using_height_comments) const {
        for (auto child_ncl_node : parent_ncl_node->GetChildren()) {
            int node_index = node_to_preorder_index.at(child_ncl_node);
            Node * new_node = &tm->_tree->_nodes.at(node_index);
            if (child_ncl_node->IsTip()) {
                // This is a leaf
                // We give leaves numbers now, but wait until
                // renumber_internals call for internal nodes
                new_node->_number = node_index;
                this->process_leaf_node_from_ncl(
                        new_node,
                        child_ncl_node,
                        taxa_block,
                        max_root_to_tip_dist,
                        min_leaf_height,
                        abs_tolerance,
                        using_height_comments);
            }
            else {
                this->process_internal_node_from_ncl(
                        tm,
                        new_node,
                        child_ncl_node,
                        max_root_to_tip_dist,
                        using_height_comments);
                if (! using_height_comments) {
                    tm->_div_events.at(node_index - num_leaves).push_back(new_node);
                }
            }
            if (parent_node->has_child()) {
                parent_node->get_rightmost_child()->_right_sib = new_node;
                new_node->_parent = parent_node;
            }
            else {
                parent_node->_left_child = new_node;
                new_node->_parent = parent_node;
            }
            if (! child_ncl_node->IsTip()) {
                // Need to continue to recurse down the tree, toward the
                // tips
                this->process_child_nodes_from_ncl(
                        tm,
                        new_node,
                        child_ncl_node,
                        taxa_block,
                        node_to_preorder_index,
                        num_leaves,
                        max_root_to_tip_dist,
                        min_leaf_height,
                        abs_tolerance,
                        using_height_comments);
            }
        }
    }

    inline void TreeIO::process_leaf_node_from_ncl(
            Node * new_node,
            const NxsSimpleNode * ncl_node,
            const NxsTaxaBlock * taxa_block,
            const double max_root_to_tip_dist,
            const double min_leaf_height,
            const double abs_tolerance,
            const bool using_height_comments) const {
        assert(ncl_node->IsTip());
        std::map<std::string, std::string> comment_map;
        this->parse_node_comments(ncl_node, comment_map);
        double height;
        if (using_height_comments) {
            std::stringstream h_converter(comment_map["height"]);
            if (! (h_converter >> height)) {
                throw NamuX("could not convert node height \'" +
                        h_converter.str() + "\'");
            }
        }
        else {
            // Need to get height from edge lengths
            height = this->get_height_from_edge_lengths(ncl_node, max_root_to_tip_dist);
            // Check to see if leaf height should be equal to zero (or the
            // min_leaf_height, but that should be zero if we are getting
            // heights from branch lengths
            if (fabs(height - min_leaf_height) <= abs_tolerance) {
                height = min_leaf_height;
            }
        }

        new_node->_height = height;
        new_node->_name = taxa_block->GetTaxonLabel(ncl_node->GetTaxonIndex());
        // new_node->extract_data_from_node_comments(comment_map);
    }

    inline void TreeIO::process_internal_node_from_ncl(
            TreeManip::SharedPtr tm,
            Node * new_node,
            const NxsSimpleNode * ncl_node,
            const double max_root_to_tip_dist,
            const bool using_height_comments) const {
        assert(! ncl_node->IsTip());
        std::map<std::string, std::string> comment_map;
        this->parse_node_comments(ncl_node, comment_map);
        double height;
        unsigned int height_index;
        if (using_height_comments) {
            std::stringstream h_converter(comment_map["height"]);
            if (! (h_converter >> height)) {
                throw NamuX("could not convert node height \'" +
                        h_converter.str() + "\'");
            }
            std::stringstream i_converter(comment_map["height_index"]);
            if (! (i_converter >> height_index)) {
                throw NamuX("could not convert node height index \'" +
                        i_converter.str() + "\'");
            }
        }
        else {
            // Need to get height from edge lengths
            height = this->get_height_from_edge_lengths(ncl_node, max_root_to_tip_dist);
        }

        if (using_height_comments) {
            if (tm->_div_events.at(height_index).size() > 0) {
                // Sanity check to make sure height indices are consistent
                // across tree
                if (height != tm->_div_events.at(height_index).at(0)->_height) {
                    std::ostringstream msg;
                    msg << "ERROR: create_internal_node: "
                            << "Height index " << height_index
                            << " has multiple heights in tree";
                    throw NamuX(msg.str());
                }
            }
            tm->_div_events.at(height_index).push_back(new_node);
        }

        new_node->_height = height;
        new_node->_name = ncl_node->GetName();
        // new_node->extract_data_from_node_comments(comment_map);
    }
}
