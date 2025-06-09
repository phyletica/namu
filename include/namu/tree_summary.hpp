#pragma once    

#include <set>
#include <map>
#include <vector>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <boost/format.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include <ncl/nxsmultiformat.h>
#include "split.hpp"
#include "tree_manip.hpp"
#include "tree_io.hpp"
#include "error.hpp"

namespace namu {

    class TreeSummary {
        public:
                                        TreeSummary();
                                        ~TreeSummary();

            void                        read_tree_file(const std::string file_name, unsigned skip);
            void                        show_summary() const;
            typename TreeManip::SharedPtr   get_tree_manip(unsigned index);
            typename Tree::SharedPtr    get_tree(unsigned index);
            std::string                 get_newick(unsigned index);
            void                        clear();

        private:

            Split::treemap_t            _treeIDs;
            std::vector<std::string>    _newicks;

        public:

            typedef std::shared_ptr< TreeSummary > SharedPtr;
    };

    inline TreeSummary::TreeSummary() { 
        // std::cout << "Constructing a TreeSummary" << std::endl;
    } 

    inline TreeSummary::~TreeSummary() {
        // std::cout << "Destroying a TreeSummary" << std::endl;
    }

    inline TreeManip::SharedPtr TreeSummary::get_tree_manip(unsigned index) {   
        if (index >= this->_newicks.size())
            throw NamuX("get_tree called with index >= number of stored trees");

        TreeIO tree_io;
        TreeManip::SharedPtrVector trees = tree_io.parse_from_newick(this->_newicks.at(index));

        return trees.at(0);
    }

    inline Tree::SharedPtr TreeSummary::get_tree(unsigned index) {   
        return this->get_tree_manip(index)->get_tree();
    }

    inline std::string TreeSummary::get_newick(unsigned index) { 
        if (index >= this->_newicks.size())
            throw NamuX("get_newick called with index >= number of stored trees");

        return this->_newicks.at(index);
    }

    inline void TreeSummary::clear() {  
        this->_treeIDs.clear();
        this->_newicks.clear();
    }

    inline void TreeSummary::read_tree_file(const std::string file_name, unsigned skip) {  
        TreeIO tree_io;
        TreeManip::SharedPtrVector trees;
        TreeManip::SharedPtr tm;
        Split::treeid_t split_set;

        trees = tree_io.parse_from_path(file_name, "nexus");
        this->clear();
        for (unsigned t = skip; t < trees.size(); ++t) {
            tm = trees.at(t);

            // store the newick tree description
            std::string newick = tm->make_newick(10, true);
            newick += ";";
            this->_newicks.push_back(newick);
            unsigned tree_index = (unsigned)this->_newicks.size() - 1;

            // store set of splits
            split_set.clear();
            tm->store_splits(split_set);

            // iterator iter will point to the value corresponding to key splitset
            // or to end (if splitset is not already a key in the map)
            Split::treemap_t::iterator iter = this->_treeIDs.lower_bound(split_set);

            if (iter == this->_treeIDs.end() || iter->first != split_set) {
                // splitset key not found in map, need to create an entry
                std::vector<unsigned> v(1, tree_index);  // vector of length 1 with only element set to tree_index
                this->_treeIDs.insert(iter, Split::treemap_t::value_type(split_set, v));
            }
            else {
                // splitset key was found in map, need to add this tree's index to vector
                iter->second.push_back(tree_index);
            }
        } // trees loop

        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for NCL documentation

        // MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
        // try {
        //     // nexus_reader.ReadFilepath(file_name.c_str(), MultiFormatReader::NEXUS_FORMAT);
        //     nexus_reader.ReadFilepath(file_name.c_str(), "relaxedphyliptree");
        // }
        // catch(...) {
        //     nexus_reader.DeleteBlocksFromFactories();
        //     throw;
        // }

        // int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
        // for (int i = 0; i < num_taxa_blocks; ++i) {
        //     this->clear();
        //     NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(i);
        //     std::string taxa_block_title = taxa_block->GetTitle();

        //     const unsigned n_trees_blocks = nexus_reader.GetNumTreesBlocks(taxa_block);
        //     for (unsigned j = 0; j < n_trees_blocks; ++j) {
        //         const NxsTreesBlock * trees_block = nexus_reader.GetTreesBlock(taxa_block, j);
        //         unsigned ntrees = trees_block->GetNumTrees();
        //         if (skip < ntrees) {
        //             //std::cout << "Trees block contains " << ntrees << " tree descriptions.\n";
        //             for (unsigned t = skip; t < ntrees; ++t) {
        //                 const NxsFullTreeDescription & d = trees_block->GetFullTreeDescription(t);

        //                 // store the newick tree description
        //                 std::string newick = d.GetNewick();
        //                 newick += ";";
        //                 this->_newicks.push_back(newick);
        //                 unsigned tree_index = (unsigned)this->_newicks.size() - 1;

        //                 // build the tree
        //                 // tm = tree_io.parse_from_newick(newick).at(0);
        //                 tm = this->get_tree_manip(tree_index);

        //                 // store set of splits
        //                 split_set.clear();
        //                 tm->store_splits(split_set);

        //                 // iterator iter will point to the value corresponding to key splitset
        //                 // or to end (if splitset is not already a key in the map)
        //                 Split::treemap_t::iterator iter = this->_treeIDs.lower_bound(split_set);

        //                 if (iter == this->_treeIDs.end() || iter->first != split_set) {
        //                     // splitset key not found in map, need to create an entry
        //                     std::vector<unsigned> v(1, tree_index);  // vector of length 1 with only element set to tree_index
        //                     this->_treeIDs.insert(iter, Split::treemap_t::value_type(split_set, v));
        //                 }
        //                 else {
        //                     // splitset key was found in map, need to add this tree's index to vector
        //                     iter->second.push_back(tree_index);
        //                 }
        //             } // trees loop
        //         } // if skip < ntrees
        //     } // TREES block loop
        // } // TAXA block loop

        // // No longer any need to store raw data from nexus file
        // nexus_reader.DeleteBlocksFromFactories();
    }

    inline void TreeSummary::show_summary() const {  
        // Produce some output to show that it works
        std::cout << boost::str(boost::format("\nRead %d trees from file") % this->_newicks.size()) << std::endl;

        // Show all unique topologies with a list of the trees that have that topology
        // Also create a map that can be used to sort topologies by their sample frequency
        typedef std::pair<unsigned, unsigned> sorted_pair_t;
        std::vector< sorted_pair_t > sorted;
        int t = 0;
        for (auto & key_value_pair : this->_treeIDs) {
            unsigned topology = ++t;
            unsigned ntrees = (unsigned)key_value_pair.second.size();
            sorted.push_back(std::pair<unsigned, unsigned>(ntrees,topology));
            std::cout << "Topology " << topology << " seen in these " << ntrees << " trees:" << std::endl << "  ";
            std::copy(key_value_pair.second.begin(), key_value_pair.second.end(), std::ostream_iterator<unsigned>(std::cout, " "));
            std::cout << std::endl;
        }

        // Show sorted histogram data
        std::sort(sorted.begin(), sorted.end());
        //unsigned npairs = (unsigned)sorted.size();
        std::cout << "\nTopologies sorted by sample frequency:" << std::endl;
        std::cout << boost::str(boost::format("%20s %20s") % "topology" % "frequency") << std::endl;
        for (auto & ntrees_topol_pair : boost::adaptors::reverse(sorted)) {
            unsigned n = ntrees_topol_pair.first;
            unsigned t = ntrees_topol_pair.second;
            std::cout << boost::str(boost::format("%20d %20d") % t % n) << std::endl;
        }
    }
}
