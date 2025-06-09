#pragma once

#include <tuple>
#include <limits>
#include <cmath>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include "namu/genetic_code.hpp"
#include "namu/data_type.hpp"
#include "namu/error.hpp"

namespace namu {

    class Partition {
        public:
            typedef std::match_results<std::string::const_iterator>::const_reference    regex_match_t;
            typedef std::tuple<unsigned, unsigned, unsigned, unsigned>                  subset_range_t;
            typedef std::vector<subset_range_t>                                         partition_t;
            typedef std::vector<DataType>                                               datatype_vect_t;
            typedef std::vector<unsigned>                                               subset_sizes_vect_t;
            typedef std::vector<std::string>                                            subset_names_vect_t;
            typedef std::shared_ptr<Partition>                                          SharedPtr;

                                                        Partition();
                                                        ~Partition();
        
            unsigned                                    get_num_sites() const;
            unsigned                                    get_num_subsets() const;
            std::string                                 get_subset_name(unsigned subset) const;
        
            const partition_t &                         get_subset_range_vect() const;
        
            unsigned                                    find_subset_by_name(const std::string & subset_name) const;
            unsigned                                    find_subset_for_site(unsigned site_index) const;
            bool                                        site_in_subset(unsigned site_index, unsigned subset_index) const;
            DataType                                    get_data_type_for_subset(unsigned subset_index) const;
            const datatype_vect_t &                     get_subset_data_types() const;
        
            unsigned                                    num_sites_in_subset(unsigned subset_index) const;
            subset_sizes_vect_t                         calc_subset_sizes() const;

            void                                        default_partition(unsigned nsites = std::numeric_limits<unsigned>::max());
            void                                        parse_subset_definition(std::string & s);
            void                                        finalize(unsigned nsites);

            void                                        clear();

        private:

            int                                         extract_int_from_regex_match(regex_match_t s, unsigned min_value);
            void                                        add_subset_range(unsigned subset_index, std::string range_definition);
            void                                        add_subset(unsigned subset_index, std::string subset_definition);

            unsigned                                    _num_sites;
            unsigned                                    _num_subsets;
            subset_names_vect_t                         _subset_names;
            partition_t                                 _subset_ranges;
            datatype_vect_t                             _subset_data_types;

            const unsigned                              _infinity;
    };
    
    inline Partition::Partition() : _infinity(std::numeric_limits<unsigned>::max()) {
        //std::cout << "Constructing a Partition" << std::endl;
        clear();
    }  

    inline Partition::~Partition() {
        //std::cout << "Destroying a Partition" << std::endl;
    }

    inline unsigned Partition::get_num_sites() const {
        return this->_num_sites;
    }
    
    inline unsigned Partition::get_num_subsets() const {
        return this->_num_subsets;
    }
    
    inline std::string Partition::get_subset_name(unsigned subset) const {
        assert(subset < this->_num_subsets);
        return this->_subset_names[subset];
    }
    
    inline const Partition::partition_t & Partition::get_subset_range_vect() const {
        return this->_subset_ranges;
    }
    
    inline DataType Partition::get_data_type_for_subset(unsigned subset_index) const {
        assert(subset_index < this->_subset_data_types.size());
        return this->_subset_data_types[subset_index];
    }

    inline const std::vector<DataType> & Partition::get_subset_data_types() const {
        return this->_subset_data_types;
    }

    inline unsigned Partition::find_subset_by_name(const std::string & subset_name) const {
        auto iter = std::find(this->_subset_names.begin(), this->_subset_names.end(), subset_name);
        if (iter == this->_subset_names.end())
            throw NamuX(boost::format("Specified subset name \"%s\" not found in partition") % subset_name);
        return (unsigned)std::distance(this->_subset_names.begin(), iter);
    }

    inline unsigned Partition::find_subset_for_site(unsigned site_index) const {
        for (auto & t : this->_subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            bool inside_range = site_index >= begin_site && site_index <= end_site;
            if (inside_range && (site_index - begin_site) % stride == 0)
                return site_subset;
        }
        throw NamuX(boost::format("Site %d not found in any subset of partition") % (site_index + 1));
    }
    
    inline bool Partition::site_in_subset(unsigned site_index, unsigned subset_index) const {
        unsigned which_subset = this->find_subset_for_site(site_index);
        return (which_subset == subset_index ? true : false);
    }
    
    inline unsigned Partition::num_sites_in_subset(unsigned subset_index) const {
        unsigned nsites = 0;
        for (auto & t : this->_subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            if (site_subset == subset_index) {
                unsigned n = end_site - begin_site + 1;
                nsites += (unsigned)(floor(n/stride)) + (n % stride == 0 ? 0 : 1);
            }
        }
        return nsites;
    }
    
    inline std::vector<unsigned> Partition::calc_subset_sizes() const {
        assert(this->_num_sites > 0); // only makes sense to call this function after subsets are defined
        std::vector<unsigned> nsites_vect(this->_num_subsets, 0);
        for (auto & t : this->_subset_ranges) {
            unsigned begin_site = std::get<0>(t);
            unsigned end_site = std::get<1>(t);
            unsigned stride = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            unsigned hull = end_site - begin_site + 1;
            unsigned n = (unsigned)(floor(hull/stride)) + (hull % stride == 0 ? 0 : 1);
            nsites_vect[site_subset] += n;
        }
        return nsites_vect;
    }
    
    inline void Partition::clear() {
        this->_num_sites = 0;
        this->_num_subsets = 1;
        this->_subset_data_types.clear();
        this->_subset_data_types.push_back(DataType());
        this->_subset_names.clear();
        this->_subset_names.push_back("default");
        this->_subset_ranges.clear();
        this->_subset_ranges.push_back(std::make_tuple(1, this->_infinity, 1, 0));
    }
    
    inline void Partition::parse_subset_definition(std::string & s) {
        std::vector<std::string> v;
        
        // first separate part before colon (stored in v[0]) from part after colon (stored in v[1])
        boost::split(v, s, boost::is_any_of(":"));
        if (v.size() != 2)
            throw NamuX("Expecting exactly one colon in partition subset definition");

        std::string before_colon = v[0];
        std::string subset_definition = v[1];

        // now see if before_colon contains a data type specification in square brackets
        const char * pattern_string = R"((.+?)\s*(\[(\S+?)\])*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(before_colon, match_obj, re);
        if (! matched) {
            throw NamuX(boost::format("Could not interpret \"%s\" as a subset label with optional data type in square brackets") % before_colon);
        }
        
        // match_obj always yields 2 strings that can be indexed using the operator[] function
        // match_obj[0] equals entire subset label/type string (e.g. "rbcL[codon:standard]")
        // match_obj[1] equals the subset label (e.g. "rbcL")
        
        // Two more elements will exist if the user has specified a data type for this partition subset
        // match_obj[2] equals data type inside square brackets (e.g. "[codon:standard]")
        // match_obj[3] equals data type only (e.g. "codon:standard")
        
        std::string subset_name = match_obj[1].str();
        DataType dt;    // nucleotide by default
        std::string data_type = "nucleotide";
        if (match_obj.size() == 4 && match_obj[3].length() > 0) {
            data_type = match_obj[3].str();
            boost::to_lower(data_type);

            // check for comma plus genetic code in case of codon
            std::regex re(R"(codon\s*,\s*(\S+))");
            std::smatch m;
            if (std::regex_match(data_type, m, re)) {
                dt.set_codon();
                std::string genetic_code_name = m[1].str();
                dt.set_genetic_code_from_name(genetic_code_name);
            }
            else if (data_type == "codon") {
                dt.set_codon();  // assumes standard genetic code
            }
            else if (data_type == "protein") {
                dt.set_protein();
            }
            else if (data_type == "nucleotide") {
                dt.set_nucleotide();
            }
            else if (data_type == "standard") {
                dt.set_standard();
            }
            else {
                throw NamuX(boost::format("Data type \"%s\" specified for subset(s) \"%s\" is invalid: must be either nucleotide, codon, protein, or standard") % data_type % subset_name);
            }
        }

        // Remove default subset if there is one
        unsigned end_site = std::get<1>(this->_subset_ranges[0]);
        if (this->_num_subsets == 1 && end_site == this->_infinity) {
            this->_subset_names.clear();
            this->_subset_data_types.clear();
            this->_subset_ranges.clear();
        }
        else if (subset_name == "default") {
            throw NamuX("Cannot specify \"default\" partition subset after already defining other subsets");
        }
        this->_subset_names.push_back(subset_name);
        this->_subset_data_types.push_back(dt);
        this->_num_subsets = (unsigned)this->_subset_names.size();
        this->add_subset(this->_num_subsets - 1, subset_definition);

        std::cout << boost::str(boost::format("Partition subset %s comprises sites %s and has type %s") % subset_name % subset_definition % data_type) << std::endl;
    }
    
    inline void Partition::add_subset(unsigned subset_index, std::string subset_definition) {
        std::vector<std::string> parts;
        boost::split(parts, subset_definition, boost::is_any_of(","));
        for (auto subset_component : parts) {
            this->add_subset_range(subset_index, subset_component);
        }
    }
    
    inline void Partition::add_subset_range(unsigned subset_index, std::string range_definition) {
        // match patterns like these: "1-.\3" "1-1000" "1001-."
        const char * pattern_string = R"((\d+)\s*(-\s*([0-9.]+)(\\\s*(\d+))*)*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(range_definition, match_obj, re);
        if (!matched) {
            throw NamuX(boost::format("Could not interpret \"%s\" as a range of site indices") % range_definition);
        }
        
        // match_obj always yields 6 strings that can be indexed using the operator[] function
        // match_obj[0] equals entire site_range (e.g. "1-.\3")
        // match_obj[1] equals beginning site index (e.g. "1")
        // match_obj[2] equals everything after beginning site index (e.g. "-.\3")
        // match_obj[3] equals "" or ending site index (e.g. ".")
        // match_obj[4] equals "" or everything after ending site index (e.g. "\3")
        // match_obj[5] equals "" or step value (e.g. "3")
        int ibegin = this->extract_int_from_regex_match(match_obj[1], 1);
        int iend   = this->extract_int_from_regex_match(match_obj[3], ibegin);
        int istep  = this->extract_int_from_regex_match(match_obj[5], 1);
        
        // record the triplet
        this->_subset_ranges.push_back(std::make_tuple(ibegin, iend, istep, subset_index));
        
        // determine last site in subset
        unsigned last_site_in_subset = iend - ((iend - ibegin) % istep);
        if (last_site_in_subset > this->_num_sites) {
            this->_num_sites = last_site_in_subset;
        }
    }
    
    inline int Partition::extract_int_from_regex_match(regex_match_t s, unsigned min_value) {
        int int_value = min_value;
        if (s.length() > 0) {
            std::string str_value = s.str();
            try {
                int_value = std::stoi(str_value);
            }
            catch(const std::invalid_argument & e) {
                throw NamuX(boost::format("Could not interpret \"%s\" as a number in partition subset definition") % s.str());
            }
            
            // sanity check
            if (int_value < (int)min_value) {
                throw NamuX(boost::format("Value specified in partition subset definition (%d) is lower than minimum value (%d)") % int_value % min_value);
            }
        }
        return int_value;
    }
    
    inline void Partition::finalize(unsigned nsites) {
        if (this->_num_sites == 0) {
            this->default_partition(nsites);
            return;
        }

        // First sanity check:
        //   nsites is the number of sites read in from a data file;
        //   _num_sites is the maximum site index specified in any partition subset.
        //   These two numbers should be the same.
        if (_num_sites != nsites) {
            throw NamuX(boost::format("Number of sites specified by the partition (%d) does not match actual number of sites (%d)") % _num_sites % nsites);
        }
        
        // Second sanity check: ensure that no sites were left out of all partition subsets
        // Third sanity check: ensure that no sites were included in more than one partition subset
        std::vector<int> tmp(nsites, -1);   // begin with -1 for all sites
        for (auto & t : this->_subset_ranges) {
            unsigned begin_site  = std::get<0>(t);
            unsigned end_site    = std::get<1>(t);
            unsigned stride  = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            for (unsigned s = begin_site; s <= end_site; s += stride) {
                if (tmp[s-1] != -1)
                    throw NamuX("Some sites were included in more than one partition subset");
                else
                    tmp[s-1] = site_subset;
            }
        }
        if (std::find(tmp.begin(), tmp.end(), -1) != tmp.end()) {
            throw NamuX("Some sites were not included in any partition subset");
        }
        tmp.clear();
    }
    
    inline void Partition::default_partition(unsigned nsites) {
        this->clear();
        this->_num_sites = nsites;
        this->_num_subsets = 1;
        this->_subset_ranges[0] = std::make_tuple(1, nsites, 1, 0);
    }
}
