#pragma once

#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/range/adaptor/reversed.hpp>
#include "libhmsbeagle/beagle.h"
#include "namu/tree.hpp"
// #include "namu/tree_manip.hpp"
#include "namu/data.hpp"
#include "namu/model.hpp"
#include "namu/error.hpp"

namespace namu {

    class Likelihood {
        public:
                                                    Likelihood();
                                                    ~Likelihood();

            void                                    set_prefer_gpu(bool prefer_gpu);
            void                                    set_ambiguity_equals_missing(bool ambig_equals_missing);

            bool                                    using_stored_data() const;
            void                                    use_stored_data(bool using_data);
            void                                    use_underflow_scaling(bool do_scaling);
            bool                                    using_underflow_scaling() const;

            std::string                             beagle_lib_version() const;
            std::string                             available_resources() const;
            std::string                             used_resources() const;

            void                                    init_beagle_lib();
            void                                    finalize_beagle_lib(bool use_exceptions);

            double                                  calc_log_likelihood(Tree::SharedPtr t);

            Data::SharedPtr                         get_data();
            void                                    set_data(Data::SharedPtr d);

            Model::SharedPtr                        get_model();
            void                                    set_model(Model::SharedPtr m);

            void                                    clear();
            
            unsigned                                calc_num_edges_in_fully_resolved_tree() const;
            unsigned                                calc_num_internals_in_fully_resolved_tree() const;

        private:
        
            struct InstanceInfo {
                int handle;
                int resourcenumber;
                std::string resourcename;
                unsigned nstates;
                unsigned nratecateg;
                unsigned npatterns;
                unsigned partial_offset;
                unsigned t_matrix_offset;
                bool invarmodel;
                std::vector<unsigned> subsets;
                
                InstanceInfo() : handle(-1), resourcenumber(-1), resourcename(""), nstates(0), nratecateg(0), npatterns(0), partial_offset(0), t_matrix_offset(0), invarmodel(false) {}
            };

            typedef std::pair<unsigned, int>        instance_pair_t;

            unsigned                                get_scaler_index(Node * nd, InstanceInfo & info) const;
            unsigned                                get_partial_index(Node * nd, InstanceInfo & info) const;
            unsigned                                get_t_matrix_index(Node * nd, InstanceInfo & info, unsigned subset_index) const;
            unsigned                                get_t_matrix_index_from_node_num(int node_num, InstanceInfo & info, unsigned subset_index) const;
            void                                    update_instance_map(instance_pair_t & p, unsigned subset);
            void                                    new_instance(unsigned nstates, int nrates, std::vector<unsigned> & subset_indices);
            void                                    set_tip_states();
            void                                    set_tip_partials();
            void                                    set_pattern_partition_assignments();
            void                                    set_pattern_weights();
            void                                    set_among_site_rate_heterogeneity();
            void                                    set_model_rate_matrix();
            void                                    add_operation(InstanceInfo & info, Node * nd, Node * lchild, Node * rchild, unsigned subset_index);
            void                                    queue_partials_recalculation(Node * nd, Node * lchild, Node * rchild);
            // void                                    queue_partials_recalculation(Node * nd, Node * lchild, Node * rchild, Node * polytomy = 0);
            void                                    queue_t_matrix_recalculation(Node * nd);
            void                                    define_operations(Tree::SharedPtr t);
            void                                    update_transition_matrices();
            void                                    calculate_partials();
            double                                  calc_instance_log_likelihood(InstanceInfo & inst, Tree::SharedPtr t);


            std::vector<InstanceInfo>               _instances;
            std::map<int, std::string>              _beagle_error;
            std::map<int, std::vector<int> >        _operations;
            std::map<int, std::vector<int> >        _pmatrix_index;
            std::map<int, std::vector<double> >     _edge_lengths;
            std::map<int, std::vector<int> >        _eigen_indices;
            std::map<int, std::vector<int> >        _category_rate_indices;
            double                                  _relrate_normalizing_constant;

            std::vector<int>                        _subset_indices;
            std::vector<int>                        _parent_indices;
            std::vector<int>                        _child_indices;
            std::vector<int>                        _t_matrix_indices;
            std::vector<int>                        _weights_indices;
            std::vector<int>                        _freqs_indices;
            std::vector<int>                        _scaling_indices;
        
            Data::SharedPtr                         _data;
            Model::SharedPtr                         _model;
            unsigned                                _ntaxa;
            bool                                    _prefer_gpu;
            bool                                    _ambiguity_equals_missing;
            bool                                    _underflow_scaling;
            bool                                    _using_data;

            // std::vector<Node *>                     _polytomy_helpers;
            // std::map<int, std::vector<int> >        _polytomy_map;
            // std::vector<double>                     _identity_matrix;

        public:
            typedef std::shared_ptr< Likelihood >   SharedPtr;
    }; 

    inline Likelihood::Likelihood() {
        //std::cout << "Constructing a Likelihood" << std::endl;
        this->clear();
    }

    inline Likelihood::~Likelihood() {
        //std::cout << "Destroying a Likelihood" << std::endl;
        this->finalize_beagle_lib(false);
        this->clear();
    }
    
    inline unsigned Likelihood::calc_num_edges_in_fully_resolved_tree() const {
        assert(this->_ntaxa > 0);
        // return (_rooted ? (2*_ntaxa - 2) : (2*_ntaxa - 3));
        return ((2 * this->_ntaxa) - 2);
    }
    
    inline unsigned Likelihood::calc_num_internals_in_fully_resolved_tree() const {
        assert(this->_ntaxa > 0);
        //return (_rooted ? (_ntaxa - 1) : (_ntaxa - 2));
        return (this->_ntaxa - 1);
    }

    inline void Likelihood::finalize_beagle_lib(bool use_exceptions) {
        // Close down all BeagleLib instances if active
        for (auto info : this->_instances) {
            if (info.handle >= 0) {
                int code = beagleFinalizeInstance(info.handle);
                if (code != 0) {
                    if (use_exceptions)
                        throw NamuX(boost::format("Likelihood failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code]);
                    else
                        std::cerr << boost::format("Likelihood destructor failed to finalize BeagleLib instance. BeagleLib error code was %d (%s).") % code % _beagle_error[code] << std::endl;
                }
            }
        }
        this->_instances.clear();
    }

    inline void Likelihood::clear() {
        this->finalize_beagle_lib(true);
        
        this->_ntaxa                      = 0;
        this->_prefer_gpu                 = false;
        this->_ambiguity_equals_missing   = true;
        this->_underflow_scaling          = false;
        this->_using_data                 = true;
        this->_data                       = nullptr;
        this->_model                      = Model::SharedPtr(new Model());
        
        this->_operations.clear();
        this->_pmatrix_index.clear();
        this->_edge_lengths.clear();
        this->_eigen_indices.clear();
        this->_category_rate_indices.clear();
        this->_relrate_normalizing_constant = 1.0;
        this->_subset_indices.assign(1, 0);
        this->_parent_indices.assign(1, 0);
        this->_child_indices.assign(1, 0);
        this->_t_matrix_indices.assign(1, 0);
        this->_weights_indices.assign(1, 0);
        this->_freqs_indices.assign(1, 0);
        this->_scaling_indices.assign(1, 0);
        // this->_identity_matrix.assign(1, 0.0);
        
        // Store BeagleLib error codes so that useful
        // error messages may be provided to the user
        this->_beagle_error.clear();
        this->_beagle_error[0]  = std::string("success");
        this->_beagle_error[-1] = std::string("unspecified error");
        this->_beagle_error[-2] = std::string("not enough memory could be allocated");
        this->_beagle_error[-3] = std::string("unspecified exception");
        this->_beagle_error[-4] = std::string("the instance index is out of range, or the instance has not been created");
        this->_beagle_error[-5] = std::string("one of the indices specified exceeded the range of the array");
        this->_beagle_error[-6] = std::string("no resource matches requirements");
        this->_beagle_error[-7] = std::string("no implementation matches requirements");
        this->_beagle_error[-8] = std::string("floating-point range exceeded");
    }

    inline std::string Likelihood::beagle_lib_version() const {
        return std::string(beagleGetVersion());
    } 
    
    inline std::string Likelihood::available_resources() const {
        BeagleResourceList * rsrcList = beagleGetResourceList();
        std::string s;
        for (int i = 0; i < rsrcList->length; ++i) {
            std::string desc = rsrcList->list[i].description;
            boost::trim(desc);
            if (desc.size() > 0)
                s += boost::str(boost::format("  resource %d: %s (%s)\n") % i % rsrcList->list[i].name % desc);
            else
                s += boost::str(boost::format("  resource %d: %s\n") % i % rsrcList->list[i].name);
        }
        boost::trim_right(s);
        return s;
    }
    
    inline std::string Likelihood::used_resources() const {
        std::string s;
        for (unsigned i = 0; i < this->_instances.size(); i++) {
            s += boost::str(boost::format("  instance %d: %s (resource %d)\n") % _instances[i].handle % _instances[i].resourcename % _instances[i].resourcenumber);
        }
        return s;
    }
    
    inline Data::SharedPtr Likelihood::get_data() {
        return this->_data;
    }
    
    inline void Likelihood::set_data(Data::SharedPtr data) {
        assert(this->_instances.size() == 0);
        assert(! data->get_data_matrix().empty());
        this->_data = data;
    }

    inline Model::SharedPtr Likelihood::get_model() {  
        return _model;
    }
    
    inline void Likelihood::set_model(Model::SharedPtr m) {
        assert(_instances.size() == 0); // can't change model after initBeagleLib called
        _model = m;
    }

    inline void Likelihood::set_ambiguity_equals_missing(bool ambig_equals_missing) {
        // Can't change GPU preference status after init_beagle_lib called
        assert(this->_instances.size() == 0 || this->_ambiguity_equals_missing == ambig_equals_missing);
        this->_ambiguity_equals_missing = ambig_equals_missing;
    }
    
    inline void Likelihood::set_prefer_gpu(bool prefer_gpu) {
        // Can't change GPU preference status after init_beagle_lib called
        assert(this->_instances.size() == 0 || this->_prefer_gpu == prefer_gpu);
        this->_prefer_gpu = prefer_gpu;
    }
    
    inline bool Likelihood::using_stored_data() const {
        return this->_using_data;
    }
    
    inline void Likelihood::use_stored_data(bool using_data) {
        this->_using_data = using_data;
    }

    inline void Likelihood::use_underflow_scaling(bool do_scaling) {
        this->_underflow_scaling = do_scaling;
    }

    inline bool Likelihood::using_underflow_scaling() const {
        return this->_underflow_scaling;
    }
    
    inline void Likelihood::init_beagle_lib() {
        assert(this->_data);
        assert(this->_model);

        // Close down any existing BeagleLib instances
        this->finalize_beagle_lib(true);

        this->_ntaxa = this->_data->get_num_taxa();
        
        unsigned nsubsets = this->_data->get_num_subsets();
        std::set<instance_pair_t> nstates_ncateg_combinations;
        std::map<instance_pair_t, std::vector<unsigned> > subsets_for_pair;
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            // Create a pair comprising number of states and number of rate categories
            unsigned nstates = this->_data->get_num_states_for_subset(subset);
            bool invar_model = _model->get_subset_is_invar_model(subset);
            int nrates = (invar_model ? -1 : 1)*_model->get_subset_num_categ(subset);
            instance_pair_t p = std::make_pair(nstates, nrates);
            
            // Add combo to set
            nstates_ncateg_combinations.insert(p);
            subsets_for_pair[p].push_back(subset);
        }

        // Create one instance for each distinct nstates-nrates combination
        this->_instances.clear();
        for (auto p : nstates_ncateg_combinations) {
            this->new_instance(p.first, p.second, subsets_for_pair[p]);
            
            InstanceInfo & info = *this->_instances.rbegin();
            std::cout << boost::str(boost::format("Created BeagleLib instance %d (%d states, %d rate%s, %d subset%s, %s invar. sites model)") % info.handle % info.nstates % info.nratecateg % (info.nratecateg == 1 ? "" : "s") % info.subsets.size() % (info.subsets.size() == 1 ? "" : "s") % (info.invarmodel ? "is" : "not")) << std::endl;
        }
        
        if (this->_ambiguity_equals_missing)
            this->set_tip_states();
        else
            this->set_tip_partials();
        this->set_pattern_weights();
        this->set_pattern_partition_assignments();
    }
    
    inline void Likelihood::new_instance(unsigned nstates, int nrates, std::vector<unsigned> & subset_indices) {
        unsigned num_subsets = (unsigned)subset_indices.size();
        
        bool is_invar_model = (nrates < 0 ? true : false);
        unsigned ngammacat = (unsigned)(is_invar_model ? -nrates : nrates);

        // Create an identity matrix used for computing partials
        // for polytomies (represents the transition matrix
        // for the zero-length edges inserted to arbitrarily
        // resolve each polytomy)
        // _identity_matrix.assign(nstates*nstates*ngammacat, 0.0);
        // for (unsigned k = 0; k < ngammacat; k++) {
        //     unsigned offset = k*nstates*nstates;
        //     _identity_matrix[0+offset]  = 1.0;
        //     _identity_matrix[5+offset]  = 1.0;
        //     _identity_matrix[10+offset] = 1.0;
        //     _identity_matrix[15+offset] = 1.0;
        // }
        
        unsigned num_patterns = 0;
        for (auto s : subset_indices) {
            num_patterns += this->_data->get_num_patterns_in_subset(s);
        }
        
        unsigned num_internals = this->calc_num_internals_in_fully_resolved_tree();

        unsigned num_edges = this->calc_num_edges_in_fully_resolved_tree();
        // add 1 to num_edges so that subroot node will have a tmatrix, root tip's tmatrix is never used
        unsigned num_nodes = num_edges;
        // unsigned num_nodes = num_edges;
        unsigned num_transition_probs = num_nodes*num_subsets;
        
        long requirement_flags = 0;

        long preference_flags = BEAGLE_FLAG_PRECISION_SINGLE | BEAGLE_FLAG_THREADING_CPP;
        if (this->_underflow_scaling) {
            preference_flags |= BEAGLE_FLAG_SCALING_MANUAL;
            preference_flags |= BEAGLE_FLAG_SCALERS_LOG;
        }
        if (this->_prefer_gpu)
            preference_flags |= BEAGLE_FLAG_PROCESSOR_GPU;
        else
            preference_flags |= BEAGLE_FLAG_PROCESSOR_CPU;
        
        BeagleInstanceDetails instance_details;
        unsigned npartials = num_internals + this->_ntaxa;
        unsigned nscalers = num_internals;
        unsigned nsequences = 0;
        if (this->_ambiguity_equals_missing) {
            npartials -= this->_ntaxa;
            nsequences += this->_ntaxa;
        }
        
        int inst = beagleCreateInstance(
             _ntaxa,                           // tips
             npartials,                        // partials
             nsequences,                       // sequences
             nstates,                          // states
             num_patterns,                     // patterns (total across all subsets that use this instance)
             num_subsets,                      // models (one for each distinct eigen decomposition)
             num_subsets*num_transition_probs, // transition matrices (one for each node in each subset)
             ngammacat,                        // rate categories
             (_underflow_scaling ? nscalers + 1 : 0), // scale buffers (+1 is for the cummulative scaler at index 0)
             NULL,                             // resource restrictions
             0,                                // length of resource list
             preference_flags,                  // preferred flags
             requirement_flags,                 // required flags
             &instance_details);               // pointer for details
        
        if (inst < 0) {
            // beagleCreateInstance returns one of the following:
            //   valid instance (0, 1, 2, ...)
            //   error code (negative integer)
            throw NamuX(boost::str(boost::format("Likelihood init function failed to create BeagleLib instance (BeagleLib error code was %d)") % this->_beagle_error[inst]));
        }
        
        InstanceInfo info;
        info.handle         = inst;
        info.resourcenumber = instance_details.resourceNumber;
        info.resourcename   = instance_details.resourceName;
        info.nstates        = nstates;
        info.nratecateg     = ngammacat;
        info.invarmodel     = is_invar_model;
        info.subsets        = subset_indices;
        info.npatterns      = num_patterns;
        info.partial_offset = num_internals;
        info.t_matrix_offset = num_nodes;
        this->_instances.push_back(info);
    }

    inline void Likelihood::set_tip_states() {
        assert(this->_instances.size() > 0);
        assert(this->_data);
        Data::state_t one = 1;

        for (auto & info : this->_instances) {
            std::vector<int> states(info.nstates*info.npatterns);
            
            // Loop through all rows of the data matrix, setting the tip states for one taxon each row
            unsigned t = 0;
            for (auto & row : this->_data->get_data_matrix()) {
            
                // Loop through all subsets assigned to this instance
                unsigned k = 0;
                for (unsigned s : info.subsets) {
                
                    // Loop through all patterns in this subset
                    auto interval = this->_data->get_subset_begin_end(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                    
                        // d is the state for taxon t, pattern p (in subset s)
                        // d is stored as a bit field (e.g., for nucleotide data, A=1, C=2, G=4, T=8, ?=15),
                        // but BeagleLib expects states to be integers (e.g. for nucleotide data,
                        // A=0, C=1, G=2, T=3, ?=4).
                        Data::state_t d = row[p];
                        
                        // Handle common nucleotide case separately
                        if (info.nstates == 4) {
                            if (d == 1)
                                states[k++] = 0;
                            else if (d == 2)
                                states[k++] = 1;
                            else if (d == 4)
                                states[k++] = 2;
                            else if (d == 8)
                                states[k++] = 3;
                            else
                                states[k++] = 4;
                        }
                        else {
                            // This case is for any other data type except nucleotide
                            int s = -1;
                            for (unsigned b = 0; b < info.nstates; b++) {
                                if (d == one << b) {
                                    s = b;
                                    break;
                                }
                            }
                            if (s == -1)
                                states[k++] = info.nstates;
                            else
                                states[k++] = s;
                        }
                    } // pattern loop
                }   // subset loop

            int code = beagleSetTipStates(
                info.handle,    // Instance number
                t,              // Index of destination compactBuffer
                &states[0]);    // Pointer to compact states vector

            if (code != 0)
                throw NamuX(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t+1) % this->_data->get_taxon_names()[t] % code % this->_beagle_error[code]);
            ++t;
            }
        }
    }

    inline void Likelihood::set_tip_partials() {
        assert(this->_instances.size() > 0);
        assert(this->_data);
        Data::state_t one = 1;
        
        for (auto & info : this->_instances) {
            std::vector<double> partials(info.nstates*info.npatterns);
            
            // Loop through all rows of data matrix, setting the tip states for one taxon each row
            unsigned t = 0;
            for (auto & row : this->_data->get_data_matrix()) {
            
                // Loop through all subsets assigned to this instance
                unsigned k = 0;
                for (unsigned s : info.subsets) {
                
                    // Loop through all patterns in this subset
                    auto interval = this->_data->get_subset_begin_end(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                    
                        // d is the state for taxon t, pattern p (in subset s)
                        Data::state_t d = row[p];
                        
                        // Handle common nucleotide case separately
                        if (info.nstates == 4) {
                            partials[k++] = d & 1 ? 1.0 : 0.0;
                            partials[k++] = d & 2 ? 1.0 : 0.0;
                            partials[k++] = d & 4 ? 1.0 : 0.0;
                            partials[k++] = d & 8 ? 1.0 : 0.0;
                        }
                        else {
                            // This case is for any other data type except nucleotide
                            for (unsigned b = 0; b < info.nstates; b++) {
                                partials[k++] = d & (one << b) ? 1.0 : 0.0;
                            }
                        }
                    }
                }
                
            int code = beagleSetTipPartials(
                info.handle,    // Instance number
                t,              // Index of destination compactBuffer
                &partials[0]);  // Pointer to compact states vector

            if (code != 0)
                throw NamuX(boost::format("failed to set tip state for taxon %d (\"%s\"; BeagleLib error code was %d)") % (t+1) % this->_data->get_taxon_names()[t] % code % this->_beagle_error[code]);
            ++t;
            }
        }
    }

    inline void Likelihood::set_pattern_partition_assignments() {
        assert(this->_instances.size() > 0);
        assert(this->_data);
        
        // beagleSetPatternPartitions does not need to be called if data are unpartitioned
        // (and, in fact, BeagleLib only supports partitioning for 4-state instances if GPU is used,
        // so not calling beagleSetPatternPartitions allows unpartitioned codon model analyses)
        if (this->_instances.size() == 1 && this->_instances[0].subsets.size() == 1)
            return;
        
        Data::partition_key_t v;

        // Loop through all instances
        for (auto & info : this->_instances) {
            unsigned nsubsets = (unsigned)info.subsets.size();
            v.resize(info.npatterns);
            unsigned pattern_index = 0;

            // Loop through all subsets assigned to this instance
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                // Loop through all patterns in this subset
                auto interval = this->_data->get_subset_begin_end(s);
                for (unsigned p = interval.first; p < interval.second; p++) {
                    v[pattern_index++] = instance_specific_subset_index;
                }
                ++instance_specific_subset_index;
            }

            int code = beagleSetPatternPartitions(
               info.handle, // instance number
               nsubsets,    // number of data subsets (equals 1 if data are unpartitioned)
               &v[0]);      // vector of subset indices: v[i] = 0 means pattern i is in subset 0

            if (code != 0) {
                throw NamuX(boost::format("failed to set pattern partition. BeagleLib error code was %d (%s)") % code % this->_beagle_error[code]);
            }
        }
    }
    
    inline void Likelihood::set_pattern_weights() {
        assert(this->_instances.size() > 0);
        assert(this->_data);
        Data::pattern_counts_t v;
        auto pattern_counts = this->_data->get_pattern_counts();
        assert(pattern_counts.size() > 0);

        // Loop through all instances
        for (auto & info : this->_instances) {
            v.resize(info.npatterns);
            unsigned pattern_index = 0;

            // Loop through all subsets assigned to this instance
            for (unsigned s : info.subsets) {
            
                // Loop through all patterns in this subset
                auto interval = this->_data->get_subset_begin_end(s);
                for (unsigned p = interval.first; p < interval.second; p++) {
                    v[pattern_index++] = pattern_counts[p];
                }
            }

            int code = beagleSetPatternWeights(
               info.handle,   // instance number
               &v[0]);        // vector of pattern counts: v[i] = 123 means pattern i was encountered 123 times

            if (code != 0)
                throw NamuX(boost::format("Failed to set pattern weights for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % this->_beagle_error[code]);
        }
    }

    inline void Likelihood::set_among_site_rate_heterogeneity() {  
        assert(_instances.size() > 0);
        int code = 0;
        
        // Loop through all instances
        for (auto & info : _instances) {

            // Loop through all subsets assigned to this instance
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                code = _model->set_beagle_among_site_rate_variation_rates(info.handle, s, instance_specific_subset_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("Failed to set category rates for BeagleLib instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]));
            
                code = _model->set_beagle_among_site_rate_variation_probs(info.handle, s, instance_specific_subset_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("Failed to set category probabilities for BeagleLib instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]));
                    
                ++instance_specific_subset_index;
            }
        }
    }

    inline void Likelihood::set_model_rate_matrix() { 
        // Loop through all instances
        for (auto & info : _instances) {

            // Loop through all subsets assigned to this instance
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                int code = _model->set_beagle_state_frequencies(info.handle, s, instance_specific_subset_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("Failed to set state frequencies for BeagleLib instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]));

                code = _model->set_beagle_eigen_decomposition(info.handle, s, instance_specific_subset_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("Failed to set eigen decomposition for BeagleLib instance %d. BeagleLib error code was %d (%s)") % info.handle % code % _beagle_error[code]));
                
                ++instance_specific_subset_index;
            }
        }
    }
    
    inline unsigned Likelihood::get_scaler_index(Node * nd, InstanceInfo & info) const {
        unsigned sindex = BEAGLE_OP_NONE;
        if (this->_underflow_scaling) {
            sindex = nd->_number - _ntaxa + 1; // +1 to skip the cumulative scaler vector
            assert((sindex >= 1) && (sindex <= this->calc_num_internals_in_fully_resolved_tree()));
        }
        return sindex;
    }
    
    inline unsigned Likelihood::get_partial_index(Node * nd, InstanceInfo & info) const {
        // Note: do not be tempted to subtract _ntaxa from pindex: BeagleLib does this itself
        assert(nd->_number >= 0);
        return nd->_number;
    }
    
    inline unsigned Likelihood::get_t_matrix_index(Node * nd, InstanceInfo & info, unsigned subset_index) const {
        assert(! nd->is_root());
        unsigned tindex = subset_index*info.t_matrix_offset + nd->_number;
        return tindex;
    }

    inline unsigned Likelihood::get_t_matrix_index_from_node_num(int node_num, InstanceInfo & info, unsigned subset_index) const {
        unsigned tindex = subset_index*info.t_matrix_offset + node_num;
        return tindex;
    }
      
    inline void Likelihood::define_operations(Tree::SharedPtr t) {
        assert(this->_instances.size() > 0);
        assert(t);
        
        this->_relrate_normalizing_constant = this->_model->calc_normalizing_constant_for_subset_rel_rates();

        // Start with a clean slate
        for (auto & info : this->_instances) {
            this->_operations[info.handle].clear();
            this->_pmatrix_index[info.handle].clear();
            this->_edge_lengths[info.handle].clear();
            this->_eigen_indices[info.handle].clear();
            this->_category_rate_indices[info.handle].clear();
        }
                
        // Loop through all nodes in reverse level order
        for (auto nd : boost::adaptors::reverse(t->_levelorder)) {
            assert(nd->_number >= 0);
            if (!nd->_left_child) {
                // This is a leaf
                this->queue_t_matrix_recalculation(nd);
            }
            else {
                // This is an internal node
                if (! nd->is_root()) {
                    this->queue_t_matrix_recalculation(nd);
                }

                // Internal nodes have partials to be calculated, so define
                // an operation to compute the partials for this node
                // TODO: Assuming no polytomies here
                Node * lchild = nd->_left_child;
                assert(lchild);
                Node * rchild = lchild->_right_sib;
                assert(rchild);
                this->queue_partials_recalculation(nd, lchild, rchild);
            }
        }
    }
    
    inline void Likelihood::queue_partials_recalculation(Node * nd, Node * lchild, Node * rchild) { // , Node * polytomy) {
        for (auto & info : this->_instances) {
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                // if (polytomy) {
                //     // nd has been pulled out of tree's _unused_nodes vector to break up the polytomy
                //     // Note that the parameter "polytomy" is the polytomous node itself
                //     
                //     // First get the transition matrix index
                //     unsigned tindex = this->get_t_matrix_index(nd, info, instance_specific_subset_index);

                //     // Set the transition matrix for nd to the identity matrix
                //     // note: last argument 1 is the value used for ambiguous states (should be 1 for transition matrices)
                //     int code = beagleSetTransitionMatrix(info.handle, tindex, &_identity_matrix[0], 1);
                //     if (code != 0)
                //         throw XStrom(boost::str(boost::format("Failed to set transition matrix for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % this->_beagle_error[code]));
                //     
                //     // Set the edgelength to 0.0 to maintain consistency with the transition matrix
                //     // nd->setEdgeLength(0.0);
                //     // This unused node should not have parent and thus should
                //     // return an edge length of zero
                //     assert(nd->get_edge_length() == 0.0);

                //     // If employing underflow scaling, the scaling factors for these fake nodes need to be
                //     // transferred to the polytomous node, as that will be the only node remaining after the
                //     // likelihood has been calculated. Save the scaling factor index, associating it with
                //     // the scaling factor index of the polytomy node.
                //     if (this->_underflow_scaling) {
                //         // Get the polytomy's scaling factor index
                //         int spolytomy = this->get_scaler_index(polytomy, info);
                //         
                //         // Get nd's scaling factor index
                //         int snd = this->get_scaler_index(nd, info);

                //         // Save nd's index in the vector associated with polytomy's index
                //         this->_polytomy_map[spolytomy].push_back(snd);
                //     }
                // }

                this->add_operation(info, nd, lchild, rchild, instance_specific_subset_index);
                ++instance_specific_subset_index;
            }
        }
    }
    
    inline void Likelihood::queue_t_matrix_recalculation(Node * nd) {
        Model::subset_relrate_vect_t & subset_relrates = _model->get_subset_rel_rates();
        for (auto & info : this->_instances) {
            unsigned instance_specific_subset_index = 0;
            for (unsigned s : info.subsets) {
                double subset_relative_rate = subset_relrates[s]/_relrate_normalizing_constant;

                unsigned tindex = this->get_t_matrix_index(nd, info, instance_specific_subset_index);
                this->_pmatrix_index[info.handle].push_back(tindex);
                this->_edge_lengths[info.handle].push_back(nd->get_edge_length()*subset_relative_rate);
                this->_eigen_indices[info.handle].push_back(s);
                this->_category_rate_indices[info.handle].push_back(s);

                ++instance_specific_subset_index;
            }
        }
    }
    
    inline void Likelihood::add_operation(InstanceInfo & info, Node * nd, Node * lchild, Node * rchild, unsigned subset_index) {
        assert(nd);
        assert(lchild);
        assert(rchild);

        // 1. destination partial to be calculated
        int partial_dest = this->get_partial_index(nd, info);
        this->_operations[info.handle].push_back(partial_dest);

        // 2. destination scaling buffer index to write to
        int scaler_index = this->get_scaler_index(nd, info);
        this->_operations[info.handle].push_back(scaler_index);

        // 3. destination scaling buffer index to read from
        this->_operations[info.handle].push_back(BEAGLE_OP_NONE);

        // 4. left child partial index
        int partial_lchild = this->get_partial_index(lchild, info);
        this->_operations[info.handle].push_back(partial_lchild);

        // 5. left child transition matrix index
        unsigned tindex_lchild = this->get_t_matrix_index(lchild, info, subset_index);
        this->_operations[info.handle].push_back(tindex_lchild);

        // 6. right child partial index
        int partial_rchild = this->get_partial_index(rchild, info);
        this->_operations[info.handle].push_back(partial_rchild);

        // 7. right child transition matrix index
        unsigned tindex_rchild = this->get_t_matrix_index(rchild, info, subset_index);
        this->_operations[info.handle].push_back(tindex_rchild);

        if (info.subsets.size() > 1) {
            // 8. index of partition subset
            this->_operations[info.handle].push_back(subset_index);
            
            // 9. cumulative scale index
            this->_operations[info.handle].push_back(BEAGLE_OP_NONE);
        }
    }
    
    inline void Likelihood::update_transition_matrices() {
        assert(this->_instances.size() > 0);
        if (this->_pmatrix_index.size() == 0)
            return;

        // Loop through all instances
        for (auto & info : this->_instances) {
            int code = 0;

            unsigned nsubsets = (unsigned)info.subsets.size();
            if (nsubsets > 1) {
                code = beagleUpdateTransitionMatricesWithMultipleModels(
                    info.handle,                                // Instance number
                    &_eigen_indices[info.handle][0],            // Index of eigen-decomposition buffer
                    &_category_rate_indices[info.handle][0],    // category rate indices
                    &_pmatrix_index[info.handle][0],            // transition probability matrices to update
                    NULL,                                       // first derivative matrices to update
                    NULL,                                       // second derivative matrices to update
                    &_edge_lengths[info.handle][0],             // List of edge lengths
                    (int)_pmatrix_index[info.handle].size());   // Length of lists
            }
            else {
                code = beagleUpdateTransitionMatrices(
                    info.handle,                                // Instance number
                    0,                                          // Index of eigen-decomposition buffer
                    &_pmatrix_index[info.handle][0],            // transition probability matrices to update
                    NULL,                                       // first derivative matrices to update
                    NULL,                                       // second derivative matrices to update
                    &_edge_lengths[info.handle][0],             // List of edge lengths
                    (int)_pmatrix_index[info.handle].size());   // Length of lists
            }

            if (code != 0)
                throw NamuX(boost::str(boost::format("Failed to update transition matrices for instance %d. BeagleLib error code was %d (%s)") % info.handle % code % this->_beagle_error[code]));
        } 
    }
    
    inline void Likelihood::calculate_partials() {
        assert(this->_instances.size() > 0);
        if (this->_operations.size() == 0)
            return;
        int code = 0;
        
        // Loop through all instances
        for (auto & info : this->_instances) {
            unsigned nsubsets = (unsigned)info.subsets.size();

            if (nsubsets > 1) {
                code = beagleUpdatePartialsByPartition(
                    info.handle,                                                    // Instance number
                    (BeagleOperationByPartition *) &_operations[info.handle][0],    // BeagleOperation list specifying operations
                    (int)(_operations[info.handle].size()/9));                      // Number of operations
                if (code != 0)
                    throw NamuX(boost::format("failed to update partials. BeagleLib error code was %d (%s)") % code % this->_beagle_error[code]);
            }
            else {
                // no partitioning, just one data subset
                code = beagleUpdatePartials(
                    info.handle,                                        // Instance number
                    (BeagleOperation *) &_operations[info.handle][0],   // BeagleOperation list specifying operations
                    (int)(_operations[info.handle].size()/7),           // Number of operations
                    BEAGLE_OP_NONE);                                    // Index number of scaleBuffer to store accumulated factors
                if (code != 0) 
                    throw NamuX(boost::format("failed to update partials. BeagleLib error code was %d (%s)") % code % this->_beagle_error[code]);
            }
        } 
    }
    
    inline double Likelihood::calc_instance_log_likelihood(InstanceInfo & info, Tree::SharedPtr t) {
        int code = 0;
        unsigned nsubsets = (unsigned)info.subsets.size();
        assert(nsubsets > 0);
        
        // Assuming there are as many transition matrices as there are edge lengths
        assert(this->_pmatrix_index[info.handle].size() == this->_edge_lengths[info.handle].size());

        int state_frequency_index  = 0;
        int category_weights_index = 0;
        int cumulative_scale_index = (this->_underflow_scaling ? 0 : BEAGLE_OP_NONE);
        int child_partials_index   = this->get_partial_index(t->_root, info);
        int parent_partials_index  = this->get_partial_index(t->_preorder[0], info);
        assert(child_partials_index == parent_partials_index);
        // int parent_t_matrix_index   = this->get_t_matrix_index(t->_preorder[0], info, 0);

        // storage for results of the likelihood calculation
        std::vector<double> subset_log_likelihoods(nsubsets, 0.0);
        double log_likelihood = 0.0;

        if (this->_underflow_scaling) {
            // Create vector of all scaling vector indices in current use
            std::vector<int> internal_node_scaler_indices;
            for (auto nd : t->_preorder) {
                if (nd->_left_child) {
                    unsigned s = get_scaler_index(nd, info);
                    internal_node_scaler_indices.push_back(s);
                }
            }

            if (nsubsets == 1) {
                code = beagleResetScaleFactors(info.handle, cumulative_scale_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("failed to reset scale factors in calcInstanceLogLikelihood. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));


                code = beagleAccumulateScaleFactors(
                     info.handle,
                     &internal_node_scaler_indices[0],
                     (int)internal_node_scaler_indices.size(),
                     cumulative_scale_index);
                if (code != 0)
                    throw NamuX(boost::str(boost::format("failed to accumulate scale factors in calcInstanceLogLikelihood. BeagleLib error code was %d (%s)") % code % _beagle_error[code]));
            }
            else {
                for (unsigned s = 0; s < nsubsets; ++s) {
                    code = beagleResetScaleFactorsByPartition(info.handle, cumulative_scale_index, s);
                    if (code != 0)
                        throw NamuX(boost::str(boost::format("failed to reset scale factors for subset %d in calcInstanceLogLikelihood. BeagleLib error code was %d (%s)") % s % code % _beagle_error[code]));
                        

                    code = beagleAccumulateScaleFactorsByPartition(
                        info.handle,
                        &internal_node_scaler_indices[0],
                        (int)internal_node_scaler_indices.size(),
                        cumulative_scale_index,
                        s);
                    if (code != 0)
                        throw NamuX(boost::str(boost::format("failed to acccumulate scale factors for subset %d in calcInstanceLogLikelihood. BeagleLib error code was %d (%s)") % s % code % _beagle_error[code]));
                }
            }
        }

        if (nsubsets > 1) {
            this->_parent_indices.assign(nsubsets, parent_partials_index);
            this->_child_indices.assign(nsubsets, child_partials_index);
            this->_weights_indices.assign(nsubsets, category_weights_index);
            this->_scaling_indices.resize(nsubsets); 
            this->_subset_indices.resize(nsubsets);
            this->_freqs_indices.resize(nsubsets);
            this->_t_matrix_indices.resize(nsubsets);

            for (unsigned s = 0; s < nsubsets; s++) {
                this->_scaling_indices[s] = (this->_underflow_scaling ? 0 : BEAGLE_OP_NONE); 
                this->_subset_indices[s]  = s;
                this->_freqs_indices[s]   = s;
                this->_t_matrix_indices[s] = this->get_t_matrix_index_from_node_num(t->_preorder[0]->_number - 1, info, s); //index_focal_child + s*t_matrix_skip;
            }
            // code = beagleCalculateEdgeLogLikelihoodsByPartition(
            //     info.handle,                 // instance number
            //     &_parent_indices[0],         // indices of parent partialsBuffers
            //     &_child_indices[0],          // indices of child partialsBuffers
            //     &_t_matrix_indices[0],        // transition probability matrices for this edge
            //     NULL,                        // first derivative matrices
            //     NULL,                        // second derivative matrices
            //     &_weights_indices[0],        // weights to apply to each partialsBuffer
            //     &_freqs_indices[0],          // state frequencies for each partialsBuffer
            //     &_scaling_indices[0],        // scaleBuffers containing accumulated factors
            //     &_subset_indices[0],         // indices of subsets
            //     nsubsets,                    // partition subset count
            //     1,                           // number of distinct eigen decompositions
            //     &subset_log_likelihoods[0],  // address of vector of log likelihoods (one for each subset)
            //     &log_likelihood,             // destination for resulting log likelihood
            //     NULL,                        // destination for vector of first derivatives (one for each subset)
            //     NULL,                        // destination for first derivative
            //     NULL,                        // destination for vector of second derivatives (one for each subset)
            //     NULL);                       // destination for second derivative
            code = beagleCalculateRootLogLikelihoodsByPartition(
                info.handle,                 // instance number
                &_parent_indices[0],         // indices of parent partialsBuffers
                &_weights_indices[0],        // weights to apply to each partialsBuffer
                &_freqs_indices[0],          // state frequencies for each partialsBuffer
                &_scaling_indices[0],        // scaleBuffers containing accumulated factors
                &_subset_indices[0],         // indices of subsets
                nsubsets,                    // partition subset count
                1,                           // number of distinct eigen decompositions
                &subset_log_likelihoods[0],  // address of vector of log likelihoods (one for each subset)
                &log_likelihood);             // destination for resulting log likelihood
        }
        else {
            // code = beagleCalculateEdgeLogLikelihoods(
            //     info.handle,                 // instance number
            //     &parent_partials_index,      // indices of parent partialsBuffers
            //     &child_partials_index,       // indices of child partialsBuffers
            //     &parent_t_matrix_index,       // transition probability matrices for this edge
            //     NULL,                        // first derivative matrices
            //     NULL,                        // second derivative matrices
            //     &category_weights_index,     // weights to apply to each partialsBuffer
            //     &state_frequency_index,      // state frequencies for each partialsBuffer
            //     &cumulative_scale_index,     // scaleBuffers containing accumulated factors
            //     1,                           // Number of partialsBuffer
            //     &log_likelihood,             // destination for log likelihood
            //     NULL,                        // destination for first derivative
            //     NULL);                       // destination for second derivative
            code = beagleCalculateRootLogLikelihoods(
                info.handle,                 // instance number
                &parent_partials_index,      // indices of parent partialsBuffers
                &category_weights_index,     // weights to apply to each partialsBuffer
                &state_frequency_index,      // state frequencies for each partialsBuffer
                &cumulative_scale_index,     // scaleBuffers containing accumulated factors
                1,                           // Number of partialsBuffer
                &log_likelihood);             // destination for log likelihood
        }

        if (code != 0)
            throw NamuX(boost::str(boost::format("failed to calculate edge logLikelihoods in CalcLogLikelihood. BeagleLib error code was %d (%s)") % code % this->_beagle_error[code]));

        if (info.invarmodel) {
            auto monomorphic = _data->get_monomorphic();
            auto counts = _data->get_pattern_counts();
            std::vector<double> site_log_likelihoods(info.npatterns, 0.0);
            double * siteLogLs = &site_log_likelihoods[0];


            beagleGetSiteLogLikelihoods(info.handle, siteLogLs);


            // Loop through all subsets assigned to this instance
            double lnL = 0.0;
            unsigned i = 0;
            for (unsigned s : info.subsets) {
                const ASRV & asrv = _model->get_asrv(s);
                const QMatrix & qmatrix = _model->get_q_matrix(s);
                const double * freq = qmatrix.get_state_freqs();
                

                double pinvar = *(asrv.get_pinvar_shared_ptr());
                assert(pinvar >= 0.0 && pinvar <= 1.0);


                if (pinvar == 0.0) {
                    // log likelihood for this subset is equal to the sum of site log-likelihoods
                    auto interval = _data->get_subset_begin_end(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                        lnL += counts[p]*site_log_likelihoods[i++];
                    }
                }
                else {
                    // Loop through all patterns in this subset
                    double log_pinvar = log(pinvar);
                    double log_one_minus_pinvar = log(1.0 - pinvar);
                    auto interval = _data->get_subset_begin_end(s);
                    for (unsigned p = interval.first; p < interval.second; p++) {
                        // Loop through all states for this pattern
                        double invar_like = 0.0;
                        if (monomorphic[p] > 0) {
                            for (unsigned k = 0; k < info.nstates; ++k) {
                                Data::state_t x = (Data::state_t)1 << k;
                                double condlike = (x & monomorphic[p] ? 1.0 : 0.0);
                                double basefreq = freq[k];
                                invar_like += condlike*basefreq;
                            }
                        }
                        double site_lnL = site_log_likelihoods[i++];
                        double log_like_term = log_one_minus_pinvar + site_lnL;
                        if (invar_like > 0.0) {
                            double log_invar_term = log_pinvar + log(invar_like);
                            double site_log_like = (log_like_term + log(1.0 + exp(log_invar_term - log_like_term)));
                            lnL += counts[p]*site_log_like;
                        }
                        else {
                            lnL += counts[p]*log_like_term;
                        }
                    }
                }
            }
            log_likelihood = lnL;
        }

        return log_likelihood;
    }
    
    inline double Likelihood::calc_log_likelihood(Tree::SharedPtr t) {
        assert(this->_instances.size() > 0);
        
        if (! this->_using_data)
            return 0.0;

        // Must call set_data and set_model before calc_log_likelihood
        assert(this->_data);
        assert(this->_model);

        // if (t->_is_rooted)
        //     throw NamuX("This version of the program can only compute likelihoods for unrooted trees");

        // Assuming "root" is leaf 0
        // assert(t->_root->_number == 0 && t->_root->_left_child == t->_preorder[0] && !t->_preorder[0]->_right_sib);
        assert(
            (t->_root->_number == (int)(t->_preorder.size() - 1))
            && (t->_root == t->_preorder[0])
            && (t->_root->_left_child == t->_preorder[1])
        );

        this->set_model_rate_matrix();
        this->set_among_site_rate_heterogeneity();
        this->define_operations(t);
        this->update_transition_matrices();
        this->calculate_partials();

        double log_likelihood = 0.0;
        for (auto & info : this->_instances) {
            log_likelihood += this->calc_instance_log_likelihood(info, t);
        }
        
        return log_likelihood;
    }
}
