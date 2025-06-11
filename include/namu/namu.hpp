#pragma once

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include "namu/data.hpp"
#include "namu/likelihood.hpp"
#include "namu/tree_summary.hpp"
#include "namu/partition.hpp"

namespace namu {

    class Namu {
        public:
                                        Namu();
                                        ~Namu();

            void                        clear();
            void                        process_command_line_options(int argc, const char * argv[]);
            void                        run();

        private:
            bool                        process_assignment_string(
                                            const std::string & which,
                                            const std::string & definition);
            void                        handle_assignment_strings(
                                            const boost::program_options::variables_map & vm,
                                            std::string label,
                                            const std::vector<std::string> & definitions,
                                            std::string default_definition);
            bool                        split_assignment_string(
                                            const std::string & definition,
                                            std::vector<std::string> & vector_of_subset_names,
                                            std::vector<double>  & vector_of_values);

            double                      _expected_log_likelihood;

            std::string                 _conf_file_path;
            std::string                 _data_file_path;
            std::string                 _tree_file_path;

            Partition::SharedPtr        _partition;
            Data::SharedPtr             _data;
            Model::SharedPtr            _model;
            Likelihood::SharedPtr       _likelihood;

            TreeSummary::SharedPtr      _tree_summary;

            bool                        _use_gpu;
            bool                        _ambig_missing;
            bool                        _use_underflow_scaling;

            static std::string          _program_name;
            static unsigned             _major_version;
            static unsigned             _minor_version;

    };

    inline Namu::Namu() {
        //std::cout << "Constructing a Namu" << std::endl;
        clear();
    }

    inline Namu::~Namu() {
        //std::cout << "Destroying a Namu" << std::endl;
    }

    inline void Namu::clear() {
        this->_conf_file_path = "";
        this->_data_file_path = "";
        this->_tree_file_path = "";
        this->_tree_summary   = nullptr;
        this->_partition.reset(new Partition());
        this->_use_gpu        = true;
        this->_ambig_missing  = true;
        this->_model.reset(new Model());
        this->_expected_log_likelihood = 0.0;
        this->_data = nullptr; 
        this->_likelihood = nullptr;
        this->_use_underflow_scaling = false;
    }

    inline void Namu::process_command_line_options(int argc, const char * argv[]) {
        std::vector<std::string> partition_statefreq;
        std::vector<std::string> partition_rmatrix;
        std::vector<std::string> partition_omega;
        std::vector<std::string> partition_ratevar;
        std::vector<std::string> partition_pinvar;
        std::vector<std::string> partition_ncateg;
        std::vector<std::string> partition_subsets;
        std::vector<std::string> partition_relrates;
        std::vector<std::string> partition_tree;
        // Need temp string variables to store options (vs member variables),
        // because member variables were being set *after* this function and so
        // any modifications to them during the function were being overridden
        // and set back the the raw command line or config file values
        std::string conf_path;
        std::string data_path;
        std::string tree_path;
        boost::program_options::options_description conf_opt("Command line only options");
        conf_opt.add_options()
            ("config,c",  boost::program_options::value(&conf_path), "Path to a namu config file")
            ("help,h", "Produce help message")
            ("version,v", "Show program version")
        ;
        boost::program_options::options_description main_options("Command line and config file options");
        main_options.add_options()
            ("datafile,d",  boost::program_options::value(&data_path)->required(), "Path to a data file in NEXUS format")
            ("treefile,t",  boost::program_options::value(&tree_path)->required(), "Path to a tree file in NEXUS format")
            ("subset",  boost::program_options::value(&partition_subsets), "A string defining a partition subset, e.g. 'first:1-1234\3' or 'default[codon:standard]:1-3702'")
            ("ncateg", boost::program_options::value(&partition_ncateg), "Number of categories in the discrete Gamma rate heterogeneity model")
            ("statefreq", boost::program_options::value(&partition_statefreq), "A string defining state frequencies for one or more data subsets, e.g. 'first,second:0.1,0.2,0.3,0.4'")
            ("omega", boost::program_options::value(&partition_omega), "A string defining the nonsynonymous/synonymous rate ratio omega for one or more data subsets, e.g. 'first,second:0.1'")
            ("rmatrix", boost::program_options::value(&partition_rmatrix), "A string defining the rmatrix for one or more data subsets, e.g. 'first,second:1,2,1,1,2,1'")
            ("ratevar", boost::program_options::value(&partition_ratevar), "A string defining the among-site rate variance for one or more data subsets, e.g. 'first,second:2.5'")
            ("pinvar", boost::program_options::value(&partition_pinvar), "A string defining the proportion of invariable sites for one or more data subsets, e.g. 'first,second:0.2'")
            ("relrate", boost::program_options::value(&partition_relrates), "A string defining the (unnormalized) relative rates for all data subsets (e.g. 'default:3,1,6').")
            ("tree", boost::program_options::value(&partition_tree), "The index of the tree in the tree file (first tree has index = 1)")
            ("expectedLnL", boost::program_options::value(&_expected_log_likelihood)->default_value(0.0), "Log likelihood expected")
            ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true), "Use GPU if available")
            ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true), "Treat all ambiguities as missing data")
            ("underflowscaling", boost::program_options::value(&_use_underflow_scaling)->default_value(false), "Scale site likelihoods to prevent underflow (slower, but safer)")
        ;

        // Aggregate options allowed on command line
        boost::program_options::options_description cmd_line_options;
        cmd_line_options.add(conf_opt).add(main_options);

        // Aggregate options allowed in config file
        // This is redundant with `main_options` currently, but is explicit and
        // makes it easy to add other groups of options in the future
        boost::program_options::options_description conf_file_options;
        conf_file_options.add(main_options);

        boost::program_options::variables_map vm;

        // parse command line options
        boost::program_options::store(boost::program_options::parse_command_line(argc, argv, cmd_line_options), vm);

        // Update config file path if provided
        boost::filesystem::path config_dir;
        if (vm.count("config") > 0) {
            std::string input_config_str = vm.at("config").as<std::string>();
            boost::filesystem::path input_config_path(input_config_str);
            // get path relative to current working directory
            boost::filesystem::path abs_config_path = boost::filesystem::canonical(input_config_path);
            this->_conf_file_path = abs_config_path.string();
            config_dir = abs_config_path.parent_path();
        }
        // Update data file path if provided
        bool data_file_provided = false;
        if (vm.count("datafile") > 0) {
            data_file_provided = true;
            std::string input_data_file_str = vm.at("datafile").as<std::string>();
            boost::filesystem::path input_data_file_path(input_data_file_str);
            // get path relative to current working directory
            boost::filesystem::path abs_data_file_path = boost::filesystem::canonical(input_data_file_path);
            this->_data_file_path = abs_data_file_path.string();
        }
        // Update tree file path if provided
        bool tree_file_provided = false;
        if (vm.count("treefile") > 0) {
            tree_file_provided = true;
            std::string input_tree_file_str = vm.at("treefile").as<std::string>();
            boost::filesystem::path input_tree_file_path(input_tree_file_str);
            // get path relative to current working directory
            boost::filesystem::path abs_tree_file_path = boost::filesystem::canonical(input_tree_file_path);
            this->_tree_file_path = abs_tree_file_path.string();
        }
        if (vm.count("config") > 0) {
            // Parse options from config file
            std::ifstream config_stream(this->_conf_file_path);
            if (! config_stream.is_open()) {
                std::cerr << "Error opening config file" << std::endl;
                std::exit(1);
            }
            try {
                // parse config file options (all options except the config
                // file option)
                const boost::program_options::parsed_options & parsed = boost::program_options::parse_config_file< char >(config_stream, conf_file_options, false);
                boost::program_options::store(parsed, vm);
            }
            catch(const boost::program_options::error & e) {
                std::cerr << "Error parsing config file: " << e.what() << std::endl;
                throw;
            }
        }

        // If user specified --help on command line, output usage summary and quit
        if (vm.count("help") > 0) {
            std::cout << cmd_line_options << "\n";
            std::exit(1);
        }

        // If user specified --version on command line, output version and quit
        if (vm.count("version") > 0) {
            std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
            std::exit(1);
        }
        
        // No command line args provided
        if (argc == 1) {
            std::cout << boost::str(boost::format("This is %s version %d.%d") % _program_name % _major_version % _minor_version) << std::endl;
            std::cout << cmd_line_options << "\n";
            std::exit(1);
        }

        // Update data file path if provided in config file
        if ((! data_file_provided) && (vm.count("datafile") > 0)) {
            data_file_provided = true;
            std::string input_data_file_str = vm.at("datafile").as<std::string>();
            boost::filesystem::path input_data_file_path(input_data_file_str);
            // get path relative to config file
            boost::filesystem::path abs_data_file_path = boost::filesystem::canonical(input_data_file_path, config_dir);
            this->_data_file_path = abs_data_file_path.string();
        }
        // Update tree file path if provided in config file
        if ((! tree_file_provided) && (vm.count("treefile") > 0)) {
            tree_file_provided = true;
            std::string input_tree_file_str = vm.at("treefile").as<std::string>();
            boost::filesystem::path input_tree_file_path(input_tree_file_str);
            // get path relative to config file
            boost::filesystem::path abs_tree_file_path = boost::filesystem::canonical(input_tree_file_path, config_dir);
            this->_tree_file_path = abs_tree_file_path.string();
        }

        // Need to run this after help and version, otherwise missing required
        // options would throw an error prior to processing help or version
        // options
        boost::program_options::notify(vm);

        // If user specified --subset on command line, break specified partition subset
        // definition into name and character set string and add to _partition
        if (vm.count("subset") > 0) {
            this->_partition.reset(new Partition());
            for (auto s : partition_subsets) {
                this->_partition->parse_subset_definition(s);
            }
        }

        this->_model->set_subset_data_types(this->_partition->get_subset_data_types());

        this->handle_assignment_strings(vm, "statefreq", partition_statefreq, "default:equal");
        this->handle_assignment_strings(vm, "rmatrix",   partition_rmatrix,   "default:equal");
        this->handle_assignment_strings(vm, "omega",     partition_omega,     "default:0.1"  );
        this->handle_assignment_strings(vm, "ncateg",    partition_ncateg,    "default:1"    );
        this->handle_assignment_strings(vm, "ratevar",   partition_ratevar,   "default:1.0"  );
        this->handle_assignment_strings(vm, "pinvar",    partition_pinvar,    "default:0.0"  );
        this->handle_assignment_strings(vm, "relrate",   partition_relrates,  "default:equal");
        this->handle_assignment_strings(vm, "tree",      partition_tree,      "default:1"    );
    }

    inline void Namu::handle_assignment_strings(const boost::program_options::variables_map & vm, std::string label, const std::vector<std::string> & definitions, std::string default_definition) { 
        if (vm.count(label) > 0) {
            bool first = true;
            for (auto s : definitions) {
                bool is_default = this->process_assignment_string(label, s);
                if (is_default && !first)
                    throw NamuX(boost::format("default specification must be first %s encountered") % label);
                first = false;
            }
        }
        else {
            this->process_assignment_string(label, default_definition);
        }
    }

    inline bool Namu::process_assignment_string(const std::string & which, const std::string & definition) { 
        unsigned num_subsets_defined = this->_partition->get_num_subsets();
        std::vector<std::string> vector_of_subset_names;
        std::vector<double> vector_of_values;
        bool fixed = this->split_assignment_string(definition, vector_of_subset_names, vector_of_values);
        
        if (vector_of_values.size() == 1 && vector_of_values[0] == -1 && !(which == "statefreq" || which == "rmatrix" || which == "relrate"))
            throw NamuX("Keyword equal is only allowed for statefreq, rmatrix, and relrate");

        // Assign values to subsets in model
        bool default_found = false;
        if (which == "statefreq") {
            QMatrix::freq_xchg_ptr_t freqs = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    this->_model->set_subset_state_freqs(freqs, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    this->_model->set_subset_state_freqs(freqs, this->_partition->find_subset_by_name(s), fixed);
                }
            }
        }
        else if (which == "rmatrix") {
            QMatrix::freq_xchg_ptr_t xchg = std::make_shared<QMatrix::freq_xchg_t>(vector_of_values);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    this->_model->set_subset_exchangeabilities(xchg, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    this->_model->set_subset_exchangeabilities(xchg, this->_partition->find_subset_by_name(s), fixed);
                }
            }
        }
        else if (which == "omega") {
            if (vector_of_values.size() > 1)
                throw NamuX(boost::format("expecting 1 value for omega, found %d values") % vector_of_values.size());
            QMatrix::omega_ptr_t omega = std::make_shared<QMatrix::omega_t>(vector_of_values[0]);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    this->_model->set_subset_omega(omega, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    this->_model->set_subset_omega(omega, this->_partition->find_subset_by_name(s), fixed);
                }
            }
        }
        else if (which == "pinvar") {
            if (vector_of_values.size() > 1)
                throw NamuX(boost::format("expecting 1 value for pinvar, found %d values") % vector_of_values.size());
            ASRV::pinvar_ptr_t p = std::make_shared<double>(vector_of_values[0]);
            bool invar_model = (*p > 0);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++) {
                    this->_model->set_subset_is_invar_model(invar_model, i);
                    this->_model->set_subset_pinvar(p, i, fixed);
                }
            }
            else {
                for (auto s : vector_of_subset_names) {
                    unsigned i = this->_partition->find_subset_by_name(s);
                    this->_model->set_subset_is_invar_model(invar_model, i);
                    this->_model->set_subset_pinvar(p, i, fixed);
                }
            }
        }
        else if (which == "ratevar") {
            if (vector_of_values.size() > 1)
                throw NamuX(boost::format("expecting 1 value for ratevar, found %d values") % vector_of_values.size());
            ASRV::ratevar_ptr_t rv = std::make_shared<double>(vector_of_values[0]);
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    this->_model->set_subset_rate_var(rv, i, fixed);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    this->_model->set_subset_rate_var(rv, this->_partition->find_subset_by_name(s), fixed);
                }
            }
        }
        else if (which == "ncateg") {
            if (vector_of_values.size() > 1)
                throw NamuX(boost::format("expecting 1 value for ncateg, found %d values") % vector_of_values.size());
            unsigned ncat = vector_of_values[0];
            if (vector_of_subset_names[0] == "default") {
                default_found = true;
                for (unsigned i = 0; i < num_subsets_defined; i++)
                    this->_model->set_subset_num_categ(ncat, i);
            }
            else {
                for (auto s : vector_of_subset_names) {
                    this->_model->set_subset_num_categ(ncat, this->_partition->find_subset_by_name(s));
                }
            }
        }
        else if (which == "tree") {
            if (vector_of_values.size() > 1)
                throw NamuX(boost::format("expecting 1 value for tree, found %d values") % vector_of_values.size());
            unsigned tree_index = vector_of_values[0];
            assert(tree_index > 0);
            this->_model->set_tree_index(tree_index - 1, fixed);
            if (vector_of_subset_names[0] != "default")
                throw NamuX("tree must be assigned to default only");
        }
        else {
            assert(which == "relrate");
            if (vector_of_subset_names[0] != "default")
                throw NamuX("relrate must be assigned to default only");
            this->_model->set_subset_rel_rates(vector_of_values, fixed);
        }

        return default_found;
    }

    inline bool Namu::split_assignment_string(const std::string & definition, std::vector<std::string> & vector_of_subset_names, std::vector<double>  & vector_of_values) {  
        // Split subset names from definition
        std::vector<std::string> twoparts;
        boost::split(twoparts, definition, boost::is_any_of(":"));
        if (twoparts.size() != 2)
            throw NamuX("Expecting exactly one colon in assignment");
        std::string comma_delimited_subset_names_string = twoparts[0];
        std::string comma_delimited_value_string = twoparts[1];
        boost::to_lower(comma_delimited_value_string);
        
        // Check for brackets indicating that parameter should be fixed     
        // now see if before_colon contains a data type specification in square brackets
        bool fixed = false;
        const char * pattern_string = R"(\s*\[(.+?)\]\s*)";
        std::regex re(pattern_string);
        std::smatch match_obj;
        bool matched = std::regex_match(comma_delimited_value_string, match_obj, re);
        if (matched) {
            comma_delimited_value_string = match_obj[1];
            fixed = true;
        }   
        
        if (comma_delimited_value_string == "equal") {
            vector_of_values.resize(1);
            vector_of_values[0] = -1;
        }
        else {
            // Convert comma_delimited_value_string to vector_of_strings
            std::vector<std::string> vector_of_strings;
            boost::split(vector_of_strings, comma_delimited_value_string, boost::is_any_of(","));

            // Convert vector_of_strings to vector_of_values (single values in case of ratevar, ncateg, and pinvar)
            vector_of_values.resize(vector_of_strings.size());
            std::transform(
                vector_of_strings.begin(),      // start of source vector
                vector_of_strings.end(),        // end of source vector
                vector_of_values.begin(),       // start of destination vector
                [](const std::string & vstr) {  // anonymous function that translates
                    return std::stod(vstr);     // each string element to a double
                }
            );
            assert(vector_of_values.size() > 0);
        }
        
        // Split comma_delimited_subset_names_string into vector_of_subset_names
        boost::split(vector_of_subset_names, comma_delimited_subset_names_string, boost::is_any_of(","));
        
        // Sanity check: at least one subset name must be provided
        if (vector_of_subset_names.size() == 0) {
            NamuX("At least 1 subset must be provided in assignments (use \"default\" if not partitioning)");
        }
        
        // Sanity check: if no partition was defined, then values should be assigned to "default" subset
        // and if "default" is in the list of subset names, it should be the only thing in that list
        unsigned num_subsets_defined = this->_partition->get_num_subsets();
        std::vector<std::string>::iterator default_iter = std::find(vector_of_subset_names.begin(), vector_of_subset_names.end(), std::string("default"));
        bool default_found = (default_iter != vector_of_subset_names.end());
        if (default_found) {
            if (vector_of_subset_names.size() > 1)
                throw NamuX("The \"default\" specification cannot be mixed with other subset-specific parameter specifications");
        }
        else if (num_subsets_defined == 0) {
            throw NamuX("Must specify partition before assigning values to particular subsets (or assign to subset \"default\")");
        }
        return fixed;
    }

    inline void Namu::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;

        try {
            std::cout << "Using underflow scaling: " << (_use_underflow_scaling ? "yes" : "no") << std::endl;
            std::cout << "\n*** Reading and storing the data in the file " << this->_data_file_path << std::endl;
            this->_data = Data::SharedPtr(new Data());
            this->_data->set_partition(this->_partition);
            this->_data->get_data_from_file(this->_data_file_path);

            this->_model->set_subset_num_patterns(this->_data->calc_num_patterns_vect());
            this->_model->set_subset_sizes(this->_partition->calc_subset_sizes());
            this->_model->activate();

            // Report information about data partition subsets
            unsigned nsubsets = this->_data->get_num_subsets();
            std::cout << "\nNumber of taxa: " << this->_data->get_num_taxa() << std::endl;
            std::cout << "Number of partition subsets: " << nsubsets << std::endl;
            for (unsigned subset = 0; subset < nsubsets; subset++) {
                DataType dt = this->_partition->get_data_type_for_subset(subset);
                std::cout << "  Subset " << (subset+1) << " (" << _data->get_subset_name(subset) << ")" << std::endl;
                std::cout << "    data type: " << dt.get_data_type_as_string() << std::endl;
                std::cout << "    sites:     " << _data->calc_seq_len_in_subset(subset) << std::endl;
                std::cout << "    patterns:  " << _data->get_num_patterns_in_subset(subset) << std::endl;
                std::cout << "    ambiguity: " << (_ambig_missing || dt.is_codon() ? "treated as missing data (faster)" : "handled appropriately (slower)") << std::endl;
            }

            std::cout << "\n*** Resources available to BeagleLib " << this->_likelihood->beagle_lib_version() << ":\n";
            std::cout << this->_likelihood->available_resources() << std::endl;

            std::cout << "\n*** Creating the likelihood calculator" << std::endl;
            this->_likelihood = Likelihood::SharedPtr(new Likelihood());
            this->_likelihood->set_prefer_gpu(this->_use_gpu);
            this->_likelihood->set_ambiguity_equals_missing(this->_ambig_missing);
            this->_likelihood->set_data(this->_data);
            this->_likelihood->use_underflow_scaling(this->_use_underflow_scaling);

            std::cout << "\n*** Model description" << std::endl;
            std::cout << this->_model->describe_model() << std::endl;
            this->_likelihood->set_model(this->_model);

            this->_likelihood->init_beagle_lib();

            std::cout << "\nUsing underflow scaling: " << (this->_likelihood->using_underflow_scaling() ? "yes" : "no") << std::endl;

            std::cout << "\n*** Reading and storing the first tree in the file " << this->_tree_file_path << std::endl;
            this->_tree_summary = TreeSummary::SharedPtr(new TreeSummary());
            this->_tree_summary->read_tree_file(this->_tree_file_path, 0);
            Tree::SharedPtr tree = this->_tree_summary->get_tree(0);

            if (tree->num_leaves() != this->_data->get_num_taxa())
                throw NamuX(boost::format("Number of taxa in tree (%d) does not equal the number of taxa in the data matrix (%d)") % tree->num_leaves() % this->_data->get_num_taxa());

            std::cout << "\n*** Calculating the likelihood of the tree" << std::endl;
            double lnL = this->_likelihood->calc_log_likelihood(tree);
            std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;

            if (this->_expected_log_likelihood != 0.0)
                std::cout << boost::str(boost::format("      (expecting %.3f)") % this->_expected_log_likelihood) << std::endl;

            std::cout << "\nUsing underflow scaling: " << (this->_likelihood->using_underflow_scaling() ? "yes" : "no") << std::endl;
        }
        catch (NamuX & x) {
            std::cerr << "Namu encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

}
