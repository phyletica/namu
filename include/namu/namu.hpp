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

            void                clear();
            void                process_command_line_options(int argc, const char * argv[]);
            void                run();

        private:

            std::string             _conf_file_path;
            std::string             _data_file_path;
            std::string             _tree_file_path;

            Partition::SharedPtr    _partition;
            Data::SharedPtr         _data;
            Likelihood::SharedPtr   _likelihood;

            TreeSummary::SharedPtr  _tree_summary;

            bool                    _use_gpu;
            bool                    _ambig_missing;

            static std::string      _program_name;
            static unsigned         _major_version;
            static unsigned         _minor_version;

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
        this->_data = nullptr; 
        this->_likelihood = nullptr;
    }

    inline void Namu::process_command_line_options(int argc, const char * argv[]) {
        std::vector<std::string> partition_subsets;
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
            ("gpu",           boost::program_options::value(&_use_gpu)->default_value(true), "Use GPU if available")
            ("ambigmissing",  boost::program_options::value(&_ambig_missing)->default_value(true), "Treat all ambiguities as missing data")
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
    }

    inline void Namu::run() {
        std::cout << "Starting..." << std::endl;
        std::cout << "Current working directory: " << boost::filesystem::current_path() << std::endl;

        try {
            std::cout << "\n*** Reading and storing the data in the file " << this->_data_file_path << std::endl;
            this->_data = Data::SharedPtr(new Data());
            this->_data->set_partition(this->_partition);
            this->_data->get_data_from_file(this->_data_file_path);

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
            this->_likelihood->set_data(_data);
            this->_likelihood->init_beagle_lib();

            std::cout << "\n*** Reading and storing the first tree in the file " << this->_tree_file_path << std::endl;
            this->_tree_summary = TreeSummary::SharedPtr(new TreeSummary());
            this->_tree_summary->read_tree_file(this->_tree_file_path, 0);
            Tree::SharedPtr tree = this->_tree_summary->get_tree(0);

            if (tree->num_leaves() != this->_data->get_num_taxa())
                throw NamuX(boost::format("Number of taxa in tree (%d) does not equal the number of taxa in the data matrix (%d)") % tree->num_leaves() % this->_data->get_num_taxa());

            std::cout << "\n*** Calculating the likelihood of the tree" << std::endl;
            double lnL = this->_likelihood->calc_log_likelihood(tree);
            std::cout << boost::str(boost::format("log likelihood = %.5f") % lnL) << std::endl;
            std::cout << "      (expecting -278.83767)" << std::endl;
        }
        catch (NamuX & x) {
            std::cerr << "Namu encountered a problem:\n  " << x.what() << std::endl;
        }

        std::cout << "\nFinished!" << std::endl;
    }

}
