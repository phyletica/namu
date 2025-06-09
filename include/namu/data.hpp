#pragma once

#include <fstream>
#include <regex>
#include <string>
#include <vector>
#include <numeric>
#include <limits>
#include <map>
#include <boost/format.hpp>
#include <boost/algorithm/string/join.hpp>
#include <ncl/nxsmultiformat.h>
#include "namu/genetic_code.hpp"
#include "namu/data_type.hpp"
#include "namu/partition.hpp"
#include "namu/error.hpp"

namespace namu {

    class Data {
        public:
            typedef std::vector<std::string>            taxon_names_t;
            typedef unsigned long long                  state_t;
            typedef std::vector<state_t>                pattern_vect_t;
            typedef std::vector<state_t>                monomorphic_vect_t;
            typedef std::vector<int>                    partition_key_t;
            typedef std::map<pattern_vect_t,unsigned>   pattern_map_t;
            typedef std::vector<pattern_vect_t>         data_matrix_t;
            typedef std::vector<pattern_map_t>          pattern_map_vect_t;
            typedef std::vector<double>                 pattern_counts_t;
            typedef std::vector<unsigned>               subset_end_t;
            typedef std::vector<unsigned>               npatterns_vect_t;
            typedef std::pair<unsigned, unsigned>       begin_end_pair_t;
            typedef std::shared_ptr<Data>               SharedPtr;

                                                        Data();
                                                        ~Data();
        
            Partition::SharedPtr                        get_partition();
            void                                        set_partition(Partition::SharedPtr partition);

            void                                        get_data_from_file(const std::string file_name);

            unsigned                                    get_num_subsets() const;
            std::string                                 get_subset_name(unsigned subset) const;

            unsigned                                    get_num_taxa() const;
            const taxon_names_t &                       get_taxon_names() const;

            unsigned                                    get_num_patterns() const;
            npatterns_vect_t                            calc_num_patterns_vect() const;
            unsigned                                    get_num_patterns_in_subset(unsigned subset) const;
            unsigned                                    get_num_states_for_subset(unsigned subset) const;
            unsigned                                    calc_seq_len() const;
            unsigned                                    calc_seq_len_in_subset(unsigned subset) const;
            const data_matrix_t &                       get_data_matrix() const;
            begin_end_pair_t                            get_subset_begin_end(unsigned subset) const;
            const pattern_counts_t &                    get_pattern_counts() const;
            const monomorphic_vect_t &                  get_monomorphic() const;
            const partition_key_t &                     get_partition_key() const;


            void                                        clear();

        private:

            unsigned                                    store_taxon_names(NxsTaxaBlock * taxaBlock, unsigned taxa_block_index);
            unsigned                                    store_data(unsigned ntax, unsigned nchar, NxsCharactersBlock * char_block, NxsCharactersBlock::DataTypesEnum data_type);
            unsigned                                    build_subset_specific_maps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets);
            void                                        update_pattern_map(Data::pattern_vect_t & pattern, unsigned subset);
            void                                        compress_patterns();

            Partition::SharedPtr                        _partition;
            pattern_counts_t                            _pattern_counts;
            monomorphic_vect_t                          _monomorphic;
            partition_key_t                             _partition_key;
            pattern_map_vect_t                          _pattern_map_vect;
            taxon_names_t                               _taxon_names;
            data_matrix_t                               _data_matrix;
            subset_end_t                                _subset_end;
    };

    inline Data::Data() {
        //std::cout << "Creating a Data object" << std::endl;
        this->clear();
    }

    inline Data::~Data() {
        //std::cout << "Destroying a Data object" << std::endl;
    }
    
    inline void Data::set_partition(Partition::SharedPtr partition) {
        this->_partition = partition;
    }

    inline Partition::SharedPtr Data::get_partition() {
        return this->_partition;
    }

    inline unsigned Data::get_num_subsets() const {
        return (this->_partition ? this->_partition->get_num_subsets() : 1);
    }
    
    inline std::string Data::get_subset_name(unsigned subset) const {
        return this->_partition ? this->_partition->get_subset_name(subset) : std::string("default");
    }

    inline const Data::partition_key_t & Data::get_partition_key() const {
        return this->_partition_key;
    }
    
    inline const Data::pattern_counts_t & Data::get_pattern_counts() const {
        return this->_pattern_counts;
    }
    
    inline const Data::monomorphic_vect_t & Data::get_monomorphic() const {
        return this->_monomorphic;
    }

    inline const Data::taxon_names_t & Data::get_taxon_names() const {
        return this->_taxon_names;
    }

    inline const Data::data_matrix_t & Data::get_data_matrix() const {
        return this->_data_matrix;
    }

    inline Data::begin_end_pair_t Data::get_subset_begin_end(unsigned subset) const {
        assert(this->_subset_end.size() > subset);
        if (subset == 0)
            return std::make_pair(0, this->_subset_end[0]);
        else
            return std::make_pair(this->_subset_end[subset-1], this->_subset_end[subset]);
    }

    inline void Data::clear() {
        this->_partition_key.clear();
        this->_pattern_counts.clear();
        this->_monomorphic.clear();
        this->_pattern_map_vect.clear();
        this->_taxon_names.clear();
        this->_data_matrix.clear();
        this->_subset_end.clear();
    }

    inline unsigned Data::get_num_patterns() const {
        if (this->_data_matrix.size() > 0)
            return (unsigned)this->_data_matrix[0].size();
        else
            return 0;
    }

    inline Data::npatterns_vect_t Data::calc_num_patterns_vect() const {
        unsigned nsubsets = (unsigned)this->_subset_end.size();
        std::vector<unsigned> num_patterns_vect(nsubsets, 0);
        for (unsigned s = 0; s < nsubsets; s++)
            num_patterns_vect[s] = this->get_num_patterns_in_subset(s);
        return num_patterns_vect;
    }
    
    inline unsigned Data::get_num_states_for_subset(unsigned subset) const {
        DataType data_type = this->_partition->get_data_type_for_subset(subset);
        return data_type.get_num_states();
    }

    inline unsigned Data::get_num_patterns_in_subset(unsigned subset) const {
        assert(this->_subset_end.size() > subset);
        return (unsigned)this->_subset_end[subset] - (subset == 0 ? 0 : this->_subset_end[subset-1]);
    }
    
    inline unsigned Data::get_num_taxa() const {
        return (unsigned)this->_taxon_names.size();
    }

    inline unsigned Data::calc_seq_len() const {
        return std::accumulate(this->_pattern_counts.begin(), this->_pattern_counts.end(), 0);
    }

    inline unsigned Data::calc_seq_len_in_subset(unsigned subset) const {
        begin_end_pair_t s = this->get_subset_begin_end(subset);
        return std::accumulate(this->_pattern_counts.begin() + s.first, this->_pattern_counts.begin() + s.second, 0);
    }
    
    inline unsigned Data::build_subset_specific_maps(unsigned ntaxa, unsigned seqlen, unsigned nsubsets) {
        pattern_vect_t pattern(ntaxa);

        this->_pattern_map_vect.clear();
        this->_pattern_map_vect.resize(nsubsets);
        
        const Partition::partition_t & tuples = this->_partition->get_subset_range_vect();
        for (auto & t : tuples) {
            unsigned site_begin  = std::get<0>(t);
            unsigned site_end    = std::get<1>(t);
            unsigned site_skip   = std::get<2>(t);
            unsigned site_subset = std::get<3>(t);
            for (unsigned site = site_begin; site <= site_end; site += site_skip) {
                // Copy site into pattern
                for (unsigned taxon = 0; taxon < ntaxa; ++taxon) {
                    pattern[taxon] = this->_data_matrix[taxon][site-1];
                }
                
                // Add this pattern to _pattern_map_vect element corresponding to subset site_subset
                this->update_pattern_map(pattern, site_subset);
            }
        }
        
        // Tally total number of patterns across all subsets
        unsigned npatterns = 0;
        for (auto & map : this->_pattern_map_vect) {
            npatterns += (unsigned)map.size();
        }
        
        return npatterns;
    }

    inline void Data::update_pattern_map(Data::pattern_vect_t & pattern, unsigned subset) {
        // If pattern is not already in pattern_map, insert it and set value to 1.
        // If it does exist, increment its current value.
        // (see item 24, p. 110, in Meyers' Efficient STL for more info on the technique used here)
        pattern_map_t::iterator lowb = this->_pattern_map_vect[subset].lower_bound(pattern);
        if (lowb != this->_pattern_map_vect[subset].end() && !(this->_pattern_map_vect[subset].key_comp()(pattern, lowb->first))) {
            // this pattern has already been seen
            lowb->second += 1;
        }
        else {
            // this pattern has not yet been seen
            this->_pattern_map_vect[subset].insert(lowb, pattern_map_t::value_type(pattern, 1));
        }
    }

    inline void Data::compress_patterns() {
        // Perform sanity checks
        if (this->_data_matrix.empty())
            throw NamuX("Attempted to compress an empty data matrix");
        
        unsigned ntaxa = (unsigned)this->_data_matrix.size();
        unsigned seqlen = (unsigned)this->_data_matrix[0].size();
        
        // Finalize partition
        unsigned nsubsets = this->get_num_subsets();
        this->_subset_end.resize(nsubsets);
        this->_partition->finalize(seqlen);

        // Compact the data, storing it in _pattern_map_vect
        unsigned npatterns = this->build_subset_specific_maps(ntaxa, seqlen, nsubsets);
        this->_pattern_counts.assign(npatterns, 0);
        this->_monomorphic.assign(npatterns, 0);
        this->_partition_key.assign(npatterns, -1);

        // Rebuild _data_matrix to hold compact data, storing counts in _pattern_counts
        this->_data_matrix.resize(ntaxa);
        for (auto & row : this->_data_matrix) {
            row.resize(npatterns);
        }

        unsigned p = 0;
        for (unsigned subset = 0; subset < nsubsets; subset++) {
            for (auto & pc : this->_pattern_map_vect[subset]) {
                this->_pattern_counts[p] = pc.second; // record how many sites have pattern p
                this->_partition_key[p] = subset;     // record the subset to which pattern p belongs
                
                state_t constant_state = pc.first[0];
                unsigned t = 0;
                for (auto sc : pc.first) {
                    assert(sc > 0);
                    constant_state &= sc;
                    this->_data_matrix[t][p] = sc;
                    ++t;
                }
                // constant_state equals 0 if polymorphic or state code of state present if monomorphic
                this->_monomorphic[p] = constant_state;
                ++p;
            }
            
            this->_subset_end[subset] = p;

            // Everything for this subset has been transferred to _data_matrix and _pattern_counts,
            // so we can now free this memory
            this->_pattern_map_vect[subset].clear();
        }
    }

    inline unsigned Data::store_taxon_names(NxsTaxaBlock * taxaBlock, unsigned taxa_block_index) {
        unsigned ntax = 0;
        if (taxa_block_index == 0) {
            // First taxa block encountered in the file
            _taxon_names.clear();
            for (auto s : taxaBlock->GetAllLabels())
                this->_taxon_names.push_back(s);
            ntax = (unsigned)this->_taxon_names.size();
            this->_data_matrix.resize(ntax);
        }
        else {
            // Second (or later) taxa block encountered in the file
            // Check to ensure taxa block is identical to the first one
            for (auto s : taxaBlock->GetAllLabels()) {
                if (this->_taxon_names[ntax++] != s)
                    throw NamuX(boost::format("Taxa block %d in data file is not identical to first taxa block read") % (taxa_block_index+1));
            }
        }
        
        return ntax;
    }
            
    inline unsigned Data::store_data(unsigned ntax, unsigned nchar_before, NxsCharactersBlock * char_block, NxsCharactersBlock::DataTypesEnum data_type) {
        unsigned seqlen = 0;
        
        // Find the data type for the partition subset containing the first site in this NxsCharactersBlock
        // Assumes that all sites in any given NxsCharactersBlock have the same type (i.e. mixed not allowed)
        assert(this->_partition);
        unsigned subset_index = this->_partition->find_subset_for_site(nchar_before + 1); // remember that sites begin at 1, not 0, in partition definitions
        DataType dt = this->_partition->get_data_type_for_subset(subset_index);

        // Determine number of states and bail out if data type not handled
        // 1 = standard, 2 = dna, 3 = rna, 4 = nucleotide, 5 = protein, 6 = continuous, 7 = codon, 8 = mixed
        NxsCharactersBlock * block = char_block;
        if (data_type == NxsCharactersBlock::dna || data_type == NxsCharactersBlock::rna || data_type == NxsCharactersBlock::nucleotide) {
            if (dt.is_codon()) {
                // Create a NxsCharactersBlock containing codons rather than nucleotides
                block = NxsCharactersBlock::NewCodonsCharactersBlock(
                    char_block,
                    true,   // map partial ambiguities to completely missing (note: false is not yet implemented in NCL)
                    true,   // gaps to missing
                    true,   // inactive characters treated as missing
                    NULL,   // if non-NULL, specifies the indices of the positions in the gene
                    NULL);  // if non-NULL, specifies a pointer to a NxsCharactersBlock that contains all non-coding positions in gene
            }
            else {
                if (!dt.is_nucleotide())
                    throw NamuX(boost::format("Partition subset has data type \"%s\" but data read from file has data type \"nucleotide\"") % dt.get_data_type_as_string());
            }
        }
        else if (data_type == NxsCharactersBlock::protein) {
            if (!dt.is_protein())
                throw NamuX(boost::format("Partition subset has data type \"%s\" but data read from file has data type \"protein\"") % dt.get_data_type_as_string());
        }
        else if (data_type == NxsCharactersBlock::standard) {
            if (!dt.is_standard())
                throw NamuX(boost::format("Partition subset has data type \"%s\" but data read from file has data type \"standard\"") % dt.get_data_type_as_string());
            assert(char_block->GetSymbols());
            std::string symbols = std::string(char_block->GetSymbols());
            dt.set_standard_num_states((unsigned)symbols.size());
        }
        else {
            // ignore block because data type is not one that is supported
            return nchar_before;
        }
        
        unsigned num_states = dt.get_num_states();
        
        // Make sure all states can be accommodated in a variable of type state_t
        unsigned bits_in_state_t = 8*sizeof(state_t);
        if (num_states > bits_in_state_t)
            throw NamuX(boost::format("This program can only process data types with fewer than %d states") % bits_in_state_t);
        
        // Copy data matrix from NxsCharactersBlock object to _data_matrix
        // Loop through all taxa, processing one row from block for each taxon
        for (unsigned t = 0; t < ntax; ++t) {

            const NxsDiscreteStateRow & row = block->GetDiscreteMatrixRow(t);
            if (seqlen == 0)
                seqlen = (unsigned)row.size();
            this->_data_matrix[t].resize(nchar_before + seqlen);
            
            // Loop through all sites/characters in row corresponding to taxon t
            unsigned k = nchar_before;
            for (int raw_state_code : row) {
                // For codon model, raw_state_code ranges from 0-63, but deletion of stop codons means fewer state codes
                state_t state = std::numeric_limits<state_t>::max(); // complete ambiguity, all bits set
                bool complete_ambiguity = (!dt.is_codon() && raw_state_code == (int)num_states);
                bool all_missing_or_gaps = (raw_state_code < 0);
                if ((!complete_ambiguity) && (!all_missing_or_gaps)) {
                    int state_code = raw_state_code;
                    if (dt.is_codon())
                        state_code = dt.get_genetic_code()->get_state_code(raw_state_code);

                    if (state_code < (int)num_states) {
                        state = (state_t)1 << state_code;
                    }
                    else {
                        // incomplete ambiguity (NCL state code > num_states)
                        const NxsDiscreteDatatypeMapper      * mapper = block->GetDatatypeMapperForChar(k - nchar_before);
                        const std::set<NxsDiscreteStateCell> & state_set = mapper->GetStateSetForCode(raw_state_code);
                        state = 0;
                        for (auto s : state_set) {
                             state |= (state_t)1 << s;
                        }
                    }
                }
                this->_data_matrix[t][k++] = state;
            }
        }
        
        return seqlen;
    }

    inline void Data::get_data_from_file(const std::string file_name) {
        // See http://phylo.bio.ku.edu/ncldocs/v2.1/funcdocs/index.html for documentation
        //
        // -1 means "process all blocks found" (this is a bit field and -1 fills the bit field with 1s)
        // Here are the bits (and nexus blocks) that are defined:
        //     enum NexusBlocksToRead
        //     {
        //         NEXUS_TAXA_BLOCK_BIT = 0x01,
        //         NEXUS_TREES_BLOCK_BIT = 0x02,
        //         NEXUS_CHARACTERS_BLOCK_BIT = 0x04,
        //         NEXUS_ASSUMPTIONS_BLOCK_BIT = 0x08,
        //         NEXUS_SETS_BLOCK_BIT = 0x10,
        //         NEXUS_UNALIGNED_BLOCK_BIT = 0x20,
        //         NEXUS_DISTANCES_BLOCK_BIT = 0x40,
        //         NEXUS_UNKNOWN_BLOCK_BIT = 0x80
        //     };
        MultiFormatReader nexus_reader(-1, NxsReader::WARNINGS_TO_STDERR);
        try {
            nexus_reader.ReadFilepath(file_name.c_str(), MultiFormatReader::NEXUS_FORMAT);
        }
        catch(...) {
            nexus_reader.DeleteBlocksFromFactories();
            throw;
        }

        // Commit to storing new data
        this->clear();

        // Ensure that Data::set_partition was called before reading data
        assert(this->_partition);

        int num_taxa_blocks = nexus_reader.GetNumTaxaBlocks();
        if (num_taxa_blocks == 0)
            throw NamuX("No taxa blocks were found in the data file");
            
        unsigned cum_nchar = 0;
        for (int i = 0; i < num_taxa_blocks; ++i) {
            NxsTaxaBlock * taxa_block = nexus_reader.GetTaxaBlock(i);
            unsigned ntax = this->store_taxon_names(taxa_block, i);
            const unsigned num_char_blocks = nexus_reader.GetNumCharactersBlocks(taxa_block);
            for (unsigned j = 0; j < num_char_blocks; ++j) {
                NxsCharactersBlock * char_block = nexus_reader.GetCharactersBlock(taxa_block, j);
                NxsCharactersBlock::DataTypesEnum data_type = char_block->GetOriginalDataType();
                cum_nchar += this->store_data(ntax, cum_nchar, char_block, data_type);
            }
        }

        // No longer any need to store raw data from nexus file
        nexus_reader.DeleteBlocksFromFactories();

        // Compress _data_matrix so that it holds only unique patterns (counts stored in _pattern_counts)
        if (this->_data_matrix.empty()) {
            std::cout << "No data were stored from the file \"" << file_name << "\"" << std::endl;
            clear();
        }
        else {
            this->compress_patterns();
        }
    }
}
