#pragma once

#include <boost/algorithm/string.hpp>
#include "namu/error.hpp"

namespace namu {

    class Data;
    class Model;
    class QMatrix;

    class GeneticCode {

        friend class Data;
        friend class Model;
        friend class QMatrix;

        public:

            typedef std::map<int, int>                  genetic_code_map_t;
            typedef std::map<char, unsigned>            amino_acid_map_t;
            typedef std::vector<unsigned>               amino_acid_vect_t;
            typedef std::vector<std::string>            codon_vect_t;
            typedef std::vector<char>                   amino_acid_symbol_vect_t;
            typedef std::map<std::string, std::string>  genetic_code_definitions_t;
            typedef std::vector<std::string>            genetic_code_names_t;
            
        
                                                GeneticCode();
                                                GeneticCode(std::string name);
                                                ~GeneticCode();
        
            std::string                         get_genetic_code_name() const;
            void                                use_genetic_code(std::string name);
        
            unsigned                            get_num_nonstop_codons() const;
            int                                 get_state_code(int triplet_index) const;
            char                                get_amino_acid_abbrev(unsigned aa_index) const;

            void                                copy_codons(std::vector<std::string> & codon_vect) const;
            void                                copy_amino_acids(std::vector<unsigned> & aa_vect) const;

            static genetic_code_names_t         get_recognized_genetic_code_names();
            static bool                         is_recognized_genetic_code_name(const std::string & name);
            static void                         ensure_genetic_code_name_is_valid(const std::string & name);

        private:

            void                                build_genetic_code_translators();

            std::string                         _genetic_code_name;
        
            genetic_code_map_t                  _genetic_code_map;
            amino_acid_map_t                    _amino_acid_map;
        
            amino_acid_vect_t                   _amino_acids;
            codon_vect_t                        _codons;
        
            const amino_acid_symbol_vect_t      _all_amino_acids = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'};
            const std::vector<std::string>      _all_codons = {"AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"};

            static genetic_code_definitions_t   _definitions;

        public:

            typedef std::shared_ptr<GeneticCode> SharedPtr;
    };

    inline GeneticCode::GeneticCode() {
        //std::cout << "Constructing a standard GeneticCode" << std::endl;
        this->use_genetic_code("standard");
    }

    inline GeneticCode::GeneticCode(std::string name) {
        //std::cout << "Constructing a " << name << " GeneticCode" << std::endl;
        this->use_genetic_code(name);
    }

    inline GeneticCode::~GeneticCode() {
        //std::cout << "Destroying a GeneticCode" << std::endl;
    }

    inline std::string GeneticCode::get_genetic_code_name() const {
        return this->_genetic_code_name;
    }

    inline void GeneticCode::use_genetic_code(std::string name) {
        this->_genetic_code_name = name;
        this->build_genetic_code_translators();
    }

    inline unsigned GeneticCode::get_num_nonstop_codons() const {
        return (unsigned)this->_codons.size();
    }

    inline int GeneticCode::get_state_code(int triplet_index) const {
        return _genetic_code_map.at(triplet_index);
    }
    
    inline char GeneticCode::get_amino_acid_abbrev(unsigned aa_index) const {
        return this->_all_amino_acids[aa_index];
    }

    inline void GeneticCode::copy_codons(std::vector<std::string> & codon_vect) const {
        codon_vect.resize(this->_codons.size());
        std::copy(this->_codons.begin(), this->_codons.end(), codon_vect.begin());
    }
    
    inline void GeneticCode::copy_amino_acids(std::vector<unsigned> & aa_vect) const {
        aa_vect.resize(this->_amino_acids.size());
        std::copy(this->_amino_acids.begin(), this->_amino_acids.end(), aa_vect.begin());
    }

    inline void GeneticCode::build_genetic_code_translators() {
        this->_amino_acid_map.clear();
        for (unsigned i = 0; i < 20; ++i) {
            char aa = this->_all_amino_acids[i];
            this->_amino_acid_map[aa] = i;
        }

        this->ensure_genetic_code_name_is_valid(this->_genetic_code_name);
        std::string gcode_aa = this->_definitions[this->_genetic_code_name];  // e.g. "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
        
        int k = 0;
        int state_code = 0;
        this->_codons.clear();
        this->_amino_acids.clear();
        this->_genetic_code_map.clear();
        for (char ch : gcode_aa) {
            if (ch != '*') {
                this->_genetic_code_map[k] = state_code++;
                this->_codons.push_back(this->_all_codons[k]);
                this->_amino_acids.push_back(this->_amino_acid_map[ch]);
            }
            ++k;
        }
    }

    inline GeneticCode::genetic_code_names_t GeneticCode::get_recognized_genetic_code_names() {
        genetic_code_names_t names;
        for (auto it = _definitions.begin(); it != _definitions.end(); ++it)
            names.push_back(it->first);
        return names;
    }
    
    inline bool GeneticCode::is_recognized_genetic_code_name(const std::string & name) {
        std::string lcname = name;
        boost::to_lower(lcname);
        return (_definitions.find(lcname) != _definitions.end());
    }
   
    inline void GeneticCode::ensure_genetic_code_name_is_valid(const std::string & name) {
        if (! is_recognized_genetic_code_name(name)) {
            auto valid_genetic_code_names = get_recognized_genetic_code_names();
            std::cout << "Recognized genetic codes:\n";
            for (std::string name : valid_genetic_code_names) {
                std::cout << "  " << name << "\n";
            }
            std::cout << std::endl;
            throw NamuX(boost::format("%s is not a recognized genetic code") % name);
        }
    }
}
