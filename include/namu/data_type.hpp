#pragma once

#include <boost/format.hpp>
#include "namu/genetic_code.hpp"

namespace namu {

    class DataType {
        public:
                                            DataType();
                                            ~DataType();

            void                            set_nucleotide();
            void                            set_codon();
            void                            set_protein();
            void                            set_standard();
        
            bool                            is_nucleotide() const;
            bool                            is_codon() const;
            bool                            is_protein() const;
            bool                            is_standard() const;

            void                            set_standard_num_states(unsigned nstates);
            void                            set_genetic_code_from_name(std::string genetic_code_name);
            void                            set_genetic_code(GeneticCode::SharedPtr gcode);

            unsigned                        get_data_type() const;
            unsigned                        get_num_states() const;
            std::string                     get_data_type_as_string() const;
            const GeneticCode::SharedPtr    get_genetic_code() const;
            
            static std::string              translate_data_type_to_string(unsigned data_type);

        private:

            unsigned                        _data_type;
            unsigned                        _num_states;
            GeneticCode::SharedPtr          _genetic_code;
    };
    
    inline DataType::DataType() : _data_type(0), _num_states(0) {
        //std::cout << "Creating a data_type object" << std::endl;
        this->set_nucleotide();
    }
    
    inline DataType::~DataType() {
        //std::cout << "Destroying a data_type object" << std::endl;
    }
    
    inline void DataType::set_nucleotide() {
        this->_data_type = 1;
        this->_num_states = 4;
        this->_genetic_code = nullptr;
    }
    
    inline void DataType::set_codon() {
        this->_data_type = 2;
        this->_genetic_code = GeneticCode::SharedPtr(new GeneticCode("standard"));
        this->_num_states = this->_genetic_code->get_num_nonstop_codons();
    }
    
    inline void DataType::set_protein() {
        this->_data_type = 3;
        this->_num_states = 20;
        this->_genetic_code = nullptr;
    }
    
    inline void DataType::set_standard() {
        this->_data_type = 4;
        this->_num_states = 2;
        this->_genetic_code = nullptr;
    }

    inline bool DataType::is_nucleotide() const {
        return (this->_data_type == 1);
    }

    inline bool DataType::is_codon() const {
        return (this->_data_type == 2);
    }

    inline bool DataType::is_protein() const {
        return (this->_data_type == 3);
    }

    inline bool DataType::is_standard() const {
        return (this->_data_type == 4);
    }

    inline void DataType::set_genetic_code_from_name(std::string genetic_code_name) {
        assert(this->is_codon());
        this->_genetic_code = GeneticCode::SharedPtr(new GeneticCode(genetic_code_name));
    }
    
    inline void DataType::set_genetic_code(GeneticCode::SharedPtr gcode) {
        assert(this->is_codon());
        assert(gcode);
        this->_genetic_code = gcode;
    }

    inline void DataType::set_standard_num_states(unsigned nstates) {
        this->_data_type = 4;
        this->_num_states = nstates;
        this->_genetic_code = nullptr;
    }

    inline unsigned DataType::get_data_type() const {
        return this->_data_type;
    }
    
    inline unsigned DataType::get_num_states() const {
        return this->_num_states;
    }
    
    inline const GeneticCode::SharedPtr DataType::get_genetic_code() const {
        assert(this->is_codon());
        return this->_genetic_code;
    }
    
    inline std::string DataType::get_data_type_as_string() const {
        std::string s = translate_data_type_to_string(this->_data_type);
        if (this->is_codon())
            s += boost::str(boost::format(",%s") % this->_genetic_code->get_genetic_code_name());
        return s;
    }
    
    inline std::string DataType::translate_data_type_to_string(unsigned data_type) {
        assert(data_type == 1 || data_type == 2 || data_type == 3 || data_type == 4);
        if (data_type == 1)
            return std::string("nucleotide");
        else if (data_type == 2)
            return std::string("codon");
        else if (data_type == 3)
            return std::string("protein");
        else
            return std::string("standard");
    }
}
