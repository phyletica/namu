#pragma once

#include <algorithm>
#include <vector>
#include <Eigen/Dense>
#include <boost/math/distributions/gamma.hpp>
#include "libhmsbeagle/beagle.h"
#include "namu/data_type.hpp"
#include "namu/qmatrix.hpp"
#include "namu/asrv.hpp"

namespace namu {

    class Likelihood;

    class Model {

        friend class Likelihood;

        public:
            typedef std::vector<ASRV::SharedPtr>      asrv_vect_t;
            typedef std::vector<QMatrix::SharedPtr>   qmat_vect_t;
            typedef std::vector<unsigned>             subset_sizes_t;
            typedef std::vector<DataType>             subset_datatype_t;
            typedef std::vector<double>               subset_relrate_vect_t;
            typedef boost::shared_ptr<Model>          SharedPtr;

                                        Model();
                                        ~Model();

            void                        activate();
            void                        inactivate();

            std::string                 describe_model();

            void                        set_subset_data_types(const subset_datatype_t & datatype_vect);

            void                        set_subset_rate_var(ASRV::ratevar_ptr_t ratevar, unsigned subset, bool fixed);
            void                        set_subset_pinvar(ASRV::pinvar_ptr_t pinvar, unsigned subset, bool fixed);
            void                        set_subset_exchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed);
            void                        set_subset_state_freqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed);
            void                        set_subset_omega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed);

            void                        set_subset_rel_rates(subset_relrate_vect_t & relrates, bool fixed);
            subset_relrate_vect_t &     get_subset_rel_rates();
            bool                        is_fixed_subset_rel_rates() const;
            double                      calc_normalizing_constant_for_subset_rel_rates() const;

            void                        set_tree_index(unsigned i, bool fixed);
            unsigned                    get_tree_index() const;
            bool                        is_fixed_tree() const;

            unsigned                    get_num_subsets() const;
            unsigned                    get_num_sites() const;
            unsigned                    get_subset_num_sites(unsigned subset) const;
            const QMatrix &             get_q_matrix(unsigned subset) const;
            const ASRV &                get_asrv(unsigned subset) const;

            void                        set_subset_is_invar_model(bool is_invar, unsigned subset);
            bool                        get_subset_is_invar_model(unsigned subset) const;

            void                        set_subset_sizes(const subset_sizes_t nsites_vect);
            subset_sizes_t &            get_subset_sizes();

            void                        set_subset_num_patterns(const subset_sizes_t npatterns_vect);
            unsigned                    get_subset_num_patterns(unsigned subset) const;

            void                        set_subset_num_categ(unsigned ncateg, unsigned subset);
            unsigned                    get_subset_num_categ(unsigned subset) const;

            int                         set_beagle_eigen_decomposition(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         set_beagle_state_frequencies(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         set_beagle_among_site_rate_variation_rates(int beagle_instance, unsigned subset, unsigned instance_subset);
            int                         set_beagle_among_site_rate_variation_probs(int beagle_instance, unsigned subset, unsigned instance_subset);

        private:

            void                        clear();

            unsigned                    _num_subsets;
            unsigned                    _num_sites;
            subset_sizes_t              _subset_sizes;
            subset_sizes_t              _subset_npatterns;
            subset_datatype_t           _subset_datatypes;
            qmat_vect_t                 _qmatrix;
            asrv_vect_t                 _asrv;

            unsigned                    _tree_index;
            bool                        _tree_fixed;

            bool                        _subset_relrates_fixed;
            subset_relrate_vect_t       _subset_relrates;
        };

    inline Model::Model() {
        //std::cout << "Constructing a Model" << std::endl;
        clear();
    }

    inline Model::~Model() {
        //std::cout << "Destroying a Model" << std::endl;
    }

    inline void Model::clear() {
        _num_subsets = 0;
        _num_sites = 0;
        _tree_index = 0;
        _tree_fixed = false;
        _subset_relrates_fixed = false;
        _subset_relrates.clear();
        _subset_sizes.clear();
        _subset_npatterns.clear();
        _subset_datatypes.clear();
        _qmatrix.clear();
        _asrv.clear();
    }

    inline std::string Model::describe_model() {
        // Creates summary such as following and returns as a string:
        //
        // Partition information:
        //
        //          data subset           1           2           3
        //    -----------------------------------------------------
        //           num. sites          20          20          20
        //        num. patterns           7           5          17
        //          num. states           4           4           4
        //      rate categories           4           1           4
        //
        // Parameter linkage:
        //
        //          data subset           1           2           3
        //    -----------------------------------------------------
        //          state freqs           1           1           1
        //    exchangeabilities           1           1           2
        //        rate variance           1           2           3
        //               pinvar           1           2           -

        // Sets used to determine which parameters are linked across subsets
        std::set<double *> freqset;
        std::set<double *> xchgset;
        std::set<double *> omegaset;
        std::set<double *> ratevarset;
        std::set<double *> pinvarset;
        std::set<double *> relrateset;

        // Vectors of pointers to distinct parameters
        std::vector<double *> unique_freq;
        std::vector<double *> unique_xchg;
        std::vector<double *> unique_omega;
        std::vector<double *> unique_ratevar;
        std::vector<double *> unique_pinvar;
        std::vector<double *> unique_relrate;

        // Map for storing strings that will contain the information for each row
        std::map<std::string, std::string> ss = {
            {"subset",    ""},
            {"dashes",    ""},
            {"freqs",     ""},
            {"xchg",      ""},
            {"omega",     ""},
            {"ratevar",   ""},
            {"pinvar",    ""},
            {"ncateg",    ""},
            {"nsites",    ""},
            {"npatterns", ""},
            {"nstates",   ""}
        };

        // Ensure that the subset relative rates are fixed if there is only one
        // subset; otherwise the subset relative rates will be added to the list
        // of free parameters that are updated, which makes no sense in this case
        if (_num_subsets == 1)
            _subset_relrates_fixed = true;

        // Loop through subsets, building up rows as we go
        for (unsigned i = 0; i < _num_subsets; i++) {
            // Ensure that for subsets in which the number of rate categories is 1 that
            // the gamma rate variance is fixed; otherwise the gamma rate variance will
            // be added to the list of free parameters that are updated, which makes
            // no sense in this case
            if (_asrv[i]->get_num_categ() == 1) {
                _asrv[i]->fix_rate_var(true);
            }

            unsigned index;
            ss["subset"] += boost::str(boost::format("%12d") % (i+1));
            ss["dashes"] += "------------";

            // Determine whether state freqs are unique for this subset
            QMatrix::freq_xchg_ptr_t pfreq = _qmatrix[i]->get_state_freqs_shared_ptr();
            QMatrix::freq_xchg_t & freq = *pfreq;
            double * freq_addr = &freq[0];
            auto f = freqset.insert(freq_addr);
            if (f.second) {
                unique_freq.push_back(freq_addr);
                index = (unsigned)unique_freq.size();
            }
            else {
                auto iter = std::find(unique_freq.begin(), unique_freq.end(), freq_addr);
                index = (unsigned)std::distance(unique_freq.begin(), iter) + 1;
            }
            ss["freqs"] += boost::str(boost::format("%12d") % index);

            // Determine whether exchangeabilities are unique for this subset
            if (_subset_datatypes[i].is_nucleotide()) {
                QMatrix::freq_xchg_ptr_t pxchg = _qmatrix[i]->get_exchangeabilities_shared_ptr();
                QMatrix::freq_xchg_t & xchg = *pxchg;
                double * xchg_addr = &xchg[0];
                auto x = xchgset.insert(xchg_addr);
                if (x.second) {
                    unique_xchg.push_back(xchg_addr);
                    index = (unsigned)unique_xchg.size();
                }
                else {
                    auto iter = std::find(unique_xchg.begin(), unique_xchg.end(), xchg_addr);
                    index = (unsigned)std::distance(unique_xchg.begin(), iter) + 1;
                }
                ss["xchg"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["xchg"] += boost::str(boost::format("%12s") % "-");
            }

            // Determine whether omega is unique for this subset
            if (_subset_datatypes[i].is_codon()) {
                QMatrix::omega_ptr_t pomega = _qmatrix[i]->get_omega_shared_ptr();
                QMatrix::omega_t omegavalue = *pomega;
                double * omega_addr = &omegavalue;
                auto o = omegaset.insert(omega_addr);
                if (o.second) {
                    unique_omega.push_back(omega_addr);
                    index = (unsigned)unique_omega.size();
                }
                else {
                    auto iter = std::find(unique_omega.begin(), unique_omega.end(), omega_addr);
                    index = (unsigned)std::distance(unique_omega.begin(), iter) + 1;
                }
                ss["omega"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["omega"] += boost::str(boost::format("%12s") % "-");
            }

            // Determine whether rate variance is unique for this subset
            ASRV::ratevar_ptr_t pratevar = _asrv[i]->get_rate_var_shared_ptr();
            double & ratevar = *pratevar;
            double * ratevar_addr = &ratevar;
            auto r = ratevarset.insert(ratevar_addr);
            if (r.second) {
                unique_ratevar.push_back(ratevar_addr);
                index = (unsigned)unique_ratevar.size();
            }
            else {
                auto iter = std::find(unique_ratevar.begin(), unique_ratevar.end(), ratevar_addr);
                index = (unsigned)std::distance(unique_ratevar.begin(), iter) + 1;
            }
            ss["ratevar"] += boost::str(boost::format("%12d") % index);

            // Determine whether pinvar is unique for this subset
            if (_asrv[i]->get_is_invar_model()) {
                ASRV::pinvar_ptr_t ppinvar = _asrv[i]->get_pinvar_shared_ptr();
                double & pinvar = *ppinvar;
                double * pinvar_addr = &pinvar;
                auto r = pinvarset.insert(pinvar_addr);
                if (r.second) {
                    unique_pinvar.push_back(pinvar_addr);
                    index = (unsigned)unique_pinvar.size();
                }
                else {
                    auto iter = std::find(unique_pinvar.begin(), unique_pinvar.end(), pinvar_addr);
                    index = (unsigned)std::distance(unique_pinvar.begin(), iter) + 1;
                }
                ss["pinvar"] += boost::str(boost::format("%12d") % index);
            }
            else {
                ss["pinvar"] += boost::str(boost::format("%12s") % "-");
            }

            // Save number of rate categories for this subset
            ss["ncateg"] += boost::str(boost::format("%12d") % _asrv[i]->get_num_categ());

            // Save number of sites for this subset
            ss["nsites"] += boost::str(boost::format("%12d") % _subset_sizes[i]);

            // Save number of patterns for this subset
            ss["npatterns"] += boost::str(boost::format("%12d") % _subset_npatterns[i]);

            // Save number of states for this subset
            if (_subset_datatypes.size() == _num_subsets)
                ss["nstates"] += boost::str(boost::format("%12d") % _subset_datatypes[i].get_num_states());
            else
                ss["nstates"] += boost::str(boost::format("%12s") % "?");

        }
        std::string s = "Partition information:\n\n";

        s += boost::str(boost::format("%20s%s\n") % "data subset" % ss["subset"]);
        s += boost::str(boost::format("%20s%s\n") % "-----------------" % ss["dashes"]);
        s += boost::str(boost::format("%20s%s\n") % "num. sites" % ss["nsites"]);
        s += boost::str(boost::format("%20s%s\n") % "num. patterns" % ss["npatterns"]);
        s += boost::str(boost::format("%20s%s\n") % "num. states" % ss["nstates"]);
        s += boost::str(boost::format("%20s%s\n") % "rate categories" % ss["ncateg"]);

        s += "\nParameter linkage:\n\n";

        s += boost::str(boost::format("%20s%s\n") % "data subset" % ss["subset"]);
        s += boost::str(boost::format("%20s%s\n") % "-----------------" % ss["dashes"]);
        s += boost::str(boost::format("%20s%s\n") % "state freqs" % ss["freqs"]);
        s += boost::str(boost::format("%20s%s\n") % "exchangeabilities" % ss["xchg"]);
        s += boost::str(boost::format("%20s%s\n") % "omega" % ss["omega"]);
        s += boost::str(boost::format("%20s%s\n") % "rate variance" % ss["ratevar"]);
        s += boost::str(boost::format("%20s%s\n") % "pinvar" % ss["pinvar"]);

        s += "\nParameter values for each subset:\n";

        s += "\n  relative rate:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            s += boost::str(boost::format("  %12d: %g\n") % (i+1) % _subset_relrates[i]);
        }

        s += "\n  state freqs:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            QMatrix::freq_xchg_t & freqs = *(_qmatrix[i]->get_state_freqs_shared_ptr());
            std::vector<std::string> freqs_as_strings(freqs.size());
            std::transform(freqs.begin(), freqs.end(), freqs_as_strings.begin(), [](double freq) {return boost::str(boost::format("%g") % freq);});
            std::string tmp = boost::algorithm::join(freqs_as_strings, ",");
            s += boost::str(boost::format("  %12d: (%s)\n") % (i+1) % tmp);
        }

        s += "\n  exchangeabilities:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_subset_datatypes[i].is_nucleotide()) {
                QMatrix::freq_xchg_t & xchg = *(_qmatrix[i]->get_exchangeabilities_shared_ptr());
                std::vector<std::string> xchg_as_strings(xchg.size());
                std::transform(xchg.begin(), xchg.end(), xchg_as_strings.begin(), [](double x) {return boost::str(boost::format("%g") % x);});
                std::string tmp = boost::algorithm::join(xchg_as_strings, ",");
                s += boost::str(boost::format("  %12d: (%s)\n") % (i+1) % tmp);
            }
            else {
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
            }
        }

        s += "\n  omega:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_subset_datatypes[i].is_codon()) {
                double omega = *(_qmatrix[i]->get_omega_shared_ptr());
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % omega);
            }
            else {
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
            }
        }

        s += "\n  rate variance:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            if (_asrv[i]->get_num_categ() > 1) {
                double ratevar = *(_asrv[i]->get_rate_var_shared_ptr());
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % ratevar);
            }
            else
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
        }

        s += "\n  pinvar:\n";
        for (unsigned i = 0; i < _num_subsets; i++) {
            double pinvar = *(_asrv[i]->get_pinvar_shared_ptr());
            bool is_invar_model = _asrv[i]->get_is_invar_model();
            if (is_invar_model)
                s += boost::str(boost::format("  %12d: %g\n") % (i+1) % pinvar);
            else
                s += boost::str(boost::format("  %12d: -\n") % (i+1));
        }

        return s;
    }

    inline unsigned Model::get_subset_num_patterns(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_npatterns[subset];
    }

    inline unsigned Model::get_subset_num_sites(unsigned subset) const {
        assert(subset < _num_subsets);
        return _subset_sizes[subset];
    }

    inline unsigned Model::get_num_sites() const {
        return _num_sites;
    }

    inline unsigned Model::get_num_subsets() const {
        return _num_subsets;
    }

    inline unsigned Model::get_subset_num_categ(unsigned subset) const {
        assert(subset < _num_subsets);
        assert(_asrv.size() == _num_subsets);
        assert(_asrv[subset]);
        return _asrv[subset]->get_num_categ();
    }

    inline bool Model::get_subset_is_invar_model(unsigned subset) const {
        assert(subset < _num_subsets);
        assert(_asrv.size() == _num_subsets);
        assert(_asrv[subset]);
        return _asrv[subset]->get_is_invar_model();
    }

    inline const QMatrix & Model::get_q_matrix(unsigned subset) const {
        assert(subset < _num_subsets);
        return *(_qmatrix[subset]);
    }

    inline const ASRV & Model::get_asrv(unsigned subset) const {
        assert(subset < _num_subsets);
        return *(_asrv[subset]);
    }

    inline double Model::calc_normalizing_constant_for_subset_rel_rates() const {
        // normalize _relrates so that expected relative rate across subsets equals 1.0
        double normalizing_constant = 0.0;
        for (unsigned s = 0; s < _num_subsets; s++) {
            normalizing_constant += _subset_sizes[s]*_subset_relrates[s]/_num_sites;
        }
        return normalizing_constant;
    }

    inline Model::subset_sizes_t & Model::get_subset_sizes() {
        return _subset_sizes;
    }

    inline void Model::set_subset_sizes(const subset_sizes_t nsites_vect) {
        assert(nsites_vect.size() == _num_subsets);
        _subset_sizes.resize(_num_subsets);
        std::copy(nsites_vect.begin(), nsites_vect.end(), _subset_sizes.begin());
        _num_sites = std::accumulate(_subset_sizes.begin(), _subset_sizes.end(), 0);
    }

    inline void Model::set_subset_num_patterns(const subset_sizes_t npatterns_vect) {
        assert(npatterns_vect.size() == _num_subsets);
        _subset_npatterns.resize(_num_subsets);
        std::copy(npatterns_vect.begin(), npatterns_vect.end(), _subset_npatterns.begin());
    }

    inline void Model::set_subset_data_types(const subset_datatype_t & datatype_vect) {
        _num_subsets = (unsigned)datatype_vect.size();

        _qmatrix.clear();
        _qmatrix.resize(_num_subsets);

        _asrv.clear();
        _asrv.resize(_num_subsets);

        _subset_datatypes.resize(_num_subsets);
        std::copy(datatype_vect.begin(), datatype_vect.end(), _subset_datatypes.begin());

		_subset_relrates.assign(_num_subsets, 1.0);

        for (unsigned s = 0; s < _num_subsets; s++) {
            _asrv[s].reset(new ASRV());
            if (_subset_datatypes[s].is_nucleotide())
                _qmatrix[s].reset(new QMatrixNucleotide());
            else if (_subset_datatypes[s].is_codon()) {
                GeneticCode::SharedPtr gcptr = _subset_datatypes[s].get_genetic_code();
                _qmatrix[s].reset(new QMatrixCodon(gcptr));
                }
            else
                throw NamuX(boost::format("Only nucleotide or codon data allowed in this version, you specified data type \"%s\" for subset %d") % _subset_datatypes[s].get_data_type_as_string() % (s+1));
        }
    }

    inline void Model::set_subset_num_categ(unsigned ncateg, unsigned subset) {
        assert(subset < _num_subsets);
        if (ncateg < 1) {
            throw NamuX(boost::str(boost::format("number of categories used for among-site rate variation must be greater than zero but the value %d was supplied") % ncateg));
        }
        _asrv[subset]->set_num_categ(ncateg);
    }

    inline void Model::set_subset_rate_var(ASRV::ratevar_ptr_t ratevar, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(ratevar);
        if (*ratevar < 0.0)
            throw NamuX(boost::str(boost::format("rate variance must be greater than or equal to zero but the value %.5f was supplied") % *ratevar));
        _asrv[subset]->set_rate_var_shared_ptr(ratevar);
        _asrv[subset]->fix_rate_var(fixed);
    }

    inline void Model::set_subset_pinvar(ASRV::pinvar_ptr_t pinvar, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(pinvar);
        if (*pinvar < 0.0)
            throw NamuX(boost::str(boost::format("proportion of invariable sites must be greater than or equal to zero but the value %.5f was supplied") % *pinvar));
        if (*pinvar >= 1.0)
            throw NamuX(boost::str(boost::format("proportion of invariable sites must be less than one but the value %.5f was supplied") % *pinvar));
        _asrv[subset]->set_pinvar_shared_ptr(pinvar);
        _asrv[subset]->fix_pinvar(fixed);
    }

    inline void Model::set_subset_is_invar_model(bool is_invar, unsigned subset) {
        assert(subset < _num_subsets);
        _asrv[subset]->set_is_invar_model(is_invar);
    }

    inline void Model::set_subset_exchangeabilities(QMatrix::freq_xchg_ptr_t exchangeabilities, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        if (!_subset_datatypes[subset].is_codon()) {
            double first_xchg = (*exchangeabilities)[0];
            if (first_xchg == -1)
                _qmatrix[subset]->set_equal_exchangeabilities(exchangeabilities);
            else
                _qmatrix[subset]->set_exchangeabilities_shared_ptr(exchangeabilities);
            _qmatrix[subset]->fix_exchangeabilities(fixed);
        }
    }

    inline void Model::set_subset_state_freqs(QMatrix::freq_xchg_ptr_t state_frequencies, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        double first_freq = (*state_frequencies)[0];
        if (first_freq == -1)
            _qmatrix[subset]->set_equal_state_freqs(state_frequencies);
        else
            _qmatrix[subset]->set_state_freqs_shared_ptr(state_frequencies);
        _qmatrix[subset]->fix_state_freqs(fixed);
    }

    inline void Model::set_subset_omega(QMatrix::omega_ptr_t omega, unsigned subset, bool fixed) {
        assert(subset < _num_subsets);
        assert(*omega > 0.0);
        if (_subset_datatypes[subset].is_codon()) {
            _qmatrix[subset]->set_omega_shared_ptr(omega);
            _qmatrix[subset]->fix_omega(fixed);
        }
    }

    inline void Model::activate() {
        for (auto q : _qmatrix)
            q->set_active(true);
    }

    inline void Model::inactivate() {
        for (auto q : _qmatrix)
            q->set_active(false);
    }

    inline int Model::set_beagle_eigen_decomposition(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _qmatrix.size());
        const double * pevec = _qmatrix[subset]->get_eigenvectors();
        const double * pivec = _qmatrix[subset]->get_inverse_eigenvectors();
        const double * pival = _qmatrix[subset]->get_eigenvalues();
        int code = beagleSetEigenDecomposition(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of eigen-decomposition buffer (input)
            pevec,              // Flattened matrix (stateCount x stateCount) of eigen-vectors (input)
            pivec,              // Flattened matrix (stateCount x stateCount) of inverse-eigen- vectors (input)
            pival);             // Vector of eigenvalues

        return code;
    }

    inline int Model::set_beagle_state_frequencies(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _qmatrix.size());
        const double * pfreq = _qmatrix[subset]->get_state_freqs();
        int code = beagleSetStateFrequencies(
             beagle_instance,   // Instance number (input)
             instance_subset,   // Index of state frequencies buffer (input)
             pfreq);            // State frequencies array (stateCount) (input)

        return code;
    }

    inline int Model::set_beagle_among_site_rate_variation_rates(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _asrv.size());
        const double * prates = _asrv[subset]->get_rates();
        int code = beagleSetCategoryRatesWithIndex(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of category rates buffer (input)
            prates);            // Array containing categoryCount rate scalers (input)

        return code;
    }

    inline int Model::set_beagle_among_site_rate_variation_probs(int beagle_instance, unsigned subset, unsigned instance_subset) {
        assert(subset < _asrv.size());
        const double * pprobs = _asrv[subset]->get_probs();
        int code = beagleSetCategoryWeights(
            beagle_instance,    // Instance number (input)
            instance_subset,    // Index of category weights buffer (input)
            pprobs);            // Category weights array (categoryCount) (input)

        return code;
    }

    inline void Model::set_subset_rel_rates(subset_relrate_vect_t & relrates, bool fixed) {
        assert(_num_subsets > 0);
        assert(relrates.size() > 0);
        if (relrates[0] == -1)
            _subset_relrates.assign(_num_subsets, 1.0);
        else
            _subset_relrates.assign(relrates.begin(), relrates.end());
        _subset_relrates_fixed = fixed;
    }

    inline Model::subset_relrate_vect_t & Model::get_subset_rel_rates() {
        return _subset_relrates;
    }

    inline bool Model::is_fixed_subset_rel_rates() const {
        return _subset_relrates_fixed;
    }

    inline void Model::set_tree_index(unsigned i, bool fixed) {
        _tree_index = i;
        _tree_fixed = fixed;
    }

    inline unsigned Model::get_tree_index() const {
        return _tree_index;
    }

    inline bool Model::is_fixed_tree() const {
        return _tree_fixed;
    }
}
