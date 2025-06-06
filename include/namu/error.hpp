#pragma once

#include <boost/format.hpp>

namespace namu {

    /**
     * Class for errors.
     */
    class NamuX: public std::exception {
        public:
                            NamuX() throw() {}
                            NamuX(const std::string s) throw() : _msg() {_msg = s;}
                            NamuX(const boost::format & f) throw() : _msg() {_msg = boost::str(f);}
            virtual         ~NamuX() throw() {}
            const char *    what() const throw() { return _msg.c_str(); }

        private:

            std::string     _msg;
    };
}
