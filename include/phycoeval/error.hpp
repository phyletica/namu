#pragma once

#include <boost/format.hpp>

namespace phycoeval {

    /**
     * Class for errors.
     */
    class PhycoevalError: public std::exception {
        public:
                            PhycoevalError() throw() {}
                            PhycoevalError(const std::string s) throw() : _msg() {_msg = s;}
                            PhycoevalError(const boost::format & f) throw() : _msg() {_msg = boost::str(f);}
            virtual         ~PhycoevalError() throw() {}
            const char *    what() const throw() { return _msg.c_str(); }

        private:

            std::string     _msg;
    };
}
