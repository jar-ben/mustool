#ifndef CustomGlucose_Clone_h
#define CustomGlucose_Clone_h


namespace CustomGlucose {

    class Clone {
        public:
          virtual Clone* clone() const = 0;
    };
};

#endif
