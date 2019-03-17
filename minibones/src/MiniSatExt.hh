/* 
 * File:   MiniSatExt.hh
 * Author: mikolas
 *
 * Created on November 29, 2010, 5:40 PM
 */

#ifndef MINISATEXT_HH
#define	MINISATEXT_HH
#include "minisat/core/Solver.h"

namespace Minisat {
  class MiniSatExt : public Solver {
  public:
    inline void bump(const Var var)        { varBumpActivity(var); }
    inline void new_variables(Var max_id) {
      const int target_number = (int)max_id+1;
      while (nVars() < target_number) newVar();
    }
  };
}
#endif	/* MINISATEXT_HH */

