//bc:jpms,mj
/*----------------------------------------------------------------------------*\
 * File:        ToolConfig.hh
\*----------------------------------------------------------------------------*/
//ec:jpms,mj

#ifndef _TOOLCONFIG_H
#define _TOOLCONFIG_H 1
#include <assert.h>
#include <string>
using std::string;
using std::ostream;

static const char* const toolname     = "MiniBones";
static const char* const authorname   = "Mikolas Janota";
static const char* const authoremail  = "mikolas.janota@gmail.com";
static const char* const contribs     = "J. Marques-Silva & I. Lynce";

//jpms:bc
/*----------------------------------------------------------------------------*\
 * Types & Defines
\*----------------------------------------------------------------------------*/
//jpms:ec

typedef enum dc_modes {
  DC_None = 0x1000,
  DC_VL = 0x1001,
  DC_SC = 0x1002
} DCModes;

#define cfg_pref  config.get_prefix()

#define cout_pref cout << cfg_pref

#define report(x) cout << cfg_pref << x << endl;


//jpms:bc
/*----------------------------------------------------------------------------*\
 * Class: ToolConfig
 *
 * Purpose: Options configuring tool execution.
\*----------------------------------------------------------------------------*/
//jpms:ec

class ToolConfig {
public:

  ToolConfig()
  : _cmdstr()
  , _use_cores(false)
  , _use_core_based(false)
  , _use_upper_bound_prog(false)
  , _use_upper_bound(false)
  , _chunk_size(0)
  , _use_variable_bumping(0)
  , _use_chunk_keeping(0)
  , _use_random_chunks(0)
  , _backbone_insertion(0)
  , _self_test(false)
  , _dc_pruning(DC_None)
  , _rotatable_pruning(false)
  {
    _verbosity = 1; _timeout = 3600; _comp_fmt = true;
    _opt_level = 2; _runs_pvar = 1;  
    _output_prefix = "c "; _stats = false;
    set_run_mode(false, true); _incr_mode = true; _phase = 2;
    _implication_order = false;
  }

  inline ostream& prefix(ostream& output) const {output << (string) get_prefix(); return output;}
  
  inline void set_input_file_name (const char* value) { _input_file_name = value; }
  inline string get_input_file_name() { return _input_file_name;}

  int get_use_upper_bound() const { return _use_upper_bound; }
  void set_use_upper_bound(bool value = true) { _use_upper_bound = value; }

  int get_use_cores() const { return _use_cores; }
  void set_use_cores(bool value = true) { _use_cores = value; }

  int get_use_core_based() const { return _use_core_based; }
  void set_use_core_based(bool value = true) { _use_core_based = value; }


  int get_use_upper_bound_prog() const { return _use_upper_bound_prog; }
  void set_use_upper_bound_prog(bool value = true) { _use_upper_bound_prog = value; }

  int get_use_variable_bumping() const { return _use_variable_bumping; }
  void set_use_variable_bumping(int value) { _use_variable_bumping = value; }

  int get_use_chunk_keeping() const { return _use_chunk_keeping; }
  void set_use_chunk_keeping(int value) { _use_chunk_keeping = value; }

  int get_use_random_chunks() const { return _use_random_chunks; }
  void set_use_random_chunks(int value) { _use_random_chunks = value; }


  int get_chunk_size() const { return _chunk_size; }
  void set_chunk_size(int chunk_size) { _chunk_size = chunk_size; }

  int get_backbone_insertion() const { return _backbone_insertion; }
  void set_backbone_insertion(int backbone_insertion) { _backbone_insertion = backbone_insertion; }


  int get_self_test() const { return _self_test; }
  void set_self_test(bool value = true) { _self_test = value; }

  int get_rotatable_pruning() const { return _rotatable_pruning; }
  void set_rotatable_pruning(bool value = true) { _rotatable_pruning = value; }


  int get_verbosity()  const  { return _verbosity; }

  void set_verbosity(int verb) { _verbosity = verb; }

  int get_timeout() { return _timeout; }

  void set_timeout(int tout) { _timeout = tout; }

  const char* get_prefix() const { return _output_prefix; }

  bool get_comp_format() { return _comp_fmt; }

  void set_comp_format() { _comp_fmt = true; _output_prefix = "c "; }

  void unset_comp_format() { _comp_fmt = false; _output_prefix = ""; }

  bool get_stats() { return _stats; }

  void set_stats() { _stats = true; }

  void unset_stats() { _stats = false; }

  bool get_run_enum_mode() { return _enum_mode; }

  void set_run_enum_mode() { set_run_mode(true, false); }

  void unset_run_enum_mode() { set_run_mode(false, true); }

  bool get_run_iter_mode() { return _iter_mode; }

  void set_run_iter_mode() { set_run_mode(false, true); }

  void unset_run_iter_mode() { set_run_mode(true, false); }

  bool get_incr_iter_mode() { assert(_iter_mode); return _incr_mode; }

  void set_incr_iter_mode() { _incr_mode = true; }

  bool get_implication_order() {return _implication_order;}
  void set_implication_order(bool value = true) {_implication_order = value;}
  
  void unset_incr_iter_mode() { _incr_mode = false; }

  int get_optimize() { return _opt_level; }

  void set_optimize(int optl) { _opt_level = optl; }

  int get_runs_per_var() { return _runs_pvar; }

  void set_runs_per_var(int nr) { _runs_pvar = nr; }

  bool get_nodc_pruning() const { return _dc_pruning == DC_None; }

  bool get_vldc_pruning() { return _dc_pruning == DC_VL; }

  bool get_scdc_pruning() const { return _dc_pruning == DC_SC; }

  void set_nodc_pruning() { _dc_pruning = DC_None; }

  void set_vldc_pruning() { _dc_pruning = DC_VL; }

  void set_scdc_pruning() { _dc_pruning = DC_SC; }

  int get_phase() { return _phase; }

  void set_phase(int nph) { _phase = nph; }
protected:

  void set_run_mode(bool en_mode, bool it_mode) {
    assert(en_mode && !it_mode || !en_mode && it_mode);
    _enum_mode = en_mode; _iter_mode = it_mode;
  }

protected:
  string _cmdstr;

  string _input_file_name;

  bool _use_cores;

  bool _use_core_based;

  bool _use_upper_bound_prog;
      
  bool _use_upper_bound;

  int _chunk_size;

  int _use_variable_bumping;

  int _use_chunk_keeping;

  int _use_random_chunks;

  int _backbone_insertion;

  bool _self_test;

  int _verbosity;

  int _timeout;

  const char* _output_prefix;

  bool _comp_fmt;

  bool _stats;

  bool _enum_mode;

  bool _iter_mode;

  bool _incr_mode;
    
  bool _implication_order;

  int _opt_level;

  int _runs_pvar;

  DCModes _dc_pruning;

  int _phase;

  bool _rotatable_pruning;
};

#endif /* _TOOLCONFIG_H */

/*----------------------------------------------------------------------------*/
