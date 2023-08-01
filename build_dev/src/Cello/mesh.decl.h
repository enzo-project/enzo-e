#ifndef _DECL_mesh_H_
#define _DECL_mesh_H_
#include "charm++.h"
#include "envelope.h"
#include <memory>
#include "sdag.h"









/* DECLS: readonly int MsgCoarsen::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int MsgAdapt::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int MsgInitial::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int MsgOutput::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int MsgRefine::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int MsgRefresh::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int DataMsg::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int FieldFace::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int ParticleData::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int InitialTrace::id0_[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly double Method::courant_global;
 */

/* DECLS: readonly Config g_config;
 */

/* DECLS: readonly Parameters g_parameters;
 */
















































#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message FieldMsg{
char a[];
}
;
 */
class FieldMsg;
class CMessage_FieldMsg:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_FieldMsg();
    static void *pack(FieldMsg *p);
    static FieldMsg* unpack(void* p);
    void *operator new(size_t, int);
    void *operator new(size_t, int, const int);
    void *operator new(size_t, int, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, int){dealloc(p);}
    void operator delete(void *p, int, const int){dealloc(p);}
    void operator delete(void *p, int, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgAdapt;
 */
class MsgAdapt;
class CMessage_MsgAdapt:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgAdapt();
    static void *pack(MsgAdapt *p);
    static MsgAdapt* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgCoarsen;
 */
class MsgCoarsen;
class CMessage_MsgCoarsen:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgCoarsen();
    static void *pack(MsgCoarsen *p);
    static MsgCoarsen* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgInitial;
 */
class MsgInitial;
class CMessage_MsgInitial:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgInitial();
    static void *pack(MsgInitial *p);
    static MsgInitial* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgOutput;
 */
class MsgOutput;
class CMessage_MsgOutput:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgOutput();
    static void *pack(MsgOutput *p);
    static MsgOutput* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgRefine;
 */
class MsgRefine;
class CMessage_MsgRefine:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgRefine();
    static void *pack(MsgRefine *p);
    static MsgRefine* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

#ifndef GROUPDEPNUM_DECLARED
# define GROUPDEPNUM_DECLARED
struct GroupDepNum
{
  int groupDepNum;
  explicit GroupDepNum(int g = 0) : groupDepNum{g} { }
  operator int() const { return groupDepNum; }
};
#endif
/* DECLS: message MsgRefresh;
 */
class MsgRefresh;
class CMessage_MsgRefresh:public CkMessage{
  public:
    static int __idx;
    void* operator new(size_t, void*p) { return p; }
    void* operator new(size_t);
    void* operator new(size_t, int*, const int);
    void* operator new(size_t, int*, const int, const GroupDepNum);
    void* operator new(size_t, int*);
#if CMK_MULTIPLE_DELETE
    void operator delete(void*p, void*){dealloc(p);}
    void operator delete(void*p){dealloc(p);}
    void operator delete(void*p, int*, const int){dealloc(p);}
    void operator delete(void*p, int*, const int, const GroupDepNum){dealloc(p);}
    void operator delete(void*p, int*){dealloc(p);}
#endif
    void operator delete(void*p, size_t){dealloc(p);}
    static void* alloc(int,size_t, int*, int, GroupDepNum);
    static void dealloc(void *p);
    CMessage_MsgRefresh();
    static void *pack(MsgRefresh *p);
    static MsgRefresh* unpack(void* p);
    void *operator new(size_t, const int);
    void *operator new(size_t, const int, const GroupDepNum);
#if CMK_MULTIPLE_DELETE
    void operator delete(void *p, const int){dealloc(p);}
    void operator delete(void *p, const int, const GroupDepNum){dealloc(p);}
#endif
    static void __register(const char *s, size_t size, CkPackFnPtr pack, CkUnpackFnPtr unpack) {
      __idx = CkRegisterMsg(s, pack, unpack, dealloc, size);
    }
};

/* DECLS: array Block: ArrayElement{
Block(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
Block();
void p_initial_exit();
void r_end_initialize(CkReductionMsg* impl_msg);
void r_initial_new_next(CkReductionMsg* impl_msg);
void r_initial_new_continue(CkReductionMsg* impl_msg);
void p_output_enter();
void r_output_enter(CkReductionMsg* impl_msg);
void p_output_end();
void p_output_exit();
void r_output_exit(CkReductionMsg* impl_msg);
void p_output_write(int index_output, int step);
void p_compute_enter();
void p_compute_continue();
void r_compute_continue(CkReductionMsg* impl_msg);
void p_compute_exit();
void r_compute_exit(CkReductionMsg* impl_msg);
void p_method_flux_correct_refresh();
void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
void r_method_order_morton_continue(CkReductionMsg* impl_msg);
void r_method_order_morton_complete(CkReductionMsg* impl_msg);
void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
void p_method_order_morton_index(int index, int count);
void p_method_output_next(MsgOutput* impl_msg);
void p_method_output_write(MsgOutput* impl_msg);
void r_method_output_continue(CkReductionMsg* impl_msg);
void r_method_output_done(CkReductionMsg* impl_msg);
void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
void p_stopping_enter();
void r_stopping_enter(CkReductionMsg* impl_msg);
void p_stopping_load_balance();
void r_stopping_load_balance(CkReductionMsg* impl_msg);
void p_stopping_exit();
void r_stopping_exit(CkReductionMsg* impl_msg);
void p_exit();
void r_exit(CkReductionMsg* impl_msg);
void p_control_sync_count(int entry_point, int id, int count);
void p_adapt_enter();
void r_adapt_enter(CkReductionMsg* impl_msg);
void p_adapt_end();
void p_adapt_update();
void r_adapt_next(CkReductionMsg* impl_msg);
void p_adapt_called();
void p_adapt_exit();
void p_adapt_delete();
void p_adapt_recv_level(MsgAdapt* impl_msg);
void p_adapt_recv_child(MsgCoarsen* impl_msg);
void r_restart_enter(CkReductionMsg* impl_msg);
void p_refresh_recv(MsgRefresh* impl_msg);
void p_refresh_child(int n, const char *a, const int *ic3);
Block(CkMigrateMessage* impl_msg);
};
 */
 class Block;
 class CkIndex_Block;
 class CProxy_Block;
 class CProxyElement_Block;
 class CProxySection_Block;
/* --------------- index object ------------------ */
class CkIndex_Block:public CkIndex_ArrayElement{
  public:
    typedef Block local_t;
    typedef CkIndex_Block index_t;
    typedef CProxy_Block proxy_t;
    typedef CProxyElement_Block element_t;
    typedef CProxySection_Block section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: Block(const process_type &ip_source, const MsgType &msg_type);
     */
    // Entry point registration at startup
    
    static int reg_Block_marshall1();
    // Entry point index lookup
    
    inline static int idx_Block_marshall1() {
      static int epidx = reg_Block_marshall1();
      return epidx;
    }

    
    static int ckNew(const process_type &ip_source, const MsgType &msg_type) { return idx_Block_marshall1(); }
    
    static void _call_Block_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Block_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_Block_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_Block_marshall1(PUP::er &p,void *msg);
    /* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_set_msg_refine_MsgRefine();
    // Entry point index lookup
    
    inline static int idx_p_set_msg_refine_MsgRefine() {
      static int epidx = reg_p_set_msg_refine_MsgRefine();
      return epidx;
    }

    
    inline static int idx_p_set_msg_refine(void (Block::*)(MsgRefine* impl_msg) ) {
      return idx_p_set_msg_refine_MsgRefine();
    }


    
    static int p_set_msg_refine(MsgRefine* impl_msg) { return idx_p_set_msg_refine_MsgRefine(); }
    
    static void _call_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj);
    /* DECLS: Block();
     */
    // Entry point registration at startup
    
    static int reg_Block_void();
    // Entry point index lookup
    
    inline static int idx_Block_void() {
      static int epidx = reg_Block_void();
      return epidx;
    }

    
    static int ckNew() { return idx_Block_void(); }
    
    static void _call_Block_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Block_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_initial_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_initial_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_initial_exit_void() {
      static int epidx = reg_p_initial_exit_void();
      return epidx;
    }

    
    inline static int idx_p_initial_exit(void (Block::*)() ) {
      return idx_p_initial_exit_void();
    }


    
    static int p_initial_exit() { return idx_p_initial_exit_void(); }
    
    static void _call_p_initial_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_initial_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_end_initialize(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_end_initialize_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_end_initialize_CkReductionMsg() {
      static int epidx = reg_r_end_initialize_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_end_initialize(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_end_initialize_CkReductionMsg();
    }


    
    static int r_end_initialize(CkReductionMsg* impl_msg) { return idx_r_end_initialize_CkReductionMsg(); }
    
    static void _call_r_end_initialize_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_end_initialize_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_initial_new_next(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_initial_new_next_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_initial_new_next_CkReductionMsg() {
      static int epidx = reg_r_initial_new_next_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_initial_new_next(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_initial_new_next_CkReductionMsg();
    }


    
    static int r_initial_new_next(CkReductionMsg* impl_msg) { return idx_r_initial_new_next_CkReductionMsg(); }
    
    static void _call_r_initial_new_next_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_initial_new_next_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_initial_new_continue(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_initial_new_continue_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_initial_new_continue_CkReductionMsg() {
      static int epidx = reg_r_initial_new_continue_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_initial_new_continue(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_initial_new_continue_CkReductionMsg();
    }


    
    static int r_initial_new_continue(CkReductionMsg* impl_msg) { return idx_r_initial_new_continue_CkReductionMsg(); }
    
    static void _call_r_initial_new_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_initial_new_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_output_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_output_enter_void() {
      static int epidx = reg_p_output_enter_void();
      return epidx;
    }

    
    inline static int idx_p_output_enter(void (Block::*)() ) {
      return idx_p_output_enter_void();
    }


    
    static int p_output_enter() { return idx_p_output_enter_void(); }
    
    static void _call_p_output_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_output_enter(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_output_enter_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_output_enter_CkReductionMsg() {
      static int epidx = reg_r_output_enter_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_output_enter(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_output_enter_CkReductionMsg();
    }


    
    static int r_output_enter(CkReductionMsg* impl_msg) { return idx_r_output_enter_CkReductionMsg(); }
    
    static void _call_r_output_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_output_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_end();
     */
    // Entry point registration at startup
    
    static int reg_p_output_end_void();
    // Entry point index lookup
    
    inline static int idx_p_output_end_void() {
      static int epidx = reg_p_output_end_void();
      return epidx;
    }

    
    inline static int idx_p_output_end(void (Block::*)() ) {
      return idx_p_output_end_void();
    }


    
    static int p_output_end() { return idx_p_output_end_void(); }
    
    static void _call_p_output_end_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_end_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_output_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_output_exit_void() {
      static int epidx = reg_p_output_exit_void();
      return epidx;
    }

    
    inline static int idx_p_output_exit(void (Block::*)() ) {
      return idx_p_output_exit_void();
    }


    
    static int p_output_exit() { return idx_p_output_exit_void(); }
    
    static void _call_p_output_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_output_exit(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_output_exit_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_output_exit_CkReductionMsg() {
      static int epidx = reg_r_output_exit_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_output_exit(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_output_exit_CkReductionMsg();
    }


    
    static int r_output_exit(CkReductionMsg* impl_msg) { return idx_r_output_exit_CkReductionMsg(); }
    
    static void _call_r_output_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_output_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_write(int index_output, int step);
     */
    // Entry point registration at startup
    
    static int reg_p_output_write_marshall13();
    // Entry point index lookup
    
    inline static int idx_p_output_write_marshall13() {
      static int epidx = reg_p_output_write_marshall13();
      return epidx;
    }

    
    inline static int idx_p_output_write(void (Block::*)(int index_output, int step) ) {
      return idx_p_output_write_marshall13();
    }


    
    static int p_output_write(int index_output, int step) { return idx_p_output_write_marshall13(); }
    
    static void _call_p_output_write_marshall13(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_write_marshall13(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_output_write_marshall13(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_output_write_marshall13(PUP::er &p,void *msg);
    /* DECLS: void p_compute_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_compute_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_compute_enter_void() {
      static int epidx = reg_p_compute_enter_void();
      return epidx;
    }

    
    inline static int idx_p_compute_enter(void (Block::*)() ) {
      return idx_p_compute_enter_void();
    }


    
    static int p_compute_enter() { return idx_p_compute_enter_void(); }
    
    static void _call_p_compute_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_compute_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_compute_continue();
     */
    // Entry point registration at startup
    
    static int reg_p_compute_continue_void();
    // Entry point index lookup
    
    inline static int idx_p_compute_continue_void() {
      static int epidx = reg_p_compute_continue_void();
      return epidx;
    }

    
    inline static int idx_p_compute_continue(void (Block::*)() ) {
      return idx_p_compute_continue_void();
    }


    
    static int p_compute_continue() { return idx_p_compute_continue_void(); }
    
    static void _call_p_compute_continue_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_compute_continue_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_compute_continue(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_compute_continue_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_compute_continue_CkReductionMsg() {
      static int epidx = reg_r_compute_continue_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_compute_continue(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_compute_continue_CkReductionMsg();
    }


    
    static int r_compute_continue(CkReductionMsg* impl_msg) { return idx_r_compute_continue_CkReductionMsg(); }
    
    static void _call_r_compute_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_compute_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_compute_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_compute_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_compute_exit_void() {
      static int epidx = reg_p_compute_exit_void();
      return epidx;
    }

    
    inline static int idx_p_compute_exit(void (Block::*)() ) {
      return idx_p_compute_exit_void();
    }


    
    static int p_compute_exit() { return idx_p_compute_exit_void(); }
    
    static void _call_p_compute_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_compute_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_compute_exit(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_compute_exit_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_compute_exit_CkReductionMsg() {
      static int epidx = reg_r_compute_exit_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_compute_exit(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_compute_exit_CkReductionMsg();
    }


    
    static int r_compute_exit(CkReductionMsg* impl_msg) { return idx_r_compute_exit_CkReductionMsg(); }
    
    static void _call_r_compute_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_compute_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_flux_correct_refresh();
     */
    // Entry point registration at startup
    
    static int reg_p_method_flux_correct_refresh_void();
    // Entry point index lookup
    
    inline static int idx_p_method_flux_correct_refresh_void() {
      static int epidx = reg_p_method_flux_correct_refresh_void();
      return epidx;
    }

    
    inline static int idx_p_method_flux_correct_refresh(void (Block::*)() ) {
      return idx_p_method_flux_correct_refresh_void();
    }


    
    static int p_method_flux_correct_refresh() { return idx_p_method_flux_correct_refresh_void(); }
    
    static void _call_p_method_flux_correct_refresh_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_flux_correct_refresh_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_flux_correct_sum_fields_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_flux_correct_sum_fields_CkReductionMsg() {
      static int epidx = reg_r_method_flux_correct_sum_fields_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_flux_correct_sum_fields(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_flux_correct_sum_fields_CkReductionMsg();
    }


    
    static int r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) { return idx_r_method_flux_correct_sum_fields_CkReductionMsg(); }
    
    static void _call_r_method_flux_correct_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_flux_correct_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_order_morton_continue_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_order_morton_continue_CkReductionMsg() {
      static int epidx = reg_r_method_order_morton_continue_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_order_morton_continue(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_order_morton_continue_CkReductionMsg();
    }


    
    static int r_method_order_morton_continue(CkReductionMsg* impl_msg) { return idx_r_method_order_morton_continue_CkReductionMsg(); }
    
    static void _call_r_method_order_morton_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_order_morton_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_order_morton_complete_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_order_morton_complete_CkReductionMsg() {
      static int epidx = reg_r_method_order_morton_complete_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_order_morton_complete(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_order_morton_complete_CkReductionMsg();
    }


    
    static int r_method_order_morton_complete(CkReductionMsg* impl_msg) { return idx_r_method_order_morton_complete_CkReductionMsg(); }
    
    static void _call_r_method_order_morton_complete_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_order_morton_complete_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
     */
    // Entry point registration at startup
    
    static int reg_p_method_order_morton_weight_marshall23();
    // Entry point index lookup
    
    inline static int idx_p_method_order_morton_weight_marshall23() {
      static int epidx = reg_p_method_order_morton_weight_marshall23();
      return epidx;
    }

    
    inline static int idx_p_method_order_morton_weight(void (Block::*)(const int *ic3, int weight, const Index &index) ) {
      return idx_p_method_order_morton_weight_marshall23();
    }


    
    static int p_method_order_morton_weight(const int *ic3, int weight, const Index &index) { return idx_p_method_order_morton_weight_marshall23(); }
    
    static void _call_p_method_order_morton_weight_marshall23(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_order_morton_weight_marshall23(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_order_morton_weight_marshall23(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_order_morton_weight_marshall23(PUP::er &p,void *msg);
    /* DECLS: void p_method_order_morton_index(int index, int count);
     */
    // Entry point registration at startup
    
    static int reg_p_method_order_morton_index_marshall24();
    // Entry point index lookup
    
    inline static int idx_p_method_order_morton_index_marshall24() {
      static int epidx = reg_p_method_order_morton_index_marshall24();
      return epidx;
    }

    
    inline static int idx_p_method_order_morton_index(void (Block::*)(int index, int count) ) {
      return idx_p_method_order_morton_index_marshall24();
    }


    
    static int p_method_order_morton_index(int index, int count) { return idx_p_method_order_morton_index_marshall24(); }
    
    static void _call_p_method_order_morton_index_marshall24(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_order_morton_index_marshall24(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_order_morton_index_marshall24(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_order_morton_index_marshall24(PUP::er &p,void *msg);
    /* DECLS: void p_method_output_next(MsgOutput* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_method_output_next_MsgOutput();
    // Entry point index lookup
    
    inline static int idx_p_method_output_next_MsgOutput() {
      static int epidx = reg_p_method_output_next_MsgOutput();
      return epidx;
    }

    
    inline static int idx_p_method_output_next(void (Block::*)(MsgOutput* impl_msg) ) {
      return idx_p_method_output_next_MsgOutput();
    }


    
    static int p_method_output_next(MsgOutput* impl_msg) { return idx_p_method_output_next_MsgOutput(); }
    
    static void _call_p_method_output_next_MsgOutput(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_output_next_MsgOutput(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_output_write(MsgOutput* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_method_output_write_MsgOutput();
    // Entry point index lookup
    
    inline static int idx_p_method_output_write_MsgOutput() {
      static int epidx = reg_p_method_output_write_MsgOutput();
      return epidx;
    }

    
    inline static int idx_p_method_output_write(void (Block::*)(MsgOutput* impl_msg) ) {
      return idx_p_method_output_write_MsgOutput();
    }


    
    static int p_method_output_write(MsgOutput* impl_msg) { return idx_p_method_output_write_MsgOutput(); }
    
    static void _call_p_method_output_write_MsgOutput(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_output_write_MsgOutput(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_output_continue(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_output_continue_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_output_continue_CkReductionMsg() {
      static int epidx = reg_r_method_output_continue_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_output_continue(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_output_continue_CkReductionMsg();
    }


    
    static int r_method_output_continue(CkReductionMsg* impl_msg) { return idx_r_method_output_continue_CkReductionMsg(); }
    
    static void _call_r_method_output_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_output_continue_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_output_done(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_output_done_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_output_done_CkReductionMsg() {
      static int epidx = reg_r_method_output_done_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_output_done(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_output_done_CkReductionMsg();
    }


    
    static int r_method_output_done(CkReductionMsg* impl_msg) { return idx_r_method_output_done_CkReductionMsg(); }
    
    static void _call_r_method_output_done_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_output_done_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_debug_sum_fields_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_debug_sum_fields_CkReductionMsg() {
      static int epidx = reg_r_method_debug_sum_fields_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_debug_sum_fields(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_debug_sum_fields_CkReductionMsg();
    }


    
    static int r_method_debug_sum_fields(CkReductionMsg* impl_msg) { return idx_r_method_debug_sum_fields_CkReductionMsg(); }
    
    static void _call_r_method_debug_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_debug_sum_fields_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_stopping_compute_timestep_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_stopping_compute_timestep_CkReductionMsg() {
      static int epidx = reg_r_stopping_compute_timestep_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_stopping_compute_timestep(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_stopping_compute_timestep_CkReductionMsg();
    }


    
    static int r_stopping_compute_timestep(CkReductionMsg* impl_msg) { return idx_r_stopping_compute_timestep_CkReductionMsg(); }
    
    static void _call_r_stopping_compute_timestep_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_stopping_compute_timestep_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_enter_void() {
      static int epidx = reg_p_stopping_enter_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_enter(void (Block::*)() ) {
      return idx_p_stopping_enter_void();
    }


    
    static int p_stopping_enter() { return idx_p_stopping_enter_void(); }
    
    static void _call_p_stopping_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_stopping_enter(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_stopping_enter_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_stopping_enter_CkReductionMsg() {
      static int epidx = reg_r_stopping_enter_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_stopping_enter(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_stopping_enter_CkReductionMsg();
    }


    
    static int r_stopping_enter(CkReductionMsg* impl_msg) { return idx_r_stopping_enter_CkReductionMsg(); }
    
    static void _call_r_stopping_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_stopping_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_load_balance();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_load_balance_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_load_balance_void() {
      static int epidx = reg_p_stopping_load_balance_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_load_balance(void (Block::*)() ) {
      return idx_p_stopping_load_balance_void();
    }


    
    static int p_stopping_load_balance() { return idx_p_stopping_load_balance_void(); }
    
    static void _call_p_stopping_load_balance_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_load_balance_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_stopping_load_balance_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_stopping_load_balance_CkReductionMsg() {
      static int epidx = reg_r_stopping_load_balance_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_stopping_load_balance(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_stopping_load_balance_CkReductionMsg();
    }


    
    static int r_stopping_load_balance(CkReductionMsg* impl_msg) { return idx_r_stopping_load_balance_CkReductionMsg(); }
    
    static void _call_r_stopping_load_balance_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_stopping_load_balance_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_exit_void() {
      static int epidx = reg_p_stopping_exit_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_exit(void (Block::*)() ) {
      return idx_p_stopping_exit_void();
    }


    
    static int p_stopping_exit() { return idx_p_stopping_exit_void(); }
    
    static void _call_p_stopping_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_stopping_exit(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_stopping_exit_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_stopping_exit_CkReductionMsg() {
      static int epidx = reg_r_stopping_exit_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_stopping_exit(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_stopping_exit_CkReductionMsg();
    }


    
    static int r_stopping_exit(CkReductionMsg* impl_msg) { return idx_r_stopping_exit_CkReductionMsg(); }
    
    static void _call_r_stopping_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_stopping_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_exit_void() {
      static int epidx = reg_p_exit_void();
      return epidx;
    }

    
    inline static int idx_p_exit(void (Block::*)() ) {
      return idx_p_exit_void();
    }


    
    static int p_exit() { return idx_p_exit_void(); }
    
    static void _call_p_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_exit(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_exit_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_exit_CkReductionMsg() {
      static int epidx = reg_r_exit_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_exit(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_exit_CkReductionMsg();
    }


    
    static int r_exit(CkReductionMsg* impl_msg) { return idx_r_exit_CkReductionMsg(); }
    
    static void _call_r_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_exit_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_control_sync_count(int entry_point, int id, int count);
     */
    // Entry point registration at startup
    
    static int reg_p_control_sync_count_marshall39();
    // Entry point index lookup
    
    inline static int idx_p_control_sync_count_marshall39() {
      static int epidx = reg_p_control_sync_count_marshall39();
      return epidx;
    }

    
    inline static int idx_p_control_sync_count(void (Block::*)(int entry_point, int id, int count) ) {
      return idx_p_control_sync_count_marshall39();
    }


    
    static int p_control_sync_count(int entry_point, int id, int count) { return idx_p_control_sync_count_marshall39(); }
    
    static void _call_p_control_sync_count_marshall39(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_control_sync_count_marshall39(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_control_sync_count_marshall39(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_control_sync_count_marshall39(PUP::er &p,void *msg);
    /* DECLS: void p_adapt_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_enter_void() {
      static int epidx = reg_p_adapt_enter_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_enter(void (Block::*)() ) {
      return idx_p_adapt_enter_void();
    }


    
    static int p_adapt_enter() { return idx_p_adapt_enter_void(); }
    
    static void _call_p_adapt_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_adapt_enter(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_adapt_enter_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_adapt_enter_CkReductionMsg() {
      static int epidx = reg_r_adapt_enter_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_adapt_enter(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_adapt_enter_CkReductionMsg();
    }


    
    static int r_adapt_enter(CkReductionMsg* impl_msg) { return idx_r_adapt_enter_CkReductionMsg(); }
    
    static void _call_r_adapt_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_adapt_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_end();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_end_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_end_void() {
      static int epidx = reg_p_adapt_end_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_end(void (Block::*)() ) {
      return idx_p_adapt_end_void();
    }


    
    static int p_adapt_end() { return idx_p_adapt_end_void(); }
    
    static void _call_p_adapt_end_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_end_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_update();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_update_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_update_void() {
      static int epidx = reg_p_adapt_update_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_update(void (Block::*)() ) {
      return idx_p_adapt_update_void();
    }


    
    static int p_adapt_update() { return idx_p_adapt_update_void(); }
    
    static void _call_p_adapt_update_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_update_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_adapt_next(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_adapt_next_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_adapt_next_CkReductionMsg() {
      static int epidx = reg_r_adapt_next_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_adapt_next(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_adapt_next_CkReductionMsg();
    }


    
    static int r_adapt_next(CkReductionMsg* impl_msg) { return idx_r_adapt_next_CkReductionMsg(); }
    
    static void _call_r_adapt_next_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_adapt_next_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_called();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_called_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_called_void() {
      static int epidx = reg_p_adapt_called_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_called(void (Block::*)() ) {
      return idx_p_adapt_called_void();
    }


    
    static int p_adapt_called() { return idx_p_adapt_called_void(); }
    
    static void _call_p_adapt_called_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_called_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_exit_void() {
      static int epidx = reg_p_adapt_exit_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_exit(void (Block::*)() ) {
      return idx_p_adapt_exit_void();
    }


    
    static int p_adapt_exit() { return idx_p_adapt_exit_void(); }
    
    static void _call_p_adapt_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_delete();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_delete_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_delete_void() {
      static int epidx = reg_p_adapt_delete_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_delete(void (Block::*)() ) {
      return idx_p_adapt_delete_void();
    }


    
    static int p_adapt_delete() { return idx_p_adapt_delete_void(); }
    
    static void _call_p_adapt_delete_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_delete_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_recv_level(MsgAdapt* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_recv_level_MsgAdapt();
    // Entry point index lookup
    
    inline static int idx_p_adapt_recv_level_MsgAdapt() {
      static int epidx = reg_p_adapt_recv_level_MsgAdapt();
      return epidx;
    }

    
    inline static int idx_p_adapt_recv_level(void (Block::*)(MsgAdapt* impl_msg) ) {
      return idx_p_adapt_recv_level_MsgAdapt();
    }


    
    static int p_adapt_recv_level(MsgAdapt* impl_msg) { return idx_p_adapt_recv_level_MsgAdapt(); }
    
    static void _call_p_adapt_recv_level_MsgAdapt(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_recv_level_MsgAdapt(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_recv_child_MsgCoarsen();
    // Entry point index lookup
    
    inline static int idx_p_adapt_recv_child_MsgCoarsen() {
      static int epidx = reg_p_adapt_recv_child_MsgCoarsen();
      return epidx;
    }

    
    inline static int idx_p_adapt_recv_child(void (Block::*)(MsgCoarsen* impl_msg) ) {
      return idx_p_adapt_recv_child_MsgCoarsen();
    }


    
    static int p_adapt_recv_child(MsgCoarsen* impl_msg) { return idx_p_adapt_recv_child_MsgCoarsen(); }
    
    static void _call_p_adapt_recv_child_MsgCoarsen(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_recv_child_MsgCoarsen(void* impl_msg, void* impl_obj);
    /* DECLS: void r_restart_enter(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_restart_enter_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_restart_enter_CkReductionMsg() {
      static int epidx = reg_r_restart_enter_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_restart_enter(void (Block::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_restart_enter_CkReductionMsg();
    }


    
    static int r_restart_enter(CkReductionMsg* impl_msg) { return idx_r_restart_enter_CkReductionMsg(); }
    
    static void _call_r_restart_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_restart_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_refresh_recv(MsgRefresh* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_refresh_recv_MsgRefresh();
    // Entry point index lookup
    
    inline static int idx_p_refresh_recv_MsgRefresh() {
      static int epidx = reg_p_refresh_recv_MsgRefresh();
      return epidx;
    }

    
    inline static int idx_p_refresh_recv(void (Block::*)(MsgRefresh* impl_msg) ) {
      return idx_p_refresh_recv_MsgRefresh();
    }


    
    static int p_refresh_recv(MsgRefresh* impl_msg) { return idx_p_refresh_recv_MsgRefresh(); }
    
    static void _call_p_refresh_recv_MsgRefresh(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_refresh_recv_MsgRefresh(void* impl_msg, void* impl_obj);
    /* DECLS: void p_refresh_child(int n, const char *a, const int *ic3);
     */
    // Entry point registration at startup
    
    static int reg_p_refresh_child_marshall52();
    // Entry point index lookup
    
    inline static int idx_p_refresh_child_marshall52() {
      static int epidx = reg_p_refresh_child_marshall52();
      return epidx;
    }

    
    inline static int idx_p_refresh_child(void (Block::*)(int n, const char *a, const int *ic3) ) {
      return idx_p_refresh_child_marshall52();
    }


    
    static int p_refresh_child(int n, const char *a, const int *ic3) { return idx_p_refresh_child_marshall52(); }
    
    static void _call_p_refresh_child_marshall52(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_refresh_child_marshall52(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_refresh_child_marshall52(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_refresh_child_marshall52(PUP::er &p,void *msg);
    /* DECLS: Block(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_Block_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_Block_CkMigrateMessage() {
      static int epidx = reg_Block_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_Block_CkMigrateMessage(); }
    
    static void _call_Block_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Block_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_Block : public CProxyElement_ArrayElement{
  public:
    typedef Block local_t;
    typedef CkIndex_Block index_t;
    typedef CProxy_Block proxy_t;
    typedef CProxyElement_Block element_t;
    typedef CProxySection_Block section_t;

    using array_index_t = CkArrayIndexIndex;

    /* TRAM aggregators */

    CProxyElement_Block(void) {
    }
    CProxyElement_Block(const ArrayElement *e) : CProxyElement_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_ArrayElement::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_ArrayElement::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_ArrayElement::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_ArrayElement::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_ArrayElement::ckGetIndex(); }

    Block *ckLocal(void) const
    { return (Block *)CProxyElement_ArrayElement::ckLocal(); }

    CProxyElement_Block(const CkArrayID &aid,const CkArrayIndexIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_Block(const CkArrayID &aid,const CkArrayIndexIndex &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}

    CProxyElement_Block(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_Block(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}
/* DECLS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
    
    void insert(const process_type &ip_source, const MsgType &msg_type, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: Block();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_initial_exit();
 */
    
    void p_initial_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
    
    void r_end_initialize(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_next(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_enter();
 */
    
    void p_output_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_enter(CkReductionMsg* impl_msg);
 */
    
    void r_output_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_end();
 */
    
    void p_output_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_output_exit();
 */
    
    void p_output_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_exit(CkReductionMsg* impl_msg);
 */
    
    void r_output_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_write(int index_output, int step);
 */
    
    void p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_enter();
 */
    
    void p_compute_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_continue();
 */
    
    void p_compute_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
    
    void r_compute_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_compute_exit();
 */
    
    void p_compute_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
    
    void r_compute_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_flux_correct_refresh();
 */
    
    void p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_complete(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
    
    void p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_order_morton_index(int index, int count);
 */
    
    void p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_output_next(MsgOutput* impl_msg);
 */
    
    void p_method_output_next(MsgOutput* impl_msg) ;

/* DECLS: void p_method_output_write(MsgOutput* impl_msg);
 */
    
    void p_method_output_write(MsgOutput* impl_msg) ;

/* DECLS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_done(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_debug_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_compute_timestep(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_enter();
 */
    
    void p_stopping_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_load_balance();
 */
    
    void p_stopping_load_balance(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_load_balance(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_exit();
 */
    
    void p_stopping_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_exit();
 */
    
    void p_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_exit(CkReductionMsg* impl_msg);
 */
    
    void r_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_control_sync_count(int entry_point, int id, int count);
 */
    
    void p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_enter();
 */
    
    void p_adapt_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_end();
 */
    
    void p_adapt_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_update();
 */
    
    void p_adapt_update(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_next(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_called();
 */
    
    void p_adapt_called(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_exit();
 */
    
    void p_adapt_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_delete();
 */
    
    void p_adapt_delete(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
    
    void p_adapt_recv_level(MsgAdapt* impl_msg) ;

/* DECLS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
    
    void p_adapt_recv_child(MsgCoarsen* impl_msg) ;

/* DECLS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
    
    void r_restart_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
    
    void p_refresh_recv(MsgRefresh* impl_msg) ;

/* DECLS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
    
    void p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: Block(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_Block : public CProxy_ArrayElement{
  public:
    typedef Block local_t;
    typedef CkIndex_Block index_t;
    typedef CProxy_Block proxy_t;
    typedef CProxyElement_Block element_t;
    typedef CProxySection_Block section_t;

    using array_index_t = CkArrayIndexIndex;
    CProxy_Block(void) {
    }
    CProxy_Block(const ArrayElement *e) : CProxy_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_ArrayElement::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_Block operator [] (const CkArrayIndexIndex &idx) const
    { return CProxyElement_Block(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_Block operator() (const CkArrayIndexIndex &idx) const
    { return CProxyElement_Block(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxy_Block(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_Block(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
    
    static CkArrayID ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: Block();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_initial_exit();
 */
    
    void p_initial_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
    
    void r_end_initialize(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_next(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_enter();
 */
    
    void p_output_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_enter(CkReductionMsg* impl_msg);
 */
    
    void r_output_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_end();
 */
    
    void p_output_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_output_exit();
 */
    
    void p_output_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_exit(CkReductionMsg* impl_msg);
 */
    
    void r_output_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_write(int index_output, int step);
 */
    
    void p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_enter();
 */
    
    void p_compute_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_continue();
 */
    
    void p_compute_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
    
    void r_compute_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_compute_exit();
 */
    
    void p_compute_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
    
    void r_compute_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_flux_correct_refresh();
 */
    
    void p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_complete(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
    
    void p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_order_morton_index(int index, int count);
 */
    
    void p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_output_next(MsgOutput* impl_msg);
 */
    
    void p_method_output_next(MsgOutput* impl_msg) ;

/* DECLS: void p_method_output_write(MsgOutput* impl_msg);
 */
    
    void p_method_output_write(MsgOutput* impl_msg) ;

/* DECLS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_done(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_debug_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_compute_timestep(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_enter();
 */
    
    void p_stopping_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_load_balance();
 */
    
    void p_stopping_load_balance(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_load_balance(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_exit();
 */
    
    void p_stopping_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_exit();
 */
    
    void p_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_exit(CkReductionMsg* impl_msg);
 */
    
    void r_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_control_sync_count(int entry_point, int id, int count);
 */
    
    void p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_enter();
 */
    
    void p_adapt_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_end();
 */
    
    void p_adapt_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_update();
 */
    
    void p_adapt_update(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_next(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_called();
 */
    
    void p_adapt_called(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_exit();
 */
    
    void p_adapt_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_delete();
 */
    
    void p_adapt_delete(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
    
    void p_adapt_recv_level(MsgAdapt* impl_msg) ;

/* DECLS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
    
    void p_adapt_recv_child(MsgCoarsen* impl_msg) ;

/* DECLS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
    
    void r_restart_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
    
    void p_refresh_recv(MsgRefresh* impl_msg) ;

/* DECLS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
    
    void p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: Block(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_Block : public CProxySection_ArrayElement{
  public:
    typedef Block local_t;
    typedef CkIndex_Block index_t;
    typedef CProxy_Block proxy_t;
    typedef CProxyElement_Block element_t;
    typedef CProxySection_Block section_t;

    using array_index_t = CkArrayIndexIndex;
    CProxySection_Block(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_ArrayElement::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_ArrayElement::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_ArrayElement::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_ArrayElement::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_ArrayElement::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_ArrayElement::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_ArrayElement::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_ArrayElement::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_ArrayElement::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_ArrayElement::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_ArrayElement::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_Block operator [] (const CkArrayIndexIndex &idx) const
        {return CProxyElement_Block(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_Block operator() (const CkArrayIndexIndex &idx) const
        {return CProxyElement_Block(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_Block(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_Block(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_Block(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems,nElems, factor) {}
    CProxySection_Block(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_Block(const CkSectionID &sid)  
        :CProxySection_ArrayElement(sid) { ckAutoDelegate(); }
    CProxySection_Block(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_Block(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_Block(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_Block(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_ArrayElement(aid,elems) { ckAutoDelegate(); }
    CProxySection_Block(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_ArrayElement(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_Block(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems) {
       return CkSectionID(aid, elems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor) {
      return CkSectionID(aid, elems, nElems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor) {
      return CkSectionID(aid, elems, factor);
    } 
    void ckAutoDelegate(int opts=1) {
      if(ckIsDelegated()) return;
      CProxySection_ArrayElement::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_ArrayElement::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_ArrayElement::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: Block(const process_type &ip_source, const MsgType &msg_type);
 */
    

/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: Block();
 */
    

/* DECLS: void p_initial_exit();
 */
    
    void p_initial_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_end_initialize(CkReductionMsg* impl_msg);
 */
    
    void r_end_initialize(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_next(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_next(CkReductionMsg* impl_msg) ;

/* DECLS: void r_initial_new_continue(CkReductionMsg* impl_msg);
 */
    
    void r_initial_new_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_enter();
 */
    
    void p_output_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_enter(CkReductionMsg* impl_msg);
 */
    
    void r_output_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_end();
 */
    
    void p_output_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_output_exit();
 */
    
    void p_output_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_output_exit(CkReductionMsg* impl_msg);
 */
    
    void r_output_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_output_write(int index_output, int step);
 */
    
    void p_output_write(int index_output, int step, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_enter();
 */
    
    void p_compute_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_compute_continue();
 */
    
    void p_compute_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_continue(CkReductionMsg* impl_msg);
 */
    
    void r_compute_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void p_compute_exit();
 */
    
    void p_compute_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_compute_exit(CkReductionMsg* impl_msg);
 */
    
    void r_compute_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_flux_correct_refresh();
 */
    
    void p_method_flux_correct_refresh(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_flux_correct_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_order_morton_complete(CkReductionMsg* impl_msg);
 */
    
    void r_method_order_morton_complete(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_order_morton_weight(const int *ic3, int weight, const Index &index);
 */
    
    void p_method_order_morton_weight(const int *ic3, int weight, const Index &index, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_order_morton_index(int index, int count);
 */
    
    void p_method_order_morton_index(int index, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_output_next(MsgOutput* impl_msg);
 */
    
    void p_method_output_next(MsgOutput* impl_msg) ;

/* DECLS: void p_method_output_write(MsgOutput* impl_msg);
 */
    
    void p_method_output_write(MsgOutput* impl_msg) ;

/* DECLS: void r_method_output_continue(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_continue(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_output_done(CkReductionMsg* impl_msg);
 */
    
    void r_method_output_done(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_debug_sum_fields(CkReductionMsg* impl_msg);
 */
    
    void r_method_debug_sum_fields(CkReductionMsg* impl_msg) ;

/* DECLS: void r_stopping_compute_timestep(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_compute_timestep(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_enter();
 */
    
    void p_stopping_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_enter(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_load_balance();
 */
    
    void p_stopping_load_balance(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_load_balance(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_load_balance(CkReductionMsg* impl_msg) ;

/* DECLS: void p_stopping_exit();
 */
    
    void p_stopping_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_stopping_exit(CkReductionMsg* impl_msg);
 */
    
    void r_stopping_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_exit();
 */
    
    void p_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_exit(CkReductionMsg* impl_msg);
 */
    
    void r_exit(CkReductionMsg* impl_msg) ;

/* DECLS: void p_control_sync_count(int entry_point, int id, int count);
 */
    
    void p_control_sync_count(int entry_point, int id, int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_enter();
 */
    
    void p_adapt_enter(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_enter(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_end();
 */
    
    void p_adapt_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_update();
 */
    
    void p_adapt_update(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_adapt_next(CkReductionMsg* impl_msg);
 */
    
    void r_adapt_next(CkReductionMsg* impl_msg) ;

/* DECLS: void p_adapt_called();
 */
    
    void p_adapt_called(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_exit();
 */
    
    void p_adapt_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_delete();
 */
    
    void p_adapt_delete(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_adapt_recv_level(MsgAdapt* impl_msg);
 */
    
    void p_adapt_recv_level(MsgAdapt* impl_msg) ;

/* DECLS: void p_adapt_recv_child(MsgCoarsen* impl_msg);
 */
    
    void p_adapt_recv_child(MsgCoarsen* impl_msg) ;

/* DECLS: void r_restart_enter(CkReductionMsg* impl_msg);
 */
    
    void r_restart_enter(CkReductionMsg* impl_msg) ;

/* DECLS: void p_refresh_recv(MsgRefresh* impl_msg);
 */
    
    void p_refresh_recv(MsgRefresh* impl_msg) ;

/* DECLS: void p_refresh_child(int n, const char *a, const int *ic3);
 */
    
    void p_refresh_child(int n, const char *a, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: Block(CkMigrateMessage* impl_msg);
 */

};
#define Block_SDAG_CODE 
typedef CBaseT1<ArrayElementT<Index>, CProxy_Block>CBase_Block;

/* DECLS: array IoReader: ArrayElement{
IoReader();
IoReader(CkMigrateMessage* impl_msg);
};
 */
 class IoReader;
 class CkIndex_IoReader;
 class CProxy_IoReader;
 class CProxyElement_IoReader;
 class CProxySection_IoReader;
/* --------------- index object ------------------ */
class CkIndex_IoReader:public CkIndex_ArrayElement{
  public:
    typedef IoReader local_t;
    typedef CkIndex_IoReader index_t;
    typedef CProxy_IoReader proxy_t;
    typedef CProxyElement_IoReader element_t;
    typedef CProxySection_IoReader section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: IoReader();
     */
    // Entry point registration at startup
    
    static int reg_IoReader_void();
    // Entry point index lookup
    
    inline static int idx_IoReader_void() {
      static int epidx = reg_IoReader_void();
      return epidx;
    }

    
    static int ckNew() { return idx_IoReader_void(); }
    
    static void _call_IoReader_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoReader_void(void* impl_msg, void* impl_obj);
    /* DECLS: IoReader(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_IoReader_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_IoReader_CkMigrateMessage() {
      static int epidx = reg_IoReader_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_IoReader_CkMigrateMessage(); }
    
    static void _call_IoReader_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoReader_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_IoReader : public CProxyElement_ArrayElement{
  public:
    typedef IoReader local_t;
    typedef CkIndex_IoReader index_t;
    typedef CProxy_IoReader proxy_t;
    typedef CProxyElement_IoReader element_t;
    typedef CProxySection_IoReader section_t;

    using array_index_t = CkArrayIndex1D;

    /* TRAM aggregators */

    CProxyElement_IoReader(void) {
    }
    CProxyElement_IoReader(const ArrayElement *e) : CProxyElement_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_ArrayElement::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_ArrayElement::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_ArrayElement::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_ArrayElement::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_ArrayElement::ckGetIndex(); }

    IoReader *ckLocal(void) const
    { return (IoReader *)CProxyElement_ArrayElement::ckLocal(); }

    CProxyElement_IoReader(const CkArrayID &aid,const CkArrayIndex1D &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoReader(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}

    CProxyElement_IoReader(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoReader(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}
/* DECLS: IoReader();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: IoReader(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_IoReader : public CProxy_ArrayElement{
  public:
    typedef IoReader local_t;
    typedef CkIndex_IoReader index_t;
    typedef CProxy_IoReader proxy_t;
    typedef CProxyElement_IoReader element_t;
    typedef CProxySection_IoReader section_t;

    using array_index_t = CkArrayIndex1D;
    CProxy_IoReader(void) {
    }
    CProxy_IoReader(const ArrayElement *e) : CProxy_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_ArrayElement::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_IoReader operator [] (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoReader operator() (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoReader operator [] (int idx) const 
        {return CProxyElement_IoReader(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxyElement_IoReader operator () (int idx) const 
        {return CProxyElement_IoReader(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxy_IoReader(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_IoReader(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: IoReader();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);
    static CkArrayID ckNew(const int s1, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: IoReader(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_IoReader : public CProxySection_ArrayElement{
  public:
    typedef IoReader local_t;
    typedef CkIndex_IoReader index_t;
    typedef CProxy_IoReader proxy_t;
    typedef CProxyElement_IoReader element_t;
    typedef CProxySection_IoReader section_t;

    using array_index_t = CkArrayIndex1D;
    CProxySection_IoReader(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_ArrayElement::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_ArrayElement::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_ArrayElement::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_ArrayElement::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_ArrayElement::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_ArrayElement::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_ArrayElement::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_ArrayElement::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_ArrayElement::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_ArrayElement::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_ArrayElement::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_IoReader operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoReader operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoReader operator [] (int idx) const 
        {return CProxyElement_IoReader(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    CProxyElement_IoReader operator () (int idx) const 
        {return CProxyElement_IoReader(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex1D *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      return CkSectionID(aid, elems, nElems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex1D> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      return CkSectionID(aid, elems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, int l, int u, int s, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      std::vector<CkArrayIndex1D> al;
      for (int i=l; i<=u; i+=s) al.emplace_back(i);
      return CkSectionID(aid, al, factor);
    } 
    CProxySection_IoReader(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoReader(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoReader(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems,nElems, factor) {}
    CProxySection_IoReader(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_IoReader(const CkSectionID &sid)  
        :CProxySection_ArrayElement(sid) { ckAutoDelegate(); }
    CProxySection_IoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_IoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_ArrayElement(aid,elems) { ckAutoDelegate(); }
    CProxySection_IoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_ArrayElement(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_IoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems) {
       return CkSectionID(aid, elems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor) {
      return CkSectionID(aid, elems, nElems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor) {
      return CkSectionID(aid, elems, factor);
    } 
    void ckAutoDelegate(int opts=1) {
      if(ckIsDelegated()) return;
      CProxySection_ArrayElement::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_ArrayElement::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_ArrayElement::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: IoReader();
 */
    

/* DECLS: IoReader(CkMigrateMessage* impl_msg);
 */

};
#define IoReader_SDAG_CODE 
typedef CBaseT1<ArrayElementT<CkIndex1D>, CProxy_IoReader>CBase_IoReader;

/* DECLS: array IoWriter: ArrayElement{
IoWriter();
IoWriter(CkMigrateMessage* impl_msg);
};
 */
 class IoWriter;
 class CkIndex_IoWriter;
 class CProxy_IoWriter;
 class CProxyElement_IoWriter;
 class CProxySection_IoWriter;
/* --------------- index object ------------------ */
class CkIndex_IoWriter:public CkIndex_ArrayElement{
  public:
    typedef IoWriter local_t;
    typedef CkIndex_IoWriter index_t;
    typedef CProxy_IoWriter proxy_t;
    typedef CProxyElement_IoWriter element_t;
    typedef CProxySection_IoWriter section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: IoWriter();
     */
    // Entry point registration at startup
    
    static int reg_IoWriter_void();
    // Entry point index lookup
    
    inline static int idx_IoWriter_void() {
      static int epidx = reg_IoWriter_void();
      return epidx;
    }

    
    static int ckNew() { return idx_IoWriter_void(); }
    
    static void _call_IoWriter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoWriter_void(void* impl_msg, void* impl_obj);
    /* DECLS: IoWriter(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_IoWriter_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_IoWriter_CkMigrateMessage() {
      static int epidx = reg_IoWriter_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_IoWriter_CkMigrateMessage(); }
    
    static void _call_IoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_IoWriter : public CProxyElement_ArrayElement{
  public:
    typedef IoWriter local_t;
    typedef CkIndex_IoWriter index_t;
    typedef CProxy_IoWriter proxy_t;
    typedef CProxyElement_IoWriter element_t;
    typedef CProxySection_IoWriter section_t;

    using array_index_t = CkArrayIndex1D;

    /* TRAM aggregators */

    CProxyElement_IoWriter(void) {
    }
    CProxyElement_IoWriter(const ArrayElement *e) : CProxyElement_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_ArrayElement::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_ArrayElement::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_ArrayElement::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_ArrayElement::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_ArrayElement::ckGetIndex(); }

    IoWriter *ckLocal(void) const
    { return (IoWriter *)CProxyElement_ArrayElement::ckLocal(); }

    CProxyElement_IoWriter(const CkArrayID &aid,const CkArrayIndex1D &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoWriter(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}

    CProxyElement_IoWriter(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoWriter(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}
/* DECLS: IoWriter();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: IoWriter(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_IoWriter : public CProxy_ArrayElement{
  public:
    typedef IoWriter local_t;
    typedef CkIndex_IoWriter index_t;
    typedef CProxy_IoWriter proxy_t;
    typedef CProxyElement_IoWriter element_t;
    typedef CProxySection_IoWriter section_t;

    using array_index_t = CkArrayIndex1D;
    CProxy_IoWriter(void) {
    }
    CProxy_IoWriter(const ArrayElement *e) : CProxy_ArrayElement(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_ArrayElement::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_IoWriter operator [] (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoWriter operator() (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoWriter operator [] (int idx) const 
        {return CProxyElement_IoWriter(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxyElement_IoWriter operator () (int idx) const 
        {return CProxyElement_IoWriter(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxy_IoWriter(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_IoWriter(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: IoWriter();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);
    static CkArrayID ckNew(const int s1, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: IoWriter(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_IoWriter : public CProxySection_ArrayElement{
  public:
    typedef IoWriter local_t;
    typedef CkIndex_IoWriter index_t;
    typedef CProxy_IoWriter proxy_t;
    typedef CProxyElement_IoWriter element_t;
    typedef CProxySection_IoWriter section_t;

    using array_index_t = CkArrayIndex1D;
    CProxySection_IoWriter(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_ArrayElement::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_ArrayElement::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_ArrayElement::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_ArrayElement::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_ArrayElement::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_ArrayElement::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_ArrayElement::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_ArrayElement::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_ArrayElement::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_ArrayElement::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_ArrayElement::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_ArrayElement::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_ArrayElement::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_ArrayElement::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_ArrayElement::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_ArrayElement::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_ArrayElement::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_ArrayElement::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_ArrayElement::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_ArrayElement::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_ArrayElement::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_ArrayElement::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_ArrayElement::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_ArrayElement::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_ArrayElement::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_ArrayElement::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_ArrayElement::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_ArrayElement::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_ArrayElement::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_IoWriter operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoWriter operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoWriter operator [] (int idx) const 
        {return CProxyElement_IoWriter(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    CProxyElement_IoWriter operator () (int idx) const 
        {return CProxyElement_IoWriter(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex1D *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      return CkSectionID(aid, elems, nElems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex1D> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      return CkSectionID(aid, elems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, int l, int u, int s, int factor=USE_DEFAULT_BRANCH_FACTOR) {
      std::vector<CkArrayIndex1D> al;
      for (int i=l; i<=u; i+=s) al.emplace_back(i);
      return CkSectionID(aid, al, factor);
    } 
    CProxySection_IoWriter(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoWriter(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoWriter(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems,nElems, factor) {}
    CProxySection_IoWriter(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_IoWriter(const CkSectionID &sid)  
        :CProxySection_ArrayElement(sid) { ckAutoDelegate(); }
    CProxySection_IoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_IoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_ArrayElement(aid,elems) { ckAutoDelegate(); }
    CProxySection_IoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_ArrayElement(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_IoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems) {
      return CkSectionID(aid, elems, nElems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems) {
       return CkSectionID(aid, elems);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor) {
      return CkSectionID(aid, elems, nElems, factor);
    } 
    static CkSectionID ckNew(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor) {
      return CkSectionID(aid, elems, factor);
    } 
    void ckAutoDelegate(int opts=1) {
      if(ckIsDelegated()) return;
      CProxySection_ArrayElement::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_ArrayElement::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_ArrayElement::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: IoWriter();
 */
    

/* DECLS: IoWriter(CkMigrateMessage* impl_msg);
 */

};
#define IoWriter_SDAG_CODE 
typedef CBaseT1<ArrayElementT<CkIndex1D>, CProxy_IoWriter>CBase_IoWriter;













































































/* ---------------- method closures -------------- */
class Closure_Block {
  public:




    struct p_initial_exit_4_closure;





    struct p_output_enter_8_closure;



    struct p_output_end_10_closure;


    struct p_output_exit_11_closure;



    struct p_output_write_13_closure;


    struct p_compute_enter_14_closure;


    struct p_compute_continue_15_closure;



    struct p_compute_exit_17_closure;



    struct p_method_flux_correct_refresh_19_closure;





    struct p_method_order_morton_weight_23_closure;


    struct p_method_order_morton_index_24_closure;








    struct p_stopping_enter_31_closure;



    struct p_stopping_load_balance_33_closure;



    struct p_stopping_exit_35_closure;



    struct p_exit_37_closure;



    struct p_control_sync_count_39_closure;


    struct p_adapt_enter_40_closure;



    struct p_adapt_end_42_closure;


    struct p_adapt_update_43_closure;



    struct p_adapt_called_45_closure;


    struct p_adapt_exit_46_closure;


    struct p_adapt_delete_47_closure;






    struct p_refresh_child_52_closure;


};

/* ---------------- method closures -------------- */
class Closure_IoReader {
  public:


};

/* ---------------- method closures -------------- */
class Closure_IoWriter {
  public:


};

extern void _registermesh(void);
#endif
