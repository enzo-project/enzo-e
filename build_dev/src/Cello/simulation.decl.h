#ifndef _DECL_simulation_H_
#define _DECL_simulation_H_
#include "charm++.h"
#include "envelope.h"
#include <memory>
#include "sdag.h"
#include "mesh.decl.h"

/* DECLS: readonly CProxy_Simulation proxy_simulation;
 */


/* DECLS: group Simulation: IrrGroup{
Simulation(const char *filename, int n);
void p_get_msg_refine(const Index &impl_noname_0);
void r_initialize_block_array(CkReductionMsg* impl_msg);
void s_write();
void r_write(CkReductionMsg* impl_msg);
void r_write_checkpoint_output();
void r_restart_start(CkReductionMsg* impl_msg);
void p_restart_enter(const std::string &dir);
void p_output_write(int n, const char *buffer);
void r_output_barrier(CkReductionMsg* impl_msg);
void p_output_start(int index_output);
void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
void p_monitor_performance();
void p_set_block_array(const CProxy_Block &block_array);
Simulation(CkMigrateMessage* impl_msg);
};
 */
 class Simulation;
 class CkIndex_Simulation;
 class CProxy_Simulation;
 class CProxyElement_Simulation;
 class CProxySection_Simulation;
/* --------------- index object ------------------ */
class CkIndex_Simulation:public CkIndex_IrrGroup{
  public:
    typedef Simulation local_t;
    typedef CkIndex_Simulation index_t;
    typedef CProxy_Simulation proxy_t;
    typedef CProxyElement_Simulation element_t;
    typedef CProxySection_Simulation section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: Simulation(const char *filename, int n);
     */
    // Entry point registration at startup
    
    static int reg_Simulation_marshall1();
    // Entry point index lookup
    
    inline static int idx_Simulation_marshall1() {
      static int epidx = reg_Simulation_marshall1();
      return epidx;
    }

    
    static int ckNew(const char *filename, int n) { return idx_Simulation_marshall1(); }
    
    static void _call_Simulation_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Simulation_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_Simulation_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_Simulation_marshall1(PUP::er &p,void *msg);
    /* DECLS: void p_get_msg_refine(const Index &impl_noname_0);
     */
    // Entry point registration at startup
    
    static int reg_p_get_msg_refine_marshall2();
    // Entry point index lookup
    
    inline static int idx_p_get_msg_refine_marshall2() {
      static int epidx = reg_p_get_msg_refine_marshall2();
      return epidx;
    }

    
    inline static int idx_p_get_msg_refine(void (Simulation::*)(const Index &impl_noname_0) ) {
      return idx_p_get_msg_refine_marshall2();
    }


    
    static int p_get_msg_refine(const Index &impl_noname_0) { return idx_p_get_msg_refine_marshall2(); }
    
    static void _call_p_get_msg_refine_marshall2(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_get_msg_refine_marshall2(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_get_msg_refine_marshall2(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_get_msg_refine_marshall2(PUP::er &p,void *msg);
    /* DECLS: void r_initialize_block_array(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_initialize_block_array_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_initialize_block_array_CkReductionMsg() {
      static int epidx = reg_r_initialize_block_array_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_initialize_block_array(void (Simulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_initialize_block_array_CkReductionMsg();
    }


    
    static int r_initialize_block_array(CkReductionMsg* impl_msg) { return idx_r_initialize_block_array_CkReductionMsg(); }
    
    static void _call_r_initialize_block_array_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_initialize_block_array_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void s_write();
     */
    // Entry point registration at startup
    
    static int reg_s_write_void();
    // Entry point index lookup
    
    inline static int idx_s_write_void() {
      static int epidx = reg_s_write_void();
      return epidx;
    }

    
    inline static int idx_s_write(void (Simulation::*)() ) {
      return idx_s_write_void();
    }


    
    static int s_write() { return idx_s_write_void(); }
    
    static void _call_s_write_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_s_write_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_write(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_write_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_write_CkReductionMsg() {
      static int epidx = reg_r_write_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_write(void (Simulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_write_CkReductionMsg();
    }


    
    static int r_write(CkReductionMsg* impl_msg) { return idx_r_write_CkReductionMsg(); }
    
    static void _call_r_write_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_write_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_write_checkpoint_output();
     */
    // Entry point registration at startup
    
    static int reg_r_write_checkpoint_output_void();
    // Entry point index lookup
    
    inline static int idx_r_write_checkpoint_output_void() {
      static int epidx = reg_r_write_checkpoint_output_void();
      return epidx;
    }

    
    inline static int idx_r_write_checkpoint_output(void (Simulation::*)() ) {
      return idx_r_write_checkpoint_output_void();
    }


    
    static int r_write_checkpoint_output() { return idx_r_write_checkpoint_output_void(); }
    
    static void _call_r_write_checkpoint_output_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_write_checkpoint_output_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_restart_start(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_restart_start_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_restart_start_CkReductionMsg() {
      static int epidx = reg_r_restart_start_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_restart_start(void (Simulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_restart_start_CkReductionMsg();
    }


    
    static int r_restart_start(CkReductionMsg* impl_msg) { return idx_r_restart_start_CkReductionMsg(); }
    
    static void _call_r_restart_start_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_restart_start_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_restart_enter(const std::string &dir);
     */
    // Entry point registration at startup
    
    static int reg_p_restart_enter_marshall8();
    // Entry point index lookup
    
    inline static int idx_p_restart_enter_marshall8() {
      static int epidx = reg_p_restart_enter_marshall8();
      return epidx;
    }

    
    inline static int idx_p_restart_enter(void (Simulation::*)(const std::string &dir) ) {
      return idx_p_restart_enter_marshall8();
    }


    
    static int p_restart_enter(const std::string &dir) { return idx_p_restart_enter_marshall8(); }
    
    static void _call_p_restart_enter_marshall8(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_enter_marshall8(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_restart_enter_marshall8(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_restart_enter_marshall8(PUP::er &p,void *msg);
    /* DECLS: void p_output_write(int n, const char *buffer);
     */
    // Entry point registration at startup
    
    static int reg_p_output_write_marshall9();
    // Entry point index lookup
    
    inline static int idx_p_output_write_marshall9() {
      static int epidx = reg_p_output_write_marshall9();
      return epidx;
    }

    
    inline static int idx_p_output_write(void (Simulation::*)(int n, const char *buffer) ) {
      return idx_p_output_write_marshall9();
    }


    
    static int p_output_write(int n, const char *buffer) { return idx_p_output_write_marshall9(); }
    
    static void _call_p_output_write_marshall9(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_write_marshall9(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_output_write_marshall9(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_output_write_marshall9(PUP::er &p,void *msg);
    /* DECLS: void r_output_barrier(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_output_barrier_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_output_barrier_CkReductionMsg() {
      static int epidx = reg_r_output_barrier_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_output_barrier(void (Simulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_output_barrier_CkReductionMsg();
    }


    
    static int r_output_barrier(CkReductionMsg* impl_msg) { return idx_r_output_barrier_CkReductionMsg(); }
    
    static void _call_r_output_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_output_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_start(int index_output);
     */
    // Entry point registration at startup
    
    static int reg_p_output_start_marshall11();
    // Entry point index lookup
    
    inline static int idx_p_output_start_marshall11() {
      static int epidx = reg_p_output_start_marshall11();
      return epidx;
    }

    
    inline static int idx_p_output_start(void (Simulation::*)(int index_output) ) {
      return idx_p_output_start_marshall11();
    }


    
    static int p_output_start(int index_output) { return idx_p_output_start_marshall11(); }
    
    static void _call_p_output_start_marshall11(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_start_marshall11(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_output_start_marshall11(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_output_start_marshall11(PUP::er &p,void *msg);
    /* DECLS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_monitor_performance_reduce_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_monitor_performance_reduce_CkReductionMsg() {
      static int epidx = reg_r_monitor_performance_reduce_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_monitor_performance_reduce(void (Simulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_monitor_performance_reduce_CkReductionMsg();
    }


    
    static int r_monitor_performance_reduce(CkReductionMsg* impl_msg) { return idx_r_monitor_performance_reduce_CkReductionMsg(); }
    
    static void _call_r_monitor_performance_reduce_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_monitor_performance_reduce_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_monitor_performance();
     */
    // Entry point registration at startup
    
    static int reg_p_monitor_performance_void();
    // Entry point index lookup
    
    inline static int idx_p_monitor_performance_void() {
      static int epidx = reg_p_monitor_performance_void();
      return epidx;
    }

    
    inline static int idx_p_monitor_performance(void (Simulation::*)() ) {
      return idx_p_monitor_performance_void();
    }


    
    static int p_monitor_performance() { return idx_p_monitor_performance_void(); }
    
    static void _call_p_monitor_performance_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_monitor_performance_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_set_block_array(const CProxy_Block &block_array);
     */
    // Entry point registration at startup
    
    static int reg_p_set_block_array_marshall14();
    // Entry point index lookup
    
    inline static int idx_p_set_block_array_marshall14() {
      static int epidx = reg_p_set_block_array_marshall14();
      return epidx;
    }

    
    inline static int idx_p_set_block_array(void (Simulation::*)(const CProxy_Block &block_array) ) {
      return idx_p_set_block_array_marshall14();
    }


    
    static int p_set_block_array(const CProxy_Block &block_array) { return idx_p_set_block_array_marshall14(); }
    
    static void _call_p_set_block_array_marshall14(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_block_array_marshall14(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_set_block_array_marshall14(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_set_block_array_marshall14(PUP::er &p,void *msg);
    /* DECLS: Simulation(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_Simulation_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_Simulation_CkMigrateMessage() {
      static int epidx = reg_Simulation_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_Simulation_CkMigrateMessage(); }
    
    static void _call_Simulation_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Simulation_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_Simulation: public CProxyElement_IrrGroup{
  public:
    typedef Simulation local_t;
    typedef CkIndex_Simulation index_t;
    typedef CProxy_Simulation proxy_t;
    typedef CProxyElement_Simulation element_t;
    typedef CProxySection_Simulation section_t;


    /* TRAM aggregators */

    CProxyElement_Simulation(void) {
    }
    CProxyElement_Simulation(const IrrGroup *g) : CProxyElement_IrrGroup(g){
    }
    CProxyElement_Simulation(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_IrrGroup(_gid,_onPE,CK_DELCTOR_ARGS){
    }
    CProxyElement_Simulation(CkGroupID _gid,int _onPE) : CProxyElement_IrrGroup(_gid,_onPE){
    }

    int ckIsDelegated(void) const
    { return CProxyElement_IrrGroup::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_IrrGroup::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_IrrGroup::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_IrrGroup::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxyElement_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IrrGroup::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IrrGroup::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_IrrGroup::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_IrrGroup::ckGetGroupPe();}

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_IrrGroup::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_IrrGroup::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_IrrGroup::ckSetGroupID(g);
    }
    Simulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static Simulation* ckLocalBranch(CkGroupID gID) {
      return (Simulation*)CkLocalBranch(gID);
    }
/* DECLS: Simulation(const char *filename, int n);
 */
    

/* DECLS: void p_get_msg_refine(const Index &impl_noname_0);
 */
    
    void p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
    
    void r_initialize_block_array(CkReductionMsg* impl_msg);

/* DECLS: void s_write();
 */
    
    void s_write(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_write(CkReductionMsg* impl_msg);
 */
    
    void r_write(CkReductionMsg* impl_msg);

/* DECLS: void r_write_checkpoint_output();
 */
    
    void r_write_checkpoint_output(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_restart_start(CkReductionMsg* impl_msg);
 */
    
    void r_restart_start(CkReductionMsg* impl_msg);

/* DECLS: void p_restart_enter(const std::string &dir);
 */
    
    void p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_output_write(int n, const char *buffer);
 */
    
    void p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_output_barrier(CkReductionMsg* impl_msg);

/* DECLS: void p_output_start(int index_output);
 */
    
    void p_output_start(int index_output, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
    
    void r_monitor_performance_reduce(CkReductionMsg* impl_msg);

/* DECLS: void p_monitor_performance();
 */
    
    void p_monitor_performance(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_block_array(const CProxy_Block &block_array);
 */
    
    void p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: Simulation(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
class CProxy_Simulation: public CProxy_IrrGroup{
  public:
    typedef Simulation local_t;
    typedef CkIndex_Simulation index_t;
    typedef CProxy_Simulation proxy_t;
    typedef CProxyElement_Simulation element_t;
    typedef CProxySection_Simulation section_t;

    CProxy_Simulation(void) {
    }
    CProxy_Simulation(const IrrGroup *g) : CProxy_IrrGroup(g){
    }
    CProxy_Simulation(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_IrrGroup(_gid,CK_DELCTOR_ARGS){  }
    CProxy_Simulation(CkGroupID _gid) : CProxy_IrrGroup(_gid){  }
    CProxyElement_Simulation operator[](int onPE) const
      {return CProxyElement_Simulation(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}

    int ckIsDelegated(void) const
    { return CProxy_IrrGroup::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_IrrGroup::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_IrrGroup::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_IrrGroup::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxy_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IrrGroup::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IrrGroup::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_IrrGroup::ckSetReductionClient(cb); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_IrrGroup::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_IrrGroup::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_IrrGroup::ckSetGroupID(g);
    }
    Simulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static Simulation* ckLocalBranch(CkGroupID gID) {
      return (Simulation*)CkLocalBranch(gID);
    }
/* DECLS: Simulation(const char *filename, int n);
 */
    
    static CkGroupID ckNew(const char *filename, int n, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_Simulation(const char *filename, int n, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_get_msg_refine(const Index &impl_noname_0);
 */
    
    void p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_refine(const Index &impl_noname_0, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_refine(const Index &impl_noname_0, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
    
    void r_initialize_block_array(CkReductionMsg* impl_msg);
    
    void r_initialize_block_array(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_initialize_block_array(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void s_write();
 */
    
    void s_write(const CkEntryOptions *impl_e_opts=NULL);
    
    void s_write(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void s_write(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_write(CkReductionMsg* impl_msg);
 */
    
    void r_write(CkReductionMsg* impl_msg);
    
    void r_write(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_write(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void r_write_checkpoint_output();
 */
    
    void r_write_checkpoint_output(const CkEntryOptions *impl_e_opts=NULL);
    
    void r_write_checkpoint_output(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void r_write_checkpoint_output(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_restart_start(CkReductionMsg* impl_msg);
 */
    
    void r_restart_start(CkReductionMsg* impl_msg);
    
    void r_restart_start(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_restart_start(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_restart_enter(const std::string &dir);
 */
    
    void p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_enter(const std::string &dir, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_enter(const std::string &dir, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_output_write(int n, const char *buffer);
 */
    
    void p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_output_write(int n, const char *buffer, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_output_write(int n, const char *buffer, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_output_barrier(CkReductionMsg* impl_msg);
    
    void r_output_barrier(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_output_barrier(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_output_start(int index_output);
 */
    
    void p_output_start(int index_output, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_output_start(int index_output, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_output_start(int index_output, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
    
    void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
    
    void r_monitor_performance_reduce(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_monitor_performance_reduce(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_monitor_performance();
 */
    
    void p_monitor_performance(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_monitor_performance(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_monitor_performance(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_block_array(const CProxy_Block &block_array);
 */
    
    void p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_block_array(const CProxy_Block &block_array, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_block_array(const CProxy_Block &block_array, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: Simulation(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
class CProxySection_Simulation: public CProxySection_IrrGroup{
  public:
    typedef Simulation local_t;
    typedef CkIndex_Simulation index_t;
    typedef CProxy_Simulation proxy_t;
    typedef CProxyElement_Simulation element_t;
    typedef CProxySection_Simulation section_t;

    CProxySection_Simulation(void) {
    }
    CProxySection_Simulation(const IrrGroup *g) : CProxySection_IrrGroup(g){
    }
    CProxySection_Simulation(const CkGroupID &_gid,const int *_pelist,int _npes, CK_DELCTOR_PARAM) : CProxySection_IrrGroup(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_Simulation(const CkGroupID &_gid,const int *_pelist,int _npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_IrrGroup(_gid,_pelist,_npes,factor){  }
    CProxySection_Simulation(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_IrrGroup(n,_gid,_pelist,_npes,factor){  }
    CProxySection_Simulation(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, CK_DELCTOR_PARAM) : CProxySection_IrrGroup(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }

    int ckIsDelegated(void) const
    { return CProxySection_IrrGroup::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_IrrGroup::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_IrrGroup::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_IrrGroup::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxySection_IrrGroup::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_IrrGroup::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_IrrGroup::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IrrGroup::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IrrGroup::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_IrrGroup::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_IrrGroup::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_IrrGroup::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_IrrGroup::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_IrrGroup::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_IrrGroup::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_IrrGroup::ckGetGroupIDn(i); }
inline const int *ckGetElements() const
{ return CProxySection_IrrGroup::ckGetElements(); }
inline const int *ckGetElements(int i) const
{ return CProxySection_IrrGroup::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_IrrGroup::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_IrrGroup::ckGetNumElements(i); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_IrrGroup::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_IrrGroup::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_IrrGroup::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_IrrGroup::ckSetGroupID(g);
    }
    Simulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static Simulation* ckLocalBranch(CkGroupID gID) {
      return (Simulation*)CkLocalBranch(gID);
    }
/* DECLS: Simulation(const char *filename, int n);
 */
    

/* DECLS: void p_get_msg_refine(const Index &impl_noname_0);
 */
    
    void p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
    
    void r_initialize_block_array(CkReductionMsg* impl_msg);

/* DECLS: void s_write();
 */
    
    void s_write(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_write(CkReductionMsg* impl_msg);
 */
    
    void r_write(CkReductionMsg* impl_msg);

/* DECLS: void r_write_checkpoint_output();
 */
    
    void r_write_checkpoint_output(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_restart_start(CkReductionMsg* impl_msg);
 */
    
    void r_restart_start(CkReductionMsg* impl_msg);

/* DECLS: void p_restart_enter(const std::string &dir);
 */
    
    void p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_output_write(int n, const char *buffer);
 */
    
    void p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_output_barrier(CkReductionMsg* impl_msg);

/* DECLS: void p_output_start(int index_output);
 */
    
    void p_output_start(int index_output, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
    
    void r_monitor_performance_reduce(CkReductionMsg* impl_msg);

/* DECLS: void p_monitor_performance();
 */
    
    void p_monitor_performance(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_block_array(const CProxy_Block &block_array);
 */
    
    void p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: Simulation(CkMigrateMessage* impl_msg);
 */

};
#define Simulation_SDAG_CODE 
typedef CBaseT1<Group, CProxy_Simulation>CBase_Simulation;

/* DECLS: group MappingArray: CkArrayMap{
MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
MappingArray(CkMigrateMessage* impl_msg);
};
 */
 class MappingArray;
 class CkIndex_MappingArray;
 class CProxy_MappingArray;
 class CProxyElement_MappingArray;
 class CProxySection_MappingArray;
/* --------------- index object ------------------ */
class CkIndex_MappingArray:public CkIndex_CkArrayMap{
  public:
    typedef MappingArray local_t;
    typedef CkIndex_MappingArray index_t;
    typedef CProxy_MappingArray proxy_t;
    typedef CProxyElement_MappingArray element_t;
    typedef CProxySection_MappingArray section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
     */
    // Entry point registration at startup
    
    static int reg_MappingArray_marshall1();
    // Entry point index lookup
    
    inline static int idx_MappingArray_marshall1() {
      static int epidx = reg_MappingArray_marshall1();
      return epidx;
    }

    
    static int ckNew(int impl_noname_1, int impl_noname_2, int impl_noname_3) { return idx_MappingArray_marshall1(); }
    
    static void _call_MappingArray_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingArray_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_MappingArray_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_MappingArray_marshall1(PUP::er &p,void *msg);
    /* DECLS: MappingArray(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_MappingArray_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_MappingArray_CkMigrateMessage() {
      static int epidx = reg_MappingArray_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_MappingArray_CkMigrateMessage(); }
    
    static void _call_MappingArray_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingArray_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_MappingArray: public CProxyElement_CkArrayMap{
  public:
    typedef MappingArray local_t;
    typedef CkIndex_MappingArray index_t;
    typedef CProxy_MappingArray proxy_t;
    typedef CProxyElement_MappingArray element_t;
    typedef CProxySection_MappingArray section_t;


    /* TRAM aggregators */

    CProxyElement_MappingArray(void) {
    }
    CProxyElement_MappingArray(const IrrGroup *g) : CProxyElement_CkArrayMap(g){
    }
    CProxyElement_MappingArray(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_CkArrayMap(_gid,_onPE,CK_DELCTOR_ARGS){
    }
    CProxyElement_MappingArray(CkGroupID _gid,int _onPE) : CProxyElement_CkArrayMap(_gid,_onPE){
    }

    int ckIsDelegated(void) const
    { return CProxyElement_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxyElement_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_CkArrayMap::ckGetGroupPe();}

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_CkArrayMap::ckSetGroupID(g);
    }
    MappingArray* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingArray* ckLocalBranch(CkGroupID gID) {
      return (MappingArray*)CkLocalBranch(gID);
    }
/* DECLS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
    

/* DECLS: MappingArray(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
class CProxy_MappingArray: public CProxy_CkArrayMap{
  public:
    typedef MappingArray local_t;
    typedef CkIndex_MappingArray index_t;
    typedef CProxy_MappingArray proxy_t;
    typedef CProxyElement_MappingArray element_t;
    typedef CProxySection_MappingArray section_t;

    CProxy_MappingArray(void) {
    }
    CProxy_MappingArray(const IrrGroup *g) : CProxy_CkArrayMap(g){
    }
    CProxy_MappingArray(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_CkArrayMap(_gid,CK_DELCTOR_ARGS){  }
    CProxy_MappingArray(CkGroupID _gid) : CProxy_CkArrayMap(_gid){  }
    CProxyElement_MappingArray operator[](int onPE) const
      {return CProxyElement_MappingArray(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}

    int ckIsDelegated(void) const
    { return CProxy_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxy_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_CkArrayMap::ckSetReductionClient(cb); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_CkArrayMap::ckSetGroupID(g);
    }
    MappingArray* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingArray* ckLocalBranch(CkGroupID gID) {
      return (MappingArray*)CkLocalBranch(gID);
    }
/* DECLS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
    
    static CkGroupID ckNew(int impl_noname_1, int impl_noname_2, int impl_noname_3, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: MappingArray(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
class CProxySection_MappingArray: public CProxySection_CkArrayMap{
  public:
    typedef MappingArray local_t;
    typedef CkIndex_MappingArray index_t;
    typedef CProxy_MappingArray proxy_t;
    typedef CProxyElement_MappingArray element_t;
    typedef CProxySection_MappingArray section_t;

    CProxySection_MappingArray(void) {
    }
    CProxySection_MappingArray(const IrrGroup *g) : CProxySection_CkArrayMap(g){
    }
    CProxySection_MappingArray(const CkGroupID &_gid,const int *_pelist,int _npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_MappingArray(const CkGroupID &_gid,const int *_pelist,int _npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(_gid,_pelist,_npes,factor){  }
    CProxySection_MappingArray(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,factor){  }
    CProxySection_MappingArray(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }

    int ckIsDelegated(void) const
    { return CProxySection_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxySection_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_CkArrayMap::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_CkArrayMap::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_CkArrayMap::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_CkArrayMap::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_CkArrayMap::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_CkArrayMap::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_CkArrayMap::ckGetGroupIDn(i); }
inline const int *ckGetElements() const
{ return CProxySection_CkArrayMap::ckGetElements(); }
inline const int *ckGetElements(int i) const
{ return CProxySection_CkArrayMap::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_CkArrayMap::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_CkArrayMap::ckGetNumElements(i); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_CkArrayMap::ckSetGroupID(g);
    }
    MappingArray* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingArray* ckLocalBranch(CkGroupID gID) {
      return (MappingArray*)CkLocalBranch(gID);
    }
/* DECLS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
    

/* DECLS: MappingArray(CkMigrateMessage* impl_msg);
 */

};
#define MappingArray_SDAG_CODE 
typedef CBaseT1<CkArrayMap, CProxy_MappingArray>CBase_MappingArray;

/* DECLS: group MappingTree: CkArrayMap{
MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
MappingTree(CkMigrateMessage* impl_msg);
};
 */
 class MappingTree;
 class CkIndex_MappingTree;
 class CProxy_MappingTree;
 class CProxyElement_MappingTree;
 class CProxySection_MappingTree;
/* --------------- index object ------------------ */
class CkIndex_MappingTree:public CkIndex_CkArrayMap{
  public:
    typedef MappingTree local_t;
    typedef CkIndex_MappingTree index_t;
    typedef CProxy_MappingTree proxy_t;
    typedef CProxyElement_MappingTree element_t;
    typedef CProxySection_MappingTree section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
     */
    // Entry point registration at startup
    
    static int reg_MappingTree_marshall1();
    // Entry point index lookup
    
    inline static int idx_MappingTree_marshall1() {
      static int epidx = reg_MappingTree_marshall1();
      return epidx;
    }

    
    static int ckNew(int impl_noname_4, int impl_noname_5, int impl_noname_6) { return idx_MappingTree_marshall1(); }
    
    static void _call_MappingTree_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingTree_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_MappingTree_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_MappingTree_marshall1(PUP::er &p,void *msg);
    /* DECLS: MappingTree(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_MappingTree_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_MappingTree_CkMigrateMessage() {
      static int epidx = reg_MappingTree_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_MappingTree_CkMigrateMessage(); }
    
    static void _call_MappingTree_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingTree_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_MappingTree: public CProxyElement_CkArrayMap{
  public:
    typedef MappingTree local_t;
    typedef CkIndex_MappingTree index_t;
    typedef CProxy_MappingTree proxy_t;
    typedef CProxyElement_MappingTree element_t;
    typedef CProxySection_MappingTree section_t;


    /* TRAM aggregators */

    CProxyElement_MappingTree(void) {
    }
    CProxyElement_MappingTree(const IrrGroup *g) : CProxyElement_CkArrayMap(g){
    }
    CProxyElement_MappingTree(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_CkArrayMap(_gid,_onPE,CK_DELCTOR_ARGS){
    }
    CProxyElement_MappingTree(CkGroupID _gid,int _onPE) : CProxyElement_CkArrayMap(_gid,_onPE){
    }

    int ckIsDelegated(void) const
    { return CProxyElement_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxyElement_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_CkArrayMap::ckGetGroupPe();}

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_CkArrayMap::ckSetGroupID(g);
    }
    MappingTree* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingTree* ckLocalBranch(CkGroupID gID) {
      return (MappingTree*)CkLocalBranch(gID);
    }
/* DECLS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
    

/* DECLS: MappingTree(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
class CProxy_MappingTree: public CProxy_CkArrayMap{
  public:
    typedef MappingTree local_t;
    typedef CkIndex_MappingTree index_t;
    typedef CProxy_MappingTree proxy_t;
    typedef CProxyElement_MappingTree element_t;
    typedef CProxySection_MappingTree section_t;

    CProxy_MappingTree(void) {
    }
    CProxy_MappingTree(const IrrGroup *g) : CProxy_CkArrayMap(g){
    }
    CProxy_MappingTree(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_CkArrayMap(_gid,CK_DELCTOR_ARGS){  }
    CProxy_MappingTree(CkGroupID _gid) : CProxy_CkArrayMap(_gid){  }
    CProxyElement_MappingTree operator[](int onPE) const
      {return CProxyElement_MappingTree(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}

    int ckIsDelegated(void) const
    { return CProxy_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxy_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_CkArrayMap::ckSetReductionClient(cb); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_CkArrayMap::ckSetGroupID(g);
    }
    MappingTree* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingTree* ckLocalBranch(CkGroupID gID) {
      return (MappingTree*)CkLocalBranch(gID);
    }
/* DECLS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
    
    static CkGroupID ckNew(int impl_noname_4, int impl_noname_5, int impl_noname_6, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: MappingTree(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
class CProxySection_MappingTree: public CProxySection_CkArrayMap{
  public:
    typedef MappingTree local_t;
    typedef CkIndex_MappingTree index_t;
    typedef CProxy_MappingTree proxy_t;
    typedef CProxyElement_MappingTree element_t;
    typedef CProxySection_MappingTree section_t;

    CProxySection_MappingTree(void) {
    }
    CProxySection_MappingTree(const IrrGroup *g) : CProxySection_CkArrayMap(g){
    }
    CProxySection_MappingTree(const CkGroupID &_gid,const int *_pelist,int _npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_MappingTree(const CkGroupID &_gid,const int *_pelist,int _npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(_gid,_pelist,_npes,factor){  }
    CProxySection_MappingTree(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,factor){  }
    CProxySection_MappingTree(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }

    int ckIsDelegated(void) const
    { return CProxySection_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxySection_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_CkArrayMap::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_CkArrayMap::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_CkArrayMap::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_CkArrayMap::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_CkArrayMap::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_CkArrayMap::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_CkArrayMap::ckGetGroupIDn(i); }
inline const int *ckGetElements() const
{ return CProxySection_CkArrayMap::ckGetElements(); }
inline const int *ckGetElements(int i) const
{ return CProxySection_CkArrayMap::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_CkArrayMap::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_CkArrayMap::ckGetNumElements(i); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_CkArrayMap::ckSetGroupID(g);
    }
    MappingTree* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingTree* ckLocalBranch(CkGroupID gID) {
      return (MappingTree*)CkLocalBranch(gID);
    }
/* DECLS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
    

/* DECLS: MappingTree(CkMigrateMessage* impl_msg);
 */

};
#define MappingTree_SDAG_CODE 
typedef CBaseT1<CkArrayMap, CProxy_MappingTree>CBase_MappingTree;

/* DECLS: group MappingIo: CkArrayMap{
MappingIo(int impl_noname_7);
MappingIo(CkMigrateMessage* impl_msg);
};
 */
 class MappingIo;
 class CkIndex_MappingIo;
 class CProxy_MappingIo;
 class CProxyElement_MappingIo;
 class CProxySection_MappingIo;
/* --------------- index object ------------------ */
class CkIndex_MappingIo:public CkIndex_CkArrayMap{
  public:
    typedef MappingIo local_t;
    typedef CkIndex_MappingIo index_t;
    typedef CProxy_MappingIo proxy_t;
    typedef CProxyElement_MappingIo element_t;
    typedef CProxySection_MappingIo section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: MappingIo(int impl_noname_7);
     */
    // Entry point registration at startup
    
    static int reg_MappingIo_marshall1();
    // Entry point index lookup
    
    inline static int idx_MappingIo_marshall1() {
      static int epidx = reg_MappingIo_marshall1();
      return epidx;
    }

    
    static int ckNew(int impl_noname_7) { return idx_MappingIo_marshall1(); }
    
    static void _call_MappingIo_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingIo_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_MappingIo_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_MappingIo_marshall1(PUP::er &p,void *msg);
    /* DECLS: MappingIo(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_MappingIo_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_MappingIo_CkMigrateMessage() {
      static int epidx = reg_MappingIo_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_MappingIo_CkMigrateMessage(); }
    
    static void _call_MappingIo_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_MappingIo_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_MappingIo: public CProxyElement_CkArrayMap{
  public:
    typedef MappingIo local_t;
    typedef CkIndex_MappingIo index_t;
    typedef CProxy_MappingIo proxy_t;
    typedef CProxyElement_MappingIo element_t;
    typedef CProxySection_MappingIo section_t;


    /* TRAM aggregators */

    CProxyElement_MappingIo(void) {
    }
    CProxyElement_MappingIo(const IrrGroup *g) : CProxyElement_CkArrayMap(g){
    }
    CProxyElement_MappingIo(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_CkArrayMap(_gid,_onPE,CK_DELCTOR_ARGS){
    }
    CProxyElement_MappingIo(CkGroupID _gid,int _onPE) : CProxyElement_CkArrayMap(_gid,_onPE){
    }

    int ckIsDelegated(void) const
    { return CProxyElement_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxyElement_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_CkArrayMap::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_CkArrayMap::ckGetGroupPe();}

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_CkArrayMap::ckSetGroupID(g);
    }
    MappingIo* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingIo* ckLocalBranch(CkGroupID gID) {
      return (MappingIo*)CkLocalBranch(gID);
    }
/* DECLS: MappingIo(int impl_noname_7);
 */
    

/* DECLS: MappingIo(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
class CProxy_MappingIo: public CProxy_CkArrayMap{
  public:
    typedef MappingIo local_t;
    typedef CkIndex_MappingIo index_t;
    typedef CProxy_MappingIo proxy_t;
    typedef CProxyElement_MappingIo element_t;
    typedef CProxySection_MappingIo section_t;

    CProxy_MappingIo(void) {
    }
    CProxy_MappingIo(const IrrGroup *g) : CProxy_CkArrayMap(g){
    }
    CProxy_MappingIo(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_CkArrayMap(_gid,CK_DELCTOR_ARGS){  }
    CProxy_MappingIo(CkGroupID _gid) : CProxy_CkArrayMap(_gid){  }
    CProxyElement_MappingIo operator[](int onPE) const
      {return CProxyElement_MappingIo(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}

    int ckIsDelegated(void) const
    { return CProxy_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxy_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_CkArrayMap::ckSetReductionClient(cb); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_CkArrayMap::ckSetGroupID(g);
    }
    MappingIo* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingIo* ckLocalBranch(CkGroupID gID) {
      return (MappingIo*)CkLocalBranch(gID);
    }
/* DECLS: MappingIo(int impl_noname_7);
 */
    
    static CkGroupID ckNew(int impl_noname_7, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_MappingIo(int impl_noname_7, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: MappingIo(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
class CProxySection_MappingIo: public CProxySection_CkArrayMap{
  public:
    typedef MappingIo local_t;
    typedef CkIndex_MappingIo index_t;
    typedef CProxy_MappingIo proxy_t;
    typedef CProxyElement_MappingIo element_t;
    typedef CProxySection_MappingIo section_t;

    CProxySection_MappingIo(void) {
    }
    CProxySection_MappingIo(const IrrGroup *g) : CProxySection_CkArrayMap(g){
    }
    CProxySection_MappingIo(const CkGroupID &_gid,const int *_pelist,int _npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_MappingIo(const CkGroupID &_gid,const int *_pelist,int _npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(_gid,_pelist,_npes,factor){  }
    CProxySection_MappingIo(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,factor){  }
    CProxySection_MappingIo(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, CK_DELCTOR_PARAM) : CProxySection_CkArrayMap(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }

    int ckIsDelegated(void) const
    { return CProxySection_CkArrayMap::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_CkArrayMap::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_CkArrayMap::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_CkArrayMap::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxySection_CkArrayMap::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_CkArrayMap::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_CkArrayMap::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_CkArrayMap::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_CkArrayMap::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_CkArrayMap::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_CkArrayMap::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_CkArrayMap::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_CkArrayMap::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_CkArrayMap::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_CkArrayMap::ckGetGroupIDn(i); }
inline const int *ckGetElements() const
{ return CProxySection_CkArrayMap::ckGetElements(); }
inline const int *ckGetElements(int i) const
{ return CProxySection_CkArrayMap::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_CkArrayMap::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_CkArrayMap::ckGetNumElements(i); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_CkArrayMap::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_CkArrayMap::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_CkArrayMap::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_CkArrayMap::ckSetGroupID(g);
    }
    MappingIo* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static MappingIo* ckLocalBranch(CkGroupID gID) {
      return (MappingIo*)CkLocalBranch(gID);
    }
/* DECLS: MappingIo(int impl_noname_7);
 */
    

/* DECLS: MappingIo(CkMigrateMessage* impl_msg);
 */

};
#define MappingIo_SDAG_CODE 
typedef CBaseT1<CkArrayMap, CProxy_MappingIo>CBase_MappingIo;




/* ---------------- method closures -------------- */
class Closure_Simulation {
  public:


    struct p_get_msg_refine_2_closure;



    struct s_write_4_closure;



    struct r_write_checkpoint_output_6_closure;



    struct p_restart_enter_8_closure;


    struct p_output_write_9_closure;



    struct p_output_start_11_closure;



    struct p_monitor_performance_13_closure;


    struct p_set_block_array_14_closure;


};

/* ---------------- method closures -------------- */
class Closure_MappingArray {
  public:


};

/* ---------------- method closures -------------- */
class Closure_MappingTree {
  public:


};

/* ---------------- method closures -------------- */
class Closure_MappingIo {
  public:


};

extern void _registersimulation(void);
#endif
