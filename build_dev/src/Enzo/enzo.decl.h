#ifndef _DECL_enzo_H_
#define _DECL_enzo_H_
#include "charm++.h"
#include "envelope.h"
#include <memory>
#include "sdag.h"



/* DECLS: readonly int EnzoMsgCheck::counter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly EnzoConfig g_enzo_config;
 */

/* DECLS: readonly int EnzoBlock::UseMinimumPressureSupport[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly enzo_float EnzoBlock::MinimumPressureSupportParameter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::MultiSpecies[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::PressureFree[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly enzo_float EnzoBlock::GravitationalConstant[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::ProblemType[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::PPMFlatteningParameter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::PPMDiffusionParameter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::PPMSteepeningParameter[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly enzo_float EnzoBlock::InitialRedshift[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly enzo_float EnzoBlock::InitialTimeInCodeUnits[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly enzo_float EnzoBlock::DomainLeftEdge[CONFIG_NODE_SIZE_3];
 */

/* DECLS: readonly enzo_float EnzoBlock::DomainRightEdge[CONFIG_NODE_SIZE_3];
 */

/* DECLS: readonly int EnzoBlock::GridRank[CONFIG_NODE_SIZE];
 */

/* DECLS: readonly int EnzoBlock::ghost_depth[CONFIG_NODE_SIZE_3];
 */

/* DECLS: readonly int EnzoBlock::NumberOfBaryonFields[CONFIG_NODE_SIZE];
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
/* DECLS: message EnzoMsgCheck;
 */
class EnzoMsgCheck;
class CMessage_EnzoMsgCheck:public CkMessage{
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
    CMessage_EnzoMsgCheck();
    static void *pack(EnzoMsgCheck *p);
    static EnzoMsgCheck* unpack(void* p);
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

#include "mesh.decl.h"

/* DECLS: readonly CProxy_EnzoSimulation proxy_enzo_simulation;
 */

/* DECLS: readonly CProxy_IoEnzoReader proxy_io_enzo_reader;
 */

/* DECLS: readonly CProxy_IoEnzoWriter proxy_io_enzo_writer;
 */

/* DECLS: readonly CProxy_EnzoLevelArray proxy_level_array;
 */

/* DECLS: group EnzoSimulation: Simulation{
EnzoSimulation(const char *filename, int n);
void r_startup_begun(CkReductionMsg* impl_msg);
void p_get_msg_refine(const Index &index);
void p_get_msg_check(const Index &index);
void r_method_balance_count(CkReductionMsg* impl_msg);
void p_method_balance_check();
void r_method_check_enter(CkReductionMsg* impl_msg);
void p_check_done();
void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
void p_infer_set_array_count(int count);
void p_infer_array_created();
void p_infer_done();
void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
void p_io_reader_created();
void p_restart_next_level();
void p_restart_level_created();
void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
void p_fbnet_concatenate_sphere_lists();
void p_fbnet_done();
EnzoSimulation(CkMigrateMessage* impl_msg);
};
 */
 class EnzoSimulation;
 class CkIndex_EnzoSimulation;
 class CProxy_EnzoSimulation;
 class CProxyElement_EnzoSimulation;
 class CProxySection_EnzoSimulation;
/* --------------- index object ------------------ */
class CkIndex_EnzoSimulation:public CkIndex_Simulation{
  public:
    typedef EnzoSimulation local_t;
    typedef CkIndex_EnzoSimulation index_t;
    typedef CProxy_EnzoSimulation proxy_t;
    typedef CProxyElement_EnzoSimulation element_t;
    typedef CProxySection_EnzoSimulation section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: EnzoSimulation(const char *filename, int n);
     */
    // Entry point registration at startup
    
    static int reg_EnzoSimulation_marshall1();
    // Entry point index lookup
    
    inline static int idx_EnzoSimulation_marshall1() {
      static int epidx = reg_EnzoSimulation_marshall1();
      return epidx;
    }

    
    static int ckNew(const char *filename, int n) { return idx_EnzoSimulation_marshall1(); }
    
    static void _call_EnzoSimulation_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoSimulation_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_EnzoSimulation_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_EnzoSimulation_marshall1(PUP::er &p,void *msg);
    /* DECLS: void r_startup_begun(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_startup_begun_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_startup_begun_CkReductionMsg() {
      static int epidx = reg_r_startup_begun_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_startup_begun(void (EnzoSimulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_startup_begun_CkReductionMsg();
    }


    
    static int r_startup_begun(CkReductionMsg* impl_msg) { return idx_r_startup_begun_CkReductionMsg(); }
    
    static void _call_r_startup_begun_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_startup_begun_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_get_msg_refine(const Index &index);
     */
    // Entry point registration at startup
    
    static int reg_p_get_msg_refine_marshall3();
    // Entry point index lookup
    
    inline static int idx_p_get_msg_refine_marshall3() {
      static int epidx = reg_p_get_msg_refine_marshall3();
      return epidx;
    }

    
    inline static int idx_p_get_msg_refine(void (EnzoSimulation::*)(const Index &index) ) {
      return idx_p_get_msg_refine_marshall3();
    }


    
    static int p_get_msg_refine(const Index &index) { return idx_p_get_msg_refine_marshall3(); }
    
    static void _call_p_get_msg_refine_marshall3(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_get_msg_refine_marshall3(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_get_msg_refine_marshall3(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_get_msg_refine_marshall3(PUP::er &p,void *msg);
    /* DECLS: void p_get_msg_check(const Index &index);
     */
    // Entry point registration at startup
    
    static int reg_p_get_msg_check_marshall4();
    // Entry point index lookup
    
    inline static int idx_p_get_msg_check_marshall4() {
      static int epidx = reg_p_get_msg_check_marshall4();
      return epidx;
    }

    
    inline static int idx_p_get_msg_check(void (EnzoSimulation::*)(const Index &index) ) {
      return idx_p_get_msg_check_marshall4();
    }


    
    static int p_get_msg_check(const Index &index) { return idx_p_get_msg_check_marshall4(); }
    
    static void _call_p_get_msg_check_marshall4(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_get_msg_check_marshall4(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_get_msg_check_marshall4(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_get_msg_check_marshall4(PUP::er &p,void *msg);
    /* DECLS: void r_method_balance_count(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_balance_count_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_balance_count_CkReductionMsg() {
      static int epidx = reg_r_method_balance_count_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_balance_count(void (EnzoSimulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_balance_count_CkReductionMsg();
    }


    
    static int r_method_balance_count(CkReductionMsg* impl_msg) { return idx_r_method_balance_count_CkReductionMsg(); }
    
    static void _call_r_method_balance_count_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_balance_count_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_balance_check();
     */
    // Entry point registration at startup
    
    static int reg_p_method_balance_check_void();
    // Entry point index lookup
    
    inline static int idx_p_method_balance_check_void() {
      static int epidx = reg_p_method_balance_check_void();
      return epidx;
    }

    
    inline static int idx_p_method_balance_check(void (EnzoSimulation::*)() ) {
      return idx_p_method_balance_check_void();
    }


    
    static int p_method_balance_check() { return idx_p_method_balance_check_void(); }
    
    static void _call_p_method_balance_check_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_balance_check_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_check_enter(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_check_enter_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_check_enter_CkReductionMsg() {
      static int epidx = reg_r_method_check_enter_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_check_enter(void (EnzoSimulation::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_check_enter_CkReductionMsg();
    }


    
    static int r_method_check_enter(CkReductionMsg* impl_msg) { return idx_r_method_check_enter_CkReductionMsg(); }
    
    static void _call_r_method_check_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_check_enter_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_check_done();
     */
    // Entry point registration at startup
    
    static int reg_p_check_done_void();
    // Entry point index lookup
    
    inline static int idx_p_check_done_void() {
      static int epidx = reg_p_check_done_void();
      return epidx;
    }

    
    inline static int idx_p_check_done(void (EnzoSimulation::*)() ) {
      return idx_p_check_done_void();
    }


    
    static int p_check_done() { return idx_p_check_done_void(); }
    
    static void _call_p_check_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_check_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
     */
    // Entry point registration at startup
    
    static int reg_p_set_io_writer_marshall9();
    // Entry point index lookup
    
    inline static int idx_p_set_io_writer_marshall9() {
      static int epidx = reg_p_set_io_writer_marshall9();
      return epidx;
    }

    
    inline static int idx_p_set_io_writer(void (EnzoSimulation::*)(const CProxy_IoEnzoWriter &proxy) ) {
      return idx_p_set_io_writer_marshall9();
    }


    
    static int p_set_io_writer(const CProxy_IoEnzoWriter &proxy) { return idx_p_set_io_writer_marshall9(); }
    
    static void _call_p_set_io_writer_marshall9(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_io_writer_marshall9(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_set_io_writer_marshall9(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_set_io_writer_marshall9(PUP::er &p,void *msg);
    /* DECLS: void p_infer_set_array_count(int count);
     */
    // Entry point registration at startup
    
    static int reg_p_infer_set_array_count_marshall10();
    // Entry point index lookup
    
    inline static int idx_p_infer_set_array_count_marshall10() {
      static int epidx = reg_p_infer_set_array_count_marshall10();
      return epidx;
    }

    
    inline static int idx_p_infer_set_array_count(void (EnzoSimulation::*)(int count) ) {
      return idx_p_infer_set_array_count_marshall10();
    }


    
    static int p_infer_set_array_count(int count) { return idx_p_infer_set_array_count_marshall10(); }
    
    static void _call_p_infer_set_array_count_marshall10(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_infer_set_array_count_marshall10(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_infer_set_array_count_marshall10(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_infer_set_array_count_marshall10(PUP::er &p,void *msg);
    /* DECLS: void p_infer_array_created();
     */
    // Entry point registration at startup
    
    static int reg_p_infer_array_created_void();
    // Entry point index lookup
    
    inline static int idx_p_infer_array_created_void() {
      static int epidx = reg_p_infer_array_created_void();
      return epidx;
    }

    
    inline static int idx_p_infer_array_created(void (EnzoSimulation::*)() ) {
      return idx_p_infer_array_created_void();
    }


    
    static int p_infer_array_created() { return idx_p_infer_array_created_void(); }
    
    static void _call_p_infer_array_created_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_infer_array_created_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_infer_done();
     */
    // Entry point registration at startup
    
    static int reg_p_infer_done_void();
    // Entry point index lookup
    
    inline static int idx_p_infer_done_void() {
      static int epidx = reg_p_infer_done_void();
      return epidx;
    }

    
    inline static int idx_p_infer_done(void (EnzoSimulation::*)() ) {
      return idx_p_infer_done_void();
    }


    
    static int p_infer_done() { return idx_p_infer_done_void(); }
    
    static void _call_p_infer_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_infer_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
     */
    // Entry point registration at startup
    
    static int reg_p_set_io_reader_marshall13();
    // Entry point index lookup
    
    inline static int idx_p_set_io_reader_marshall13() {
      static int epidx = reg_p_set_io_reader_marshall13();
      return epidx;
    }

    
    inline static int idx_p_set_io_reader(void (EnzoSimulation::*)(const CProxy_IoEnzoReader &proxy) ) {
      return idx_p_set_io_reader_marshall13();
    }


    
    static int p_set_io_reader(const CProxy_IoEnzoReader &proxy) { return idx_p_set_io_reader_marshall13(); }
    
    static void _call_p_set_io_reader_marshall13(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_io_reader_marshall13(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_set_io_reader_marshall13(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_set_io_reader_marshall13(PUP::er &p,void *msg);
    /* DECLS: void p_io_reader_created();
     */
    // Entry point registration at startup
    
    static int reg_p_io_reader_created_void();
    // Entry point index lookup
    
    inline static int idx_p_io_reader_created_void() {
      static int epidx = reg_p_io_reader_created_void();
      return epidx;
    }

    
    inline static int idx_p_io_reader_created(void (EnzoSimulation::*)() ) {
      return idx_p_io_reader_created_void();
    }


    
    static int p_io_reader_created() { return idx_p_io_reader_created_void(); }
    
    static void _call_p_io_reader_created_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_io_reader_created_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_restart_next_level();
     */
    // Entry point registration at startup
    
    static int reg_p_restart_next_level_void();
    // Entry point index lookup
    
    inline static int idx_p_restart_next_level_void() {
      static int epidx = reg_p_restart_next_level_void();
      return epidx;
    }

    
    inline static int idx_p_restart_next_level(void (EnzoSimulation::*)() ) {
      return idx_p_restart_next_level_void();
    }


    
    static int p_restart_next_level() { return idx_p_restart_next_level_void(); }
    
    static void _call_p_restart_next_level_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_next_level_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_restart_level_created();
     */
    // Entry point registration at startup
    
    static int reg_p_restart_level_created_void();
    // Entry point index lookup
    
    inline static int idx_p_restart_level_created_void() {
      static int epidx = reg_p_restart_level_created_void();
      return epidx;
    }

    
    inline static int idx_p_restart_level_created(void (EnzoSimulation::*)() ) {
      return idx_p_restart_level_created_void();
    }


    
    static int p_restart_level_created() { return idx_p_restart_level_created_void(); }
    
    static void _call_p_restart_level_created_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_level_created_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
     */
    // Entry point registration at startup
    
    static int reg_p_set_level_array_marshall17();
    // Entry point index lookup
    
    inline static int idx_p_set_level_array_marshall17() {
      static int epidx = reg_p_set_level_array_marshall17();
      return epidx;
    }

    
    inline static int idx_p_set_level_array(void (EnzoSimulation::*)(const CProxy_EnzoLevelArray &proxy) ) {
      return idx_p_set_level_array_marshall17();
    }


    
    static int p_set_level_array(const CProxy_EnzoLevelArray &proxy) { return idx_p_set_level_array_marshall17(); }
    
    static void _call_p_set_level_array_marshall17(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_level_array_marshall17(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_set_level_array_marshall17(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_set_level_array_marshall17(PUP::er &p,void *msg);
    /* DECLS: void p_fbnet_concatenate_sphere_lists();
     */
    // Entry point registration at startup
    
    static int reg_p_fbnet_concatenate_sphere_lists_void();
    // Entry point index lookup
    
    inline static int idx_p_fbnet_concatenate_sphere_lists_void() {
      static int epidx = reg_p_fbnet_concatenate_sphere_lists_void();
      return epidx;
    }

    
    inline static int idx_p_fbnet_concatenate_sphere_lists(void (EnzoSimulation::*)() ) {
      return idx_p_fbnet_concatenate_sphere_lists_void();
    }


    
    static int p_fbnet_concatenate_sphere_lists() { return idx_p_fbnet_concatenate_sphere_lists_void(); }
    
    static void _call_p_fbnet_concatenate_sphere_lists_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_fbnet_concatenate_sphere_lists_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_fbnet_done();
     */
    // Entry point registration at startup
    
    static int reg_p_fbnet_done_void();
    // Entry point index lookup
    
    inline static int idx_p_fbnet_done_void() {
      static int epidx = reg_p_fbnet_done_void();
      return epidx;
    }

    
    inline static int idx_p_fbnet_done(void (EnzoSimulation::*)() ) {
      return idx_p_fbnet_done_void();
    }


    
    static int p_fbnet_done() { return idx_p_fbnet_done_void(); }
    
    static void _call_p_fbnet_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_fbnet_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: EnzoSimulation(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_EnzoSimulation_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_EnzoSimulation_CkMigrateMessage() {
      static int epidx = reg_EnzoSimulation_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_EnzoSimulation_CkMigrateMessage(); }
    
    static void _call_EnzoSimulation_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoSimulation_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxyElement_EnzoSimulation: public CProxyElement_Simulation{
  public:
    typedef EnzoSimulation local_t;
    typedef CkIndex_EnzoSimulation index_t;
    typedef CProxy_EnzoSimulation proxy_t;
    typedef CProxyElement_EnzoSimulation element_t;
    typedef CProxySection_EnzoSimulation section_t;


    /* TRAM aggregators */

    CProxyElement_EnzoSimulation(void) {
    }
    CProxyElement_EnzoSimulation(const IrrGroup *g) : CProxyElement_Simulation(g){
    }
    CProxyElement_EnzoSimulation(CkGroupID _gid,int _onPE,CK_DELCTOR_PARAM) : CProxyElement_Simulation(_gid,_onPE,CK_DELCTOR_ARGS){
    }
    CProxyElement_EnzoSimulation(CkGroupID _gid,int _onPE) : CProxyElement_Simulation(_gid,_onPE){
    }

    int ckIsDelegated(void) const
    { return CProxyElement_Simulation::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_Simulation::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_Simulation::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_Simulation::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxyElement_Simulation::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxyElement_Simulation::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxyElement_Simulation::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_Simulation::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_Simulation::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_Simulation::ckSetReductionClient(cb); }
int ckGetGroupPe(void) const
{return CProxyElement_Simulation::ckGetGroupPe();}

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_Simulation::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_Simulation::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_Simulation::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxyElement_Simulation::ckSetGroupID(g);
    }
    EnzoSimulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static EnzoSimulation* ckLocalBranch(CkGroupID gID) {
      return (EnzoSimulation*)CkLocalBranch(gID);
    }
/* DECLS: EnzoSimulation(const char *filename, int n);
 */
    

/* DECLS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
    
    void r_startup_begun(CkReductionMsg* impl_msg);

/* DECLS: void p_get_msg_refine(const Index &index);
 */
    
    void p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_get_msg_check(const Index &index);
 */
    
    void p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
    
    void r_method_balance_count(CkReductionMsg* impl_msg);

/* DECLS: void p_method_balance_check();
 */
    
    void p_method_balance_check(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
    
    void r_method_check_enter(CkReductionMsg* impl_msg);

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
    
    void p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_set_array_count(int count);
 */
    
    void p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_array_created();
 */
    
    void p_infer_array_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_done();
 */
    
    void p_infer_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
    
    void p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_io_reader_created();
 */
    
    void p_io_reader_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_next_level();
 */
    
    void p_restart_next_level(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_level_created();
 */
    
    void p_restart_level_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
    
    void p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_concatenate_sphere_lists();
 */
    
    void p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_done();
 */
    
    void p_fbnet_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
class CProxy_EnzoSimulation: public CProxy_Simulation{
  public:
    typedef EnzoSimulation local_t;
    typedef CkIndex_EnzoSimulation index_t;
    typedef CProxy_EnzoSimulation proxy_t;
    typedef CProxyElement_EnzoSimulation element_t;
    typedef CProxySection_EnzoSimulation section_t;

    CProxy_EnzoSimulation(void) {
    }
    CProxy_EnzoSimulation(const IrrGroup *g) : CProxy_Simulation(g){
    }
    CProxy_EnzoSimulation(CkGroupID _gid,CK_DELCTOR_PARAM) : CProxy_Simulation(_gid,CK_DELCTOR_ARGS){  }
    CProxy_EnzoSimulation(CkGroupID _gid) : CProxy_Simulation(_gid){  }
    CProxyElement_EnzoSimulation operator[](int onPE) const
      {return CProxyElement_EnzoSimulation(ckGetGroupID(),onPE,CK_DELCTOR_CALL);}

    int ckIsDelegated(void) const
    { return CProxy_Simulation::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_Simulation::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_Simulation::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_Simulation::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxy_Simulation::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxy_Simulation::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxy_Simulation::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_Simulation::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_Simulation::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_Simulation::ckSetReductionClient(cb); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_Simulation::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_Simulation::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_Simulation::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxy_Simulation::ckSetGroupID(g);
    }
    EnzoSimulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static EnzoSimulation* ckLocalBranch(CkGroupID gID) {
      return (EnzoSimulation*)CkLocalBranch(gID);
    }
/* DECLS: EnzoSimulation(const char *filename, int n);
 */
    
    static CkGroupID ckNew(const char *filename, int n, const CkEntryOptions *impl_e_opts=NULL);
    CProxy_EnzoSimulation(const char *filename, int n, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
    
    void r_startup_begun(CkReductionMsg* impl_msg);
    
    void r_startup_begun(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_startup_begun(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_get_msg_refine(const Index &index);
 */
    
    void p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_refine(const Index &index, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_refine(const Index &index, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_get_msg_check(const Index &index);
 */
    
    void p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_check(const Index &index, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_get_msg_check(const Index &index, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
    
    void r_method_balance_count(CkReductionMsg* impl_msg);
    
    void r_method_balance_count(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_method_balance_count(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_method_balance_check();
 */
    
    void p_method_balance_check(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_method_balance_check(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_method_balance_check(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
    
    void r_method_check_enter(CkReductionMsg* impl_msg);
    
    void r_method_check_enter(CkReductionMsg* impl_msg, int npes, int *pes);
    
    void r_method_check_enter(CkReductionMsg* impl_msg, CmiGroup &grp);

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_check_done(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_check_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
    
    void p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_io_writer(const CProxy_IoEnzoWriter &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_io_writer(const CProxy_IoEnzoWriter &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_set_array_count(int count);
 */
    
    void p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_set_array_count(int count, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_set_array_count(int count, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_array_created();
 */
    
    void p_infer_array_created(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_array_created(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_array_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_done();
 */
    
    void p_infer_done(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_done(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_infer_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
    
    void p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_io_reader(const CProxy_IoEnzoReader &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_io_reader(const CProxy_IoEnzoReader &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_io_reader_created();
 */
    
    void p_io_reader_created(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_io_reader_created(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_io_reader_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_next_level();
 */
    
    void p_restart_next_level(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_next_level(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_next_level(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_level_created();
 */
    
    void p_restart_level_created(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_level_created(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_restart_level_created(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
    
    void p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_level_array(const CProxy_EnzoLevelArray &proxy, int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_set_level_array(const CProxy_EnzoLevelArray &proxy, CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_concatenate_sphere_lists();
 */
    
    void p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_fbnet_concatenate_sphere_lists(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_fbnet_concatenate_sphere_lists(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_done();
 */
    
    void p_fbnet_done(const CkEntryOptions *impl_e_opts=NULL);
    
    void p_fbnet_done(int npes, int *pes, const CkEntryOptions *impl_e_opts=NULL);
    
    void p_fbnet_done(CmiGroup &grp, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
class CProxySection_EnzoSimulation: public CProxySection_Simulation{
  public:
    typedef EnzoSimulation local_t;
    typedef CkIndex_EnzoSimulation index_t;
    typedef CProxy_EnzoSimulation proxy_t;
    typedef CProxyElement_EnzoSimulation element_t;
    typedef CProxySection_EnzoSimulation section_t;

    CProxySection_EnzoSimulation(void) {
    }
    CProxySection_EnzoSimulation(const IrrGroup *g) : CProxySection_Simulation(g){
    }
    CProxySection_EnzoSimulation(const CkGroupID &_gid,const int *_pelist,int _npes, CK_DELCTOR_PARAM) : CProxySection_Simulation(_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }
    CProxySection_EnzoSimulation(const CkGroupID &_gid,const int *_pelist,int _npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_Simulation(_gid,_pelist,_npes,factor){  }
    CProxySection_EnzoSimulation(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, int factor = USE_DEFAULT_BRANCH_FACTOR) : CProxySection_Simulation(n,_gid,_pelist,_npes,factor){  }
    CProxySection_EnzoSimulation(int n,const CkGroupID *_gid, int const * const *_pelist,const int *_npes, CK_DELCTOR_PARAM) : CProxySection_Simulation(n,_gid,_pelist,_npes,CK_DELCTOR_ARGS){  }

    int ckIsDelegated(void) const
    { return CProxySection_Simulation::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_Simulation::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_Simulation::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_Simulation::ckDelegatedIdx(); }
inline void ckCheck(void) const {CProxySection_Simulation::ckCheck();}
CkChareID ckGetChareID(void) const
   {return CProxySection_Simulation::ckGetChareID();}
CkGroupID ckGetGroupID(void) const
   {return CProxySection_Simulation::ckGetGroupID();}
operator CkGroupID () const { return ckGetGroupID(); }

    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_Simulation::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_Simulation::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_Simulation::ckSetReductionClient(cb); }
inline int ckGetNumSections() const
{ return CProxySection_Simulation::ckGetNumSections(); }
inline CkSectionInfo &ckGetSectionInfo()
{ return CProxySection_Simulation::ckGetSectionInfo(); }
inline CkSectionID *ckGetSectionIDs()
{ return CProxySection_Simulation::ckGetSectionIDs(); }
inline CkSectionID &ckGetSectionID()
{ return CProxySection_Simulation::ckGetSectionID(); }
inline CkSectionID &ckGetSectionID(int i)
{ return CProxySection_Simulation::ckGetSectionID(i); }
inline CkGroupID ckGetGroupIDn(int i) const
{ return CProxySection_Simulation::ckGetGroupIDn(i); }
inline const int *ckGetElements() const
{ return CProxySection_Simulation::ckGetElements(); }
inline const int *ckGetElements(int i) const
{ return CProxySection_Simulation::ckGetElements(i); }
inline int ckGetNumElements() const
{ return CProxySection_Simulation::ckGetNumElements(); } 
inline int ckGetNumElements(int i) const
{ return CProxySection_Simulation::ckGetNumElements(i); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_Simulation::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_Simulation::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_Simulation::pup(p);
    }
    void ckSetGroupID(CkGroupID g) {
      CProxySection_Simulation::ckSetGroupID(g);
    }
    EnzoSimulation* ckLocalBranch(void) const {
      return ckLocalBranch(ckGetGroupID());
    }
    static EnzoSimulation* ckLocalBranch(CkGroupID gID) {
      return (EnzoSimulation*)CkLocalBranch(gID);
    }
/* DECLS: EnzoSimulation(const char *filename, int n);
 */
    

/* DECLS: void r_startup_begun(CkReductionMsg* impl_msg);
 */
    
    void r_startup_begun(CkReductionMsg* impl_msg);

/* DECLS: void p_get_msg_refine(const Index &index);
 */
    
    void p_get_msg_refine(const Index &index, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_get_msg_check(const Index &index);
 */
    
    void p_get_msg_check(const Index &index, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_balance_count(CkReductionMsg* impl_msg);
 */
    
    void r_method_balance_count(CkReductionMsg* impl_msg);

/* DECLS: void p_method_balance_check();
 */
    
    void p_method_balance_check(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void r_method_check_enter(CkReductionMsg* impl_msg);
 */
    
    void r_method_check_enter(CkReductionMsg* impl_msg);

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_writer(const CProxy_IoEnzoWriter &proxy);
 */
    
    void p_set_io_writer(const CProxy_IoEnzoWriter &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_set_array_count(int count);
 */
    
    void p_infer_set_array_count(int count, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_array_created();
 */
    
    void p_infer_array_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_infer_done();
 */
    
    void p_infer_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_io_reader(const CProxy_IoEnzoReader &proxy);
 */
    
    void p_set_io_reader(const CProxy_IoEnzoReader &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_io_reader_created();
 */
    
    void p_io_reader_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_next_level();
 */
    
    void p_restart_next_level(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_restart_level_created();
 */
    
    void p_restart_level_created(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_level_array(const CProxy_EnzoLevelArray &proxy);
 */
    
    void p_set_level_array(const CProxy_EnzoLevelArray &proxy, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_concatenate_sphere_lists();
 */
    
    void p_fbnet_concatenate_sphere_lists(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_fbnet_done();
 */
    
    void p_fbnet_done(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: EnzoSimulation(CkMigrateMessage* impl_msg);
 */

};
#define EnzoSimulation_SDAG_CODE 
typedef CBaseT1<Simulation, CProxy_EnzoSimulation>CBase_EnzoSimulation;

/* DECLS: array EnzoBlock: Block{
EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
void p_set_msg_refine(MsgRefine* impl_msg);
void p_set_msg_check(EnzoMsgCheck* impl_msg);
EnzoBlock();
void p_method_feedback_starss_end();
void p_method_m1_closure_solve_transport_eqn();
void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
void r_method_turbulence_end(CkReductionMsg* impl_msg);
void p_initial_hdf5_recv(MsgInitial* impl_msg);
void p_method_balance_migrate();
void p_method_balance_done();
void p_method_gravity_continue();
void p_method_gravity_end();
void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
void p_method_infer_count_arrays(int count);
void p_method_infer_request_data(const int *il3);
void p_method_infer_update(int n, const char *buffer, const int *il3);
void p_method_infer_exit();
void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
void p_method_fbnet_exit();
void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
void p_check_write_next(int num_files, const std::string &ordering);
void p_check_done();
void p_restart_refine(const int *ic3, int io_reader);
void p_restart_set_data(EnzoMsgCheck* impl_msg);
void p_restart_done();
void p_method_accretion_end();
void p_solver_cg_matvec();
void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
void p_solver_cg_loop_2();
void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
void p_solver_bicgstab_loop_2();
void p_solver_bicgstab_loop_3();
void p_solver_bicgstab_loop_8();
void p_solver_bicgstab_loop_9();
void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
void p_solver_dd_solve_coarse();
void p_solver_dd_solve_domain();
void p_solver_dd_last_smooth();
void r_solver_dd_barrier(CkReductionMsg* impl_msg);
void r_solver_dd_end(CkReductionMsg* impl_msg);
void p_solver_jacobi_continue();
void p_solver_mg0_restrict();
void p_solver_mg0_solve_coarse();
void p_solver_mg0_post_smooth();
void p_solver_mg0_last_smooth();
void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
EnzoBlock(CkMigrateMessage* impl_msg);
};
 */
 class EnzoBlock;
 class CkIndex_EnzoBlock;
 class CProxy_EnzoBlock;
 class CProxyElement_EnzoBlock;
 class CProxySection_EnzoBlock;
/* --------------- index object ------------------ */
class CkIndex_EnzoBlock:public CkIndex_Block{
  public:
    typedef EnzoBlock local_t;
    typedef CkIndex_EnzoBlock index_t;
    typedef CProxy_EnzoBlock proxy_t;
    typedef CProxyElement_EnzoBlock element_t;
    typedef CProxySection_EnzoBlock section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
     */
    // Entry point registration at startup
    
    static int reg_EnzoBlock_marshall1();
    // Entry point index lookup
    
    inline static int idx_EnzoBlock_marshall1() {
      static int epidx = reg_EnzoBlock_marshall1();
      return epidx;
    }

    
    static int ckNew(const process_type &ip_source, const MsgType &msg_type) { return idx_EnzoBlock_marshall1(); }
    
    static void _call_EnzoBlock_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoBlock_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_EnzoBlock_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_EnzoBlock_marshall1(PUP::er &p,void *msg);
    /* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_set_msg_refine_MsgRefine();
    // Entry point index lookup
    
    inline static int idx_p_set_msg_refine_MsgRefine() {
      static int epidx = reg_p_set_msg_refine_MsgRefine();
      return epidx;
    }

    
    inline static int idx_p_set_msg_refine(void (EnzoBlock::*)(MsgRefine* impl_msg) ) {
      return idx_p_set_msg_refine_MsgRefine();
    }


    
    static int p_set_msg_refine(MsgRefine* impl_msg) { return idx_p_set_msg_refine_MsgRefine(); }
    
    static void _call_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_msg_refine_MsgRefine(void* impl_msg, void* impl_obj);
    /* DECLS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_set_msg_check_EnzoMsgCheck();
    // Entry point index lookup
    
    inline static int idx_p_set_msg_check_EnzoMsgCheck() {
      static int epidx = reg_p_set_msg_check_EnzoMsgCheck();
      return epidx;
    }

    
    inline static int idx_p_set_msg_check(void (EnzoBlock::*)(EnzoMsgCheck* impl_msg) ) {
      return idx_p_set_msg_check_EnzoMsgCheck();
    }


    
    static int p_set_msg_check(EnzoMsgCheck* impl_msg) { return idx_p_set_msg_check_EnzoMsgCheck(); }
    
    static void _call_p_set_msg_check_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_set_msg_check_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    /* DECLS: EnzoBlock();
     */
    // Entry point registration at startup
    
    static int reg_EnzoBlock_void();
    // Entry point index lookup
    
    inline static int idx_EnzoBlock_void() {
      static int epidx = reg_EnzoBlock_void();
      return epidx;
    }

    
    static int ckNew() { return idx_EnzoBlock_void(); }
    
    static void _call_EnzoBlock_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoBlock_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_feedback_starss_end();
     */
    // Entry point registration at startup
    
    static int reg_p_method_feedback_starss_end_void();
    // Entry point index lookup
    
    inline static int idx_p_method_feedback_starss_end_void() {
      static int epidx = reg_p_method_feedback_starss_end_void();
      return epidx;
    }

    
    inline static int idx_p_method_feedback_starss_end(void (EnzoBlock::*)() ) {
      return idx_p_method_feedback_starss_end_void();
    }


    
    static int p_method_feedback_starss_end() { return idx_p_method_feedback_starss_end_void(); }
    
    static void _call_p_method_feedback_starss_end_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_feedback_starss_end_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_m1_closure_solve_transport_eqn();
     */
    // Entry point registration at startup
    
    static int reg_p_method_m1_closure_solve_transport_eqn_void();
    // Entry point index lookup
    
    inline static int idx_p_method_m1_closure_solve_transport_eqn_void() {
      static int epidx = reg_p_method_m1_closure_solve_transport_eqn_void();
      return epidx;
    }

    
    inline static int idx_p_method_m1_closure_solve_transport_eqn(void (EnzoBlock::*)() ) {
      return idx_p_method_m1_closure_solve_transport_eqn_void();
    }


    
    static int p_method_m1_closure_solve_transport_eqn() { return idx_p_method_m1_closure_solve_transport_eqn_void(); }
    
    static void _call_p_method_m1_closure_solve_transport_eqn_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_m1_closure_solve_transport_eqn_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_method_m1_closure_set_global_averages_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_p_method_m1_closure_set_global_averages_CkReductionMsg() {
      static int epidx = reg_p_method_m1_closure_set_global_averages_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_p_method_m1_closure_set_global_averages(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_p_method_m1_closure_set_global_averages_CkReductionMsg();
    }


    
    static int p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) { return idx_p_method_m1_closure_set_global_averages_CkReductionMsg(); }
    
    static void _call_p_method_m1_closure_set_global_averages_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_m1_closure_set_global_averages_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_method_turbulence_end_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_method_turbulence_end_CkReductionMsg() {
      static int epidx = reg_r_method_turbulence_end_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_method_turbulence_end(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_method_turbulence_end_CkReductionMsg();
    }


    
    static int r_method_turbulence_end(CkReductionMsg* impl_msg) { return idx_r_method_turbulence_end_CkReductionMsg(); }
    
    static void _call_r_method_turbulence_end_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_method_turbulence_end_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_initial_hdf5_recv_MsgInitial();
    // Entry point index lookup
    
    inline static int idx_p_initial_hdf5_recv_MsgInitial() {
      static int epidx = reg_p_initial_hdf5_recv_MsgInitial();
      return epidx;
    }

    
    inline static int idx_p_initial_hdf5_recv(void (EnzoBlock::*)(MsgInitial* impl_msg) ) {
      return idx_p_initial_hdf5_recv_MsgInitial();
    }


    
    static int p_initial_hdf5_recv(MsgInitial* impl_msg) { return idx_p_initial_hdf5_recv_MsgInitial(); }
    
    static void _call_p_initial_hdf5_recv_MsgInitial(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_initial_hdf5_recv_MsgInitial(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_balance_migrate();
     */
    // Entry point registration at startup
    
    static int reg_p_method_balance_migrate_void();
    // Entry point index lookup
    
    inline static int idx_p_method_balance_migrate_void() {
      static int epidx = reg_p_method_balance_migrate_void();
      return epidx;
    }

    
    inline static int idx_p_method_balance_migrate(void (EnzoBlock::*)() ) {
      return idx_p_method_balance_migrate_void();
    }


    
    static int p_method_balance_migrate() { return idx_p_method_balance_migrate_void(); }
    
    static void _call_p_method_balance_migrate_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_balance_migrate_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_balance_done();
     */
    // Entry point registration at startup
    
    static int reg_p_method_balance_done_void();
    // Entry point index lookup
    
    inline static int idx_p_method_balance_done_void() {
      static int epidx = reg_p_method_balance_done_void();
      return epidx;
    }

    
    inline static int idx_p_method_balance_done(void (EnzoBlock::*)() ) {
      return idx_p_method_balance_done_void();
    }


    
    static int p_method_balance_done() { return idx_p_method_balance_done_void(); }
    
    static void _call_p_method_balance_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_balance_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_gravity_continue();
     */
    // Entry point registration at startup
    
    static int reg_p_method_gravity_continue_void();
    // Entry point index lookup
    
    inline static int idx_p_method_gravity_continue_void() {
      static int epidx = reg_p_method_gravity_continue_void();
      return epidx;
    }

    
    inline static int idx_p_method_gravity_continue(void (EnzoBlock::*)() ) {
      return idx_p_method_gravity_continue_void();
    }


    
    static int p_method_gravity_continue() { return idx_p_method_gravity_continue_void(); }
    
    static void _call_p_method_gravity_continue_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_gravity_continue_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_gravity_end();
     */
    // Entry point registration at startup
    
    static int reg_p_method_gravity_end_void();
    // Entry point index lookup
    
    inline static int idx_p_method_gravity_end_void() {
      static int epidx = reg_p_method_gravity_end_void();
      return epidx;
    }

    
    inline static int idx_p_method_gravity_end(void (EnzoBlock::*)() ) {
      return idx_p_method_gravity_end_void();
    }


    
    static int p_method_gravity_end() { return idx_p_method_gravity_end_void(); }
    
    static void _call_p_method_gravity_end_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_gravity_end_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
     */
    // Entry point registration at startup
    
    static int reg_p_method_infer_merge_masks_marshall14();
    // Entry point index lookup
    
    inline static int idx_p_method_infer_merge_masks_marshall14() {
      static int epidx = reg_p_method_infer_merge_masks_marshall14();
      return epidx;
    }

    
    inline static int idx_p_method_infer_merge_masks(void (EnzoBlock::*)(int n, const char *mask, const int *ic3) ) {
      return idx_p_method_infer_merge_masks_marshall14();
    }


    
    static int p_method_infer_merge_masks(int n, const char *mask, const int *ic3) { return idx_p_method_infer_merge_masks_marshall14(); }
    
    static void _call_p_method_infer_merge_masks_marshall14(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_infer_merge_masks_marshall14(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_infer_merge_masks_marshall14(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_infer_merge_masks_marshall14(PUP::er &p,void *msg);
    /* DECLS: void p_method_infer_count_arrays(int count);
     */
    // Entry point registration at startup
    
    static int reg_p_method_infer_count_arrays_marshall15();
    // Entry point index lookup
    
    inline static int idx_p_method_infer_count_arrays_marshall15() {
      static int epidx = reg_p_method_infer_count_arrays_marshall15();
      return epidx;
    }

    
    inline static int idx_p_method_infer_count_arrays(void (EnzoBlock::*)(int count) ) {
      return idx_p_method_infer_count_arrays_marshall15();
    }


    
    static int p_method_infer_count_arrays(int count) { return idx_p_method_infer_count_arrays_marshall15(); }
    
    static void _call_p_method_infer_count_arrays_marshall15(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_infer_count_arrays_marshall15(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_infer_count_arrays_marshall15(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_infer_count_arrays_marshall15(PUP::er &p,void *msg);
    /* DECLS: void p_method_infer_request_data(const int *il3);
     */
    // Entry point registration at startup
    
    static int reg_p_method_infer_request_data_marshall16();
    // Entry point index lookup
    
    inline static int idx_p_method_infer_request_data_marshall16() {
      static int epidx = reg_p_method_infer_request_data_marshall16();
      return epidx;
    }

    
    inline static int idx_p_method_infer_request_data(void (EnzoBlock::*)(const int *il3) ) {
      return idx_p_method_infer_request_data_marshall16();
    }


    
    static int p_method_infer_request_data(const int *il3) { return idx_p_method_infer_request_data_marshall16(); }
    
    static void _call_p_method_infer_request_data_marshall16(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_infer_request_data_marshall16(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_infer_request_data_marshall16(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_infer_request_data_marshall16(PUP::er &p,void *msg);
    /* DECLS: void p_method_infer_update(int n, const char *buffer, const int *il3);
     */
    // Entry point registration at startup
    
    static int reg_p_method_infer_update_marshall17();
    // Entry point index lookup
    
    inline static int idx_p_method_infer_update_marshall17() {
      static int epidx = reg_p_method_infer_update_marshall17();
      return epidx;
    }

    
    inline static int idx_p_method_infer_update(void (EnzoBlock::*)(int n, const char *buffer, const int *il3) ) {
      return idx_p_method_infer_update_marshall17();
    }


    
    static int p_method_infer_update(int n, const char *buffer, const int *il3) { return idx_p_method_infer_update_marshall17(); }
    
    static void _call_p_method_infer_update_marshall17(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_infer_update_marshall17(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_method_infer_update_marshall17(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_method_infer_update_marshall17(PUP::er &p,void *msg);
    /* DECLS: void p_method_infer_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_method_infer_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_method_infer_exit_void() {
      static int epidx = reg_p_method_infer_exit_void();
      return epidx;
    }

    
    inline static int idx_p_method_infer_exit(void (EnzoBlock::*)() ) {
      return idx_p_method_infer_exit_void();
    }


    
    static int p_method_infer_exit() { return idx_p_method_infer_exit_void(); }
    
    static void _call_p_method_infer_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_infer_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_method_fbnet_update_mesh_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_p_method_fbnet_update_mesh_CkReductionMsg() {
      static int epidx = reg_p_method_fbnet_update_mesh_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_p_method_fbnet_update_mesh(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_p_method_fbnet_update_mesh_CkReductionMsg();
    }


    
    static int p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) { return idx_p_method_fbnet_update_mesh_CkReductionMsg(); }
    
    static void _call_p_method_fbnet_update_mesh_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_fbnet_update_mesh_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_fbnet_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_method_fbnet_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_method_fbnet_exit_void() {
      static int epidx = reg_p_method_fbnet_exit_void();
      return epidx;
    }

    
    inline static int idx_p_method_fbnet_exit(void (EnzoBlock::*)() ) {
      return idx_p_method_fbnet_exit_void();
    }


    
    static int p_method_fbnet_exit() { return idx_p_method_fbnet_exit_void(); }
    
    static void _call_p_method_fbnet_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_fbnet_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
     */
    // Entry point registration at startup
    
    static int reg_p_check_write_first_marshall21();
    // Entry point index lookup
    
    inline static int idx_p_check_write_first_marshall21() {
      static int epidx = reg_p_check_write_first_marshall21();
      return epidx;
    }

    
    inline static int idx_p_check_write_first(void (EnzoBlock::*)(int num_files, const std::string &ordering, const std::string &name_dir) ) {
      return idx_p_check_write_first_marshall21();
    }


    
    static int p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir) { return idx_p_check_write_first_marshall21(); }
    
    static void _call_p_check_write_first_marshall21(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_check_write_first_marshall21(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_check_write_first_marshall21(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_check_write_first_marshall21(PUP::er &p,void *msg);
    /* DECLS: void p_check_write_next(int num_files, const std::string &ordering);
     */
    // Entry point registration at startup
    
    static int reg_p_check_write_next_marshall22();
    // Entry point index lookup
    
    inline static int idx_p_check_write_next_marshall22() {
      static int epidx = reg_p_check_write_next_marshall22();
      return epidx;
    }

    
    inline static int idx_p_check_write_next(void (EnzoBlock::*)(int num_files, const std::string &ordering) ) {
      return idx_p_check_write_next_marshall22();
    }


    
    static int p_check_write_next(int num_files, const std::string &ordering) { return idx_p_check_write_next_marshall22(); }
    
    static void _call_p_check_write_next_marshall22(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_check_write_next_marshall22(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_check_write_next_marshall22(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_check_write_next_marshall22(PUP::er &p,void *msg);
    /* DECLS: void p_check_done();
     */
    // Entry point registration at startup
    
    static int reg_p_check_done_void();
    // Entry point index lookup
    
    inline static int idx_p_check_done_void() {
      static int epidx = reg_p_check_done_void();
      return epidx;
    }

    
    inline static int idx_p_check_done(void (EnzoBlock::*)() ) {
      return idx_p_check_done_void();
    }


    
    static int p_check_done() { return idx_p_check_done_void(); }
    
    static void _call_p_check_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_check_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_restart_refine(const int *ic3, int io_reader);
     */
    // Entry point registration at startup
    
    static int reg_p_restart_refine_marshall24();
    // Entry point index lookup
    
    inline static int idx_p_restart_refine_marshall24() {
      static int epidx = reg_p_restart_refine_marshall24();
      return epidx;
    }

    
    inline static int idx_p_restart_refine(void (EnzoBlock::*)(const int *ic3, int io_reader) ) {
      return idx_p_restart_refine_marshall24();
    }


    
    static int p_restart_refine(const int *ic3, int io_reader) { return idx_p_restart_refine_marshall24(); }
    
    static void _call_p_restart_refine_marshall24(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_refine_marshall24(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_restart_refine_marshall24(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_restart_refine_marshall24(PUP::er &p,void *msg);
    /* DECLS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_restart_set_data_EnzoMsgCheck();
    // Entry point index lookup
    
    inline static int idx_p_restart_set_data_EnzoMsgCheck() {
      static int epidx = reg_p_restart_set_data_EnzoMsgCheck();
      return epidx;
    }

    
    inline static int idx_p_restart_set_data(void (EnzoBlock::*)(EnzoMsgCheck* impl_msg) ) {
      return idx_p_restart_set_data_EnzoMsgCheck();
    }


    
    static int p_restart_set_data(EnzoMsgCheck* impl_msg) { return idx_p_restart_set_data_EnzoMsgCheck(); }
    
    static void _call_p_restart_set_data_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_set_data_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    /* DECLS: void p_restart_done();
     */
    // Entry point registration at startup
    
    static int reg_p_restart_done_void();
    // Entry point index lookup
    
    inline static int idx_p_restart_done_void() {
      static int epidx = reg_p_restart_done_void();
      return epidx;
    }

    
    inline static int idx_p_restart_done(void (EnzoBlock::*)() ) {
      return idx_p_restart_done_void();
    }


    
    static int p_restart_done() { return idx_p_restart_done_void(); }
    
    static void _call_p_restart_done_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_restart_done_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_method_accretion_end();
     */
    // Entry point registration at startup
    
    static int reg_p_method_accretion_end_void();
    // Entry point index lookup
    
    inline static int idx_p_method_accretion_end_void() {
      static int epidx = reg_p_method_accretion_end_void();
      return epidx;
    }

    
    inline static int idx_p_method_accretion_end(void (EnzoBlock::*)() ) {
      return idx_p_method_accretion_end_void();
    }


    
    static int p_method_accretion_end() { return idx_p_method_accretion_end_void(); }
    
    static void _call_p_method_accretion_end_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_method_accretion_end_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_cg_matvec();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_cg_matvec_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_cg_matvec_void() {
      static int epidx = reg_p_solver_cg_matvec_void();
      return epidx;
    }

    
    inline static int idx_p_solver_cg_matvec(void (EnzoBlock::*)() ) {
      return idx_p_solver_cg_matvec_void();
    }


    
    static int p_solver_cg_matvec() { return idx_p_solver_cg_matvec_void(); }
    
    static void _call_p_solver_cg_matvec_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_cg_matvec_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_cg_loop_0a_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_cg_loop_0a_CkReductionMsg() {
      static int epidx = reg_r_solver_cg_loop_0a_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_cg_loop_0a(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_cg_loop_0a_CkReductionMsg();
    }


    
    static int r_solver_cg_loop_0a(CkReductionMsg* impl_msg) { return idx_r_solver_cg_loop_0a_CkReductionMsg(); }
    
    static void _call_r_solver_cg_loop_0a_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_cg_loop_0a_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_cg_loop_0b_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_cg_loop_0b_CkReductionMsg() {
      static int epidx = reg_r_solver_cg_loop_0b_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_cg_loop_0b(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_cg_loop_0b_CkReductionMsg();
    }


    
    static int r_solver_cg_loop_0b(CkReductionMsg* impl_msg) { return idx_r_solver_cg_loop_0b_CkReductionMsg(); }
    
    static void _call_r_solver_cg_loop_0b_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_cg_loop_0b_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_cg_shift_1_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_cg_shift_1_CkReductionMsg() {
      static int epidx = reg_r_solver_cg_shift_1_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_cg_shift_1(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_cg_shift_1_CkReductionMsg();
    }


    
    static int r_solver_cg_shift_1(CkReductionMsg* impl_msg) { return idx_r_solver_cg_shift_1_CkReductionMsg(); }
    
    static void _call_r_solver_cg_shift_1_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_cg_shift_1_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_cg_loop_2();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_cg_loop_2_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_cg_loop_2_void() {
      static int epidx = reg_p_solver_cg_loop_2_void();
      return epidx;
    }

    
    inline static int idx_p_solver_cg_loop_2(void (EnzoBlock::*)() ) {
      return idx_p_solver_cg_loop_2_void();
    }


    
    static int p_solver_cg_loop_2() { return idx_p_solver_cg_loop_2_void(); }
    
    static void _call_p_solver_cg_loop_2_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_cg_loop_2_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_cg_loop_3_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_cg_loop_3_CkReductionMsg() {
      static int epidx = reg_r_solver_cg_loop_3_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_cg_loop_3(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_cg_loop_3_CkReductionMsg();
    }


    
    static int r_solver_cg_loop_3(CkReductionMsg* impl_msg) { return idx_r_solver_cg_loop_3_CkReductionMsg(); }
    
    static void _call_r_solver_cg_loop_3_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_cg_loop_3_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_cg_loop_5_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_cg_loop_5_CkReductionMsg() {
      static int epidx = reg_r_solver_cg_loop_5_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_cg_loop_5(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_cg_loop_5_CkReductionMsg();
    }


    
    static int r_solver_cg_loop_5(CkReductionMsg* impl_msg) { return idx_r_solver_cg_loop_5_CkReductionMsg(); }
    
    static void _call_r_solver_cg_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_cg_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_start_1_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_start_1_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_start_1_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_start_1(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_start_1_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_start_1_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_start_1_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_start_1_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_start_3_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_start_3_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_start_3_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_start_3(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_start_3_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_start_3_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_start_3_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_start_3_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_loop_5_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_loop_5_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_loop_5_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_loop_5(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_loop_5_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_loop_5_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_loop_5_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_loop_11_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_loop_11_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_loop_11_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_loop_11(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_loop_11_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_loop_11_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_loop_11_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_loop_11_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_loop_13_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_loop_13_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_loop_13_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_loop_13(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_loop_13_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_loop_13_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_loop_13_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_loop_13_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_bicgstab_loop_15_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_bicgstab_loop_15_CkReductionMsg() {
      static int epidx = reg_r_solver_bicgstab_loop_15_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_bicgstab_loop_15(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_bicgstab_loop_15_CkReductionMsg();
    }


    
    static int r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) { return idx_r_solver_bicgstab_loop_15_CkReductionMsg(); }
    
    static void _call_r_solver_bicgstab_loop_15_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_bicgstab_loop_15_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_bicgstab_loop_2();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_bicgstab_loop_2_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_bicgstab_loop_2_void() {
      static int epidx = reg_p_solver_bicgstab_loop_2_void();
      return epidx;
    }

    
    inline static int idx_p_solver_bicgstab_loop_2(void (EnzoBlock::*)() ) {
      return idx_p_solver_bicgstab_loop_2_void();
    }


    
    static int p_solver_bicgstab_loop_2() { return idx_p_solver_bicgstab_loop_2_void(); }
    
    static void _call_p_solver_bicgstab_loop_2_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_bicgstab_loop_2_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_bicgstab_loop_3();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_bicgstab_loop_3_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_bicgstab_loop_3_void() {
      static int epidx = reg_p_solver_bicgstab_loop_3_void();
      return epidx;
    }

    
    inline static int idx_p_solver_bicgstab_loop_3(void (EnzoBlock::*)() ) {
      return idx_p_solver_bicgstab_loop_3_void();
    }


    
    static int p_solver_bicgstab_loop_3() { return idx_p_solver_bicgstab_loop_3_void(); }
    
    static void _call_p_solver_bicgstab_loop_3_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_bicgstab_loop_3_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_bicgstab_loop_8();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_bicgstab_loop_8_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_bicgstab_loop_8_void() {
      static int epidx = reg_p_solver_bicgstab_loop_8_void();
      return epidx;
    }

    
    inline static int idx_p_solver_bicgstab_loop_8(void (EnzoBlock::*)() ) {
      return idx_p_solver_bicgstab_loop_8_void();
    }


    
    static int p_solver_bicgstab_loop_8() { return idx_p_solver_bicgstab_loop_8_void(); }
    
    static void _call_p_solver_bicgstab_loop_8_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_bicgstab_loop_8_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_bicgstab_loop_9();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_bicgstab_loop_9_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_bicgstab_loop_9_void() {
      static int epidx = reg_p_solver_bicgstab_loop_9_void();
      return epidx;
    }

    
    inline static int idx_p_solver_bicgstab_loop_9(void (EnzoBlock::*)() ) {
      return idx_p_solver_bicgstab_loop_9_void();
    }


    
    static int p_solver_bicgstab_loop_9() { return idx_p_solver_bicgstab_loop_9_void(); }
    
    static void _call_p_solver_bicgstab_loop_9_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_bicgstab_loop_9_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
     */
    // Entry point registration at startup
    
    static int reg_p_dot_recv_parent_marshall45();
    // Entry point index lookup
    
    inline static int idx_p_dot_recv_parent_marshall45() {
      static int epidx = reg_p_dot_recv_parent_marshall45();
      return epidx;
    }

    
    inline static int idx_p_dot_recv_parent(void (EnzoBlock::*)(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter) ) {
      return idx_p_dot_recv_parent_marshall45();
    }


    
    static int p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter) { return idx_p_dot_recv_parent_marshall45(); }
    
    static void _call_p_dot_recv_parent_marshall45(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_dot_recv_parent_marshall45(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_dot_recv_parent_marshall45(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_dot_recv_parent_marshall45(PUP::er &p,void *msg);
    /* DECLS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
     */
    // Entry point registration at startup
    
    static int reg_p_dot_recv_children_marshall46();
    // Entry point index lookup
    
    inline static int idx_p_dot_recv_children_marshall46() {
      static int epidx = reg_p_dot_recv_children_marshall46();
      return epidx;
    }

    
    inline static int idx_p_dot_recv_children(void (EnzoBlock::*)(int n, const long double *dot, const std::vector<int> &isa, int i_function) ) {
      return idx_p_dot_recv_children_marshall46();
    }


    
    static int p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function) { return idx_p_dot_recv_children_marshall46(); }
    
    static void _call_p_dot_recv_children_marshall46(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_dot_recv_children_marshall46(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_dot_recv_children_marshall46(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_dot_recv_children_marshall46(PUP::er &p,void *msg);
    /* DECLS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_solver_dd_restrict_recv_FieldMsg();
    // Entry point index lookup
    
    inline static int idx_p_solver_dd_restrict_recv_FieldMsg() {
      static int epidx = reg_p_solver_dd_restrict_recv_FieldMsg();
      return epidx;
    }

    
    inline static int idx_p_solver_dd_restrict_recv(void (EnzoBlock::*)(FieldMsg* impl_msg) ) {
      return idx_p_solver_dd_restrict_recv_FieldMsg();
    }


    
    static int p_solver_dd_restrict_recv(FieldMsg* impl_msg) { return idx_p_solver_dd_restrict_recv_FieldMsg(); }
    
    static void _call_p_solver_dd_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_dd_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_solver_dd_prolong_recv_FieldMsg();
    // Entry point index lookup
    
    inline static int idx_p_solver_dd_prolong_recv_FieldMsg() {
      static int epidx = reg_p_solver_dd_prolong_recv_FieldMsg();
      return epidx;
    }

    
    inline static int idx_p_solver_dd_prolong_recv(void (EnzoBlock::*)(FieldMsg* impl_msg) ) {
      return idx_p_solver_dd_prolong_recv_FieldMsg();
    }


    
    static int p_solver_dd_prolong_recv(FieldMsg* impl_msg) { return idx_p_solver_dd_prolong_recv_FieldMsg(); }
    
    static void _call_p_solver_dd_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_dd_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_dd_solve_coarse();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_dd_solve_coarse_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_dd_solve_coarse_void() {
      static int epidx = reg_p_solver_dd_solve_coarse_void();
      return epidx;
    }

    
    inline static int idx_p_solver_dd_solve_coarse(void (EnzoBlock::*)() ) {
      return idx_p_solver_dd_solve_coarse_void();
    }


    
    static int p_solver_dd_solve_coarse() { return idx_p_solver_dd_solve_coarse_void(); }
    
    static void _call_p_solver_dd_solve_coarse_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_dd_solve_coarse_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_dd_solve_domain();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_dd_solve_domain_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_dd_solve_domain_void() {
      static int epidx = reg_p_solver_dd_solve_domain_void();
      return epidx;
    }

    
    inline static int idx_p_solver_dd_solve_domain(void (EnzoBlock::*)() ) {
      return idx_p_solver_dd_solve_domain_void();
    }


    
    static int p_solver_dd_solve_domain() { return idx_p_solver_dd_solve_domain_void(); }
    
    static void _call_p_solver_dd_solve_domain_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_dd_solve_domain_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_dd_last_smooth();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_dd_last_smooth_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_dd_last_smooth_void() {
      static int epidx = reg_p_solver_dd_last_smooth_void();
      return epidx;
    }

    
    inline static int idx_p_solver_dd_last_smooth(void (EnzoBlock::*)() ) {
      return idx_p_solver_dd_last_smooth_void();
    }


    
    static int p_solver_dd_last_smooth() { return idx_p_solver_dd_last_smooth_void(); }
    
    static void _call_p_solver_dd_last_smooth_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_dd_last_smooth_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_dd_barrier_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_dd_barrier_CkReductionMsg() {
      static int epidx = reg_r_solver_dd_barrier_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_dd_barrier(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_dd_barrier_CkReductionMsg();
    }


    
    static int r_solver_dd_barrier(CkReductionMsg* impl_msg) { return idx_r_solver_dd_barrier_CkReductionMsg(); }
    
    static void _call_r_solver_dd_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_dd_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_dd_end(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_dd_end_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_dd_end_CkReductionMsg() {
      static int epidx = reg_r_solver_dd_end_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_dd_end(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_dd_end_CkReductionMsg();
    }


    
    static int r_solver_dd_end(CkReductionMsg* impl_msg) { return idx_r_solver_dd_end_CkReductionMsg(); }
    
    static void _call_r_solver_dd_end_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_dd_end_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_jacobi_continue();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_jacobi_continue_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_jacobi_continue_void() {
      static int epidx = reg_p_solver_jacobi_continue_void();
      return epidx;
    }

    
    inline static int idx_p_solver_jacobi_continue(void (EnzoBlock::*)() ) {
      return idx_p_solver_jacobi_continue_void();
    }


    
    static int p_solver_jacobi_continue() { return idx_p_solver_jacobi_continue_void(); }
    
    static void _call_p_solver_jacobi_continue_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_jacobi_continue_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_restrict();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_restrict_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_restrict_void() {
      static int epidx = reg_p_solver_mg0_restrict_void();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_restrict(void (EnzoBlock::*)() ) {
      return idx_p_solver_mg0_restrict_void();
    }


    
    static int p_solver_mg0_restrict() { return idx_p_solver_mg0_restrict_void(); }
    
    static void _call_p_solver_mg0_restrict_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_restrict_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_solve_coarse();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_solve_coarse_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_solve_coarse_void() {
      static int epidx = reg_p_solver_mg0_solve_coarse_void();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_solve_coarse(void (EnzoBlock::*)() ) {
      return idx_p_solver_mg0_solve_coarse_void();
    }


    
    static int p_solver_mg0_solve_coarse() { return idx_p_solver_mg0_solve_coarse_void(); }
    
    static void _call_p_solver_mg0_solve_coarse_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_solve_coarse_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_post_smooth();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_post_smooth_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_post_smooth_void() {
      static int epidx = reg_p_solver_mg0_post_smooth_void();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_post_smooth(void (EnzoBlock::*)() ) {
      return idx_p_solver_mg0_post_smooth_void();
    }


    
    static int p_solver_mg0_post_smooth() { return idx_p_solver_mg0_post_smooth_void(); }
    
    static void _call_p_solver_mg0_post_smooth_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_post_smooth_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_last_smooth();
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_last_smooth_void();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_last_smooth_void() {
      static int epidx = reg_p_solver_mg0_last_smooth_void();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_last_smooth(void (EnzoBlock::*)() ) {
      return idx_p_solver_mg0_last_smooth_void();
    }


    
    static int p_solver_mg0_last_smooth() { return idx_p_solver_mg0_last_smooth_void(); }
    
    static void _call_p_solver_mg0_last_smooth_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_last_smooth_void(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_mg0_begin_solve_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_mg0_begin_solve_CkReductionMsg() {
      static int epidx = reg_r_solver_mg0_begin_solve_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_mg0_begin_solve(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_mg0_begin_solve_CkReductionMsg();
    }


    
    static int r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) { return idx_r_solver_mg0_begin_solve_CkReductionMsg(); }
    
    static void _call_r_solver_mg0_begin_solve_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_mg0_begin_solve_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_r_solver_mg0_barrier_CkReductionMsg();
    // Entry point index lookup
    
    inline static int idx_r_solver_mg0_barrier_CkReductionMsg() {
      static int epidx = reg_r_solver_mg0_barrier_CkReductionMsg();
      return epidx;
    }

    
    inline static int idx_r_solver_mg0_barrier(void (EnzoBlock::*)(CkReductionMsg* impl_msg) ) {
      return idx_r_solver_mg0_barrier_CkReductionMsg();
    }


    
    static int r_solver_mg0_barrier(CkReductionMsg* impl_msg) { return idx_r_solver_mg0_barrier_CkReductionMsg(); }
    
    static void _call_r_solver_mg0_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_r_solver_mg0_barrier_CkReductionMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_prolong_recv_FieldMsg();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_prolong_recv_FieldMsg() {
      static int epidx = reg_p_solver_mg0_prolong_recv_FieldMsg();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_prolong_recv(void (EnzoBlock::*)(FieldMsg* impl_msg) ) {
      return idx_p_solver_mg0_prolong_recv_FieldMsg();
    }


    
    static int p_solver_mg0_prolong_recv(FieldMsg* impl_msg) { return idx_p_solver_mg0_prolong_recv_FieldMsg(); }
    
    static void _call_p_solver_mg0_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_prolong_recv_FieldMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_solver_mg0_restrict_recv_FieldMsg();
    // Entry point index lookup
    
    inline static int idx_p_solver_mg0_restrict_recv_FieldMsg() {
      static int epidx = reg_p_solver_mg0_restrict_recv_FieldMsg();
      return epidx;
    }

    
    inline static int idx_p_solver_mg0_restrict_recv(void (EnzoBlock::*)(FieldMsg* impl_msg) ) {
      return idx_p_solver_mg0_restrict_recv_FieldMsg();
    }


    
    static int p_solver_mg0_restrict_recv(FieldMsg* impl_msg) { return idx_p_solver_mg0_restrict_recv_FieldMsg(); }
    
    static void _call_p_solver_mg0_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_solver_mg0_restrict_recv_FieldMsg(void* impl_msg, void* impl_obj);
    /* DECLS: EnzoBlock(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_EnzoBlock_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_EnzoBlock_CkMigrateMessage() {
      static int epidx = reg_EnzoBlock_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_EnzoBlock_CkMigrateMessage(); }
    
    static void _call_EnzoBlock_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoBlock_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_EnzoBlock : public CProxyElement_Block{
  public:
    typedef EnzoBlock local_t;
    typedef CkIndex_EnzoBlock index_t;
    typedef CProxy_EnzoBlock proxy_t;
    typedef CProxyElement_EnzoBlock element_t;
    typedef CProxySection_EnzoBlock section_t;

    using array_index_t = CkArrayIndexIndex;

    /* TRAM aggregators */

    CProxyElement_EnzoBlock(void) {
    }
    CProxyElement_EnzoBlock(const ArrayElement *e) : CProxyElement_Block(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_Block::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_Block::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_Block::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_Block::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_Block::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_Block::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_Block::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_Block::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_Block::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_Block::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_Block::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_Block::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_Block::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_Block::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_Block::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_Block::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_Block::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_Block::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_Block::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_Block::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_Block::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_Block::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_Block::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_Block::ckGetIndex(); }

    EnzoBlock *ckLocal(void) const
    { return (EnzoBlock *)CProxyElement_Block::ckLocal(); }

    CProxyElement_EnzoBlock(const CkArrayID &aid,const CkArrayIndexIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_Block(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_EnzoBlock(const CkArrayID &aid,const CkArrayIndexIndex &idx)
        :CProxyElement_Block(aid,idx)
    {
}

    CProxyElement_EnzoBlock(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_Block(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_EnzoBlock(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_Block(aid,idx)
    {
}
/* DECLS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
    
    void insert(const process_type &ip_source, const MsgType &msg_type, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
    
    void p_set_msg_check(EnzoMsgCheck* impl_msg) ;

/* DECLS: EnzoBlock();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_method_feedback_starss_end();
 */
    
    void p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_solve_transport_eqn();
 */
    
    void p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
    
    void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
    
    void r_method_turbulence_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
    
    void p_initial_hdf5_recv(MsgInitial* impl_msg) ;

/* DECLS: void p_method_balance_migrate();
 */
    
    void p_method_balance_migrate(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_balance_done();
 */
    
    void p_method_balance_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_continue();
 */
    
    void p_method_gravity_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_end();
 */
    
    void p_method_gravity_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
    
    void p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_count_arrays(int count);
 */
    
    void p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_request_data(const int *il3);
 */
    
    void p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
    
    void p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_exit();
 */
    
    void p_method_infer_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
    
    void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_fbnet_exit();
 */
    
    void p_method_fbnet_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
    
    void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_next(int num_files, const std::string &ordering);
 */
    
    void p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_refine(const int *ic3, int io_reader);
 */
    
    void p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
    
    void p_restart_set_data(EnzoMsgCheck* impl_msg) ;

/* DECLS: void p_restart_done();
 */
    
    void p_restart_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_accretion_end();
 */
    
    void p_method_accretion_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_cg_matvec();
 */
    
    void p_solver_cg_matvec(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0a(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0b(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_shift_1(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_cg_loop_2();
 */
    
    void p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_bicgstab_loop_2();
 */
    
    void p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_3();
 */
    
    void p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_8();
 */
    
    void p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_9();
 */
    
    void p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
    
    void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
    
    void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_solve_coarse();
 */
    
    void p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_solve_domain();
 */
    
    void p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_last_smooth();
 */
    
    void p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_jacobi_continue();
 */
    
    void p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_restrict();
 */
    
    void p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_solve_coarse();
 */
    
    void p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_post_smooth();
 */
    
    void p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_last_smooth();
 */
    
    void p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: EnzoBlock(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_EnzoBlock : public CProxy_Block{
  public:
    typedef EnzoBlock local_t;
    typedef CkIndex_EnzoBlock index_t;
    typedef CProxy_EnzoBlock proxy_t;
    typedef CProxyElement_EnzoBlock element_t;
    typedef CProxySection_EnzoBlock section_t;

    using array_index_t = CkArrayIndexIndex;
    CProxy_EnzoBlock(void) {
    }
    CProxy_EnzoBlock(const ArrayElement *e) : CProxy_Block(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_Block::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_Block::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_Block::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_Block::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_Block::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_Block::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_Block::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_Block::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_Block::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_Block::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_Block::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_Block::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_Block::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_Block::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_Block::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_Block::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_Block::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_Block::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_Block::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_Block::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_EnzoBlock operator [] (const CkArrayIndexIndex &idx) const
    { return CProxyElement_EnzoBlock(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_EnzoBlock operator() (const CkArrayIndexIndex &idx) const
    { return CProxyElement_EnzoBlock(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxy_EnzoBlock(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_Block(aid,CK_DELCTOR_ARGS) {}
    CProxy_EnzoBlock(const CkArrayID &aid) 
        :CProxy_Block(aid) {}
/* DECLS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
    
    static CkArrayID ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const process_type &ip_source, const MsgType &msg_type, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
    
    void p_set_msg_check(EnzoMsgCheck* impl_msg) ;

/* DECLS: EnzoBlock();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_method_feedback_starss_end();
 */
    
    void p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_solve_transport_eqn();
 */
    
    void p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
    
    void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
    
    void r_method_turbulence_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
    
    void p_initial_hdf5_recv(MsgInitial* impl_msg) ;

/* DECLS: void p_method_balance_migrate();
 */
    
    void p_method_balance_migrate(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_balance_done();
 */
    
    void p_method_balance_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_continue();
 */
    
    void p_method_gravity_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_end();
 */
    
    void p_method_gravity_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
    
    void p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_count_arrays(int count);
 */
    
    void p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_request_data(const int *il3);
 */
    
    void p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
    
    void p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_exit();
 */
    
    void p_method_infer_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
    
    void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_fbnet_exit();
 */
    
    void p_method_fbnet_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
    
    void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_next(int num_files, const std::string &ordering);
 */
    
    void p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_refine(const int *ic3, int io_reader);
 */
    
    void p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
    
    void p_restart_set_data(EnzoMsgCheck* impl_msg) ;

/* DECLS: void p_restart_done();
 */
    
    void p_restart_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_accretion_end();
 */
    
    void p_method_accretion_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_cg_matvec();
 */
    
    void p_solver_cg_matvec(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0a(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0b(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_shift_1(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_cg_loop_2();
 */
    
    void p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_bicgstab_loop_2();
 */
    
    void p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_3();
 */
    
    void p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_8();
 */
    
    void p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_9();
 */
    
    void p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
    
    void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
    
    void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_solve_coarse();
 */
    
    void p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_solve_domain();
 */
    
    void p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_last_smooth();
 */
    
    void p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_jacobi_continue();
 */
    
    void p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_restrict();
 */
    
    void p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_solve_coarse();
 */
    
    void p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_post_smooth();
 */
    
    void p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_last_smooth();
 */
    
    void p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: EnzoBlock(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_EnzoBlock : public CProxySection_Block{
  public:
    typedef EnzoBlock local_t;
    typedef CkIndex_EnzoBlock index_t;
    typedef CProxy_EnzoBlock proxy_t;
    typedef CProxyElement_EnzoBlock element_t;
    typedef CProxySection_EnzoBlock section_t;

    using array_index_t = CkArrayIndexIndex;
    CProxySection_EnzoBlock(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_Block::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_Block::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_Block::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_Block::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_Block::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_Block::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_Block::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_Block::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_Block::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_Block::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_Block::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_Block::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_Block::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_Block::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_Block::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_Block::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_Block::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_Block::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_Block::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_Block::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_Block::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_Block::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_Block::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_Block::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_Block::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_Block::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_Block::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_Block::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_Block::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_Block::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_EnzoBlock operator [] (const CkArrayIndexIndex &idx) const
        {return CProxyElement_EnzoBlock(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_EnzoBlock operator() (const CkArrayIndexIndex &idx) const
        {return CProxyElement_EnzoBlock(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_EnzoBlock(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_Block(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoBlock(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_Block(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoBlock(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_Block(aid,elems,nElems, factor) {}
    CProxySection_EnzoBlock(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_Block(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_EnzoBlock(const CkSectionID &sid)  
        :CProxySection_Block(sid) { ckAutoDelegate(); }
    CProxySection_EnzoBlock(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_Block(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoBlock(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_Block(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoBlock(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_Block(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_EnzoBlock(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_Block(aid,elems) { ckAutoDelegate(); }
    CProxySection_EnzoBlock(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_Block(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_EnzoBlock(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_Block(aid,elems, factor) { ckAutoDelegate(); }
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
      CProxySection_Block::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_Block::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_Block::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: EnzoBlock(const process_type &ip_source, const MsgType &msg_type);
 */
    

/* DECLS: void p_set_msg_refine(MsgRefine* impl_msg);
 */
    
    void p_set_msg_refine(MsgRefine* impl_msg) ;

/* DECLS: void p_set_msg_check(EnzoMsgCheck* impl_msg);
 */
    
    void p_set_msg_check(EnzoMsgCheck* impl_msg) ;

/* DECLS: EnzoBlock();
 */
    

/* DECLS: void p_method_feedback_starss_end();
 */
    
    void p_method_feedback_starss_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_solve_transport_eqn();
 */
    
    void p_method_m1_closure_solve_transport_eqn(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg);
 */
    
    void p_method_m1_closure_set_global_averages(CkReductionMsg* impl_msg) ;

/* DECLS: void r_method_turbulence_end(CkReductionMsg* impl_msg);
 */
    
    void r_method_turbulence_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_initial_hdf5_recv(MsgInitial* impl_msg);
 */
    
    void p_initial_hdf5_recv(MsgInitial* impl_msg) ;

/* DECLS: void p_method_balance_migrate();
 */
    
    void p_method_balance_migrate(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_balance_done();
 */
    
    void p_method_balance_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_continue();
 */
    
    void p_method_gravity_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_gravity_end();
 */
    
    void p_method_gravity_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_merge_masks(int n, const char *mask, const int *ic3);
 */
    
    void p_method_infer_merge_masks(int n, const char *mask, const int *ic3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_count_arrays(int count);
 */
    
    void p_method_infer_count_arrays(int count, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_request_data(const int *il3);
 */
    
    void p_method_infer_request_data(const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_update(int n, const char *buffer, const int *il3);
 */
    
    void p_method_infer_update(int n, const char *buffer, const int *il3, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_infer_exit();
 */
    
    void p_method_infer_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg);
 */
    
    void p_method_fbnet_update_mesh(CkReductionMsg* impl_msg) ;

/* DECLS: void p_method_fbnet_exit();
 */
    
    void p_method_fbnet_exit(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir);
 */
    
    void p_check_write_first(int num_files, const std::string &ordering, const std::string &name_dir, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_write_next(int num_files, const std::string &ordering);
 */
    
    void p_check_write_next(int num_files, const std::string &ordering, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_check_done();
 */
    
    void p_check_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_refine(const int *ic3, int io_reader);
 */
    
    void p_restart_refine(const int *ic3, int io_reader, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_restart_set_data(EnzoMsgCheck* impl_msg);
 */
    
    void p_restart_set_data(EnzoMsgCheck* impl_msg) ;

/* DECLS: void p_restart_done();
 */
    
    void p_restart_done(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_method_accretion_end();
 */
    
    void p_method_accretion_end(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_cg_matvec();
 */
    
    void p_solver_cg_matvec(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_0a(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0a(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_0b(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_0b(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_shift_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_shift_1(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_cg_loop_2();
 */
    
    void p_solver_cg_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_cg_loop_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_cg_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_cg_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_1(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_start_3(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_5(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_11(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_13(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg);
 */
    
    void r_solver_bicgstab_loop_15(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_bicgstab_loop_2();
 */
    
    void p_solver_bicgstab_loop_2(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_3();
 */
    
    void p_solver_bicgstab_loop_3(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_8();
 */
    
    void p_solver_bicgstab_loop_8(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_bicgstab_loop_9();
 */
    
    void p_solver_bicgstab_loop_9(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter);
 */
    
    void p_dot_recv_parent(int n, const long double *dot, const std::vector<int> &isa, int i_function, int iter, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function);
 */
    
    void p_dot_recv_children(int n, const long double *dot, const std::vector<int> &isa, int i_function, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_dd_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_dd_solve_coarse();
 */
    
    void p_solver_dd_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_solve_domain();
 */
    
    void p_solver_dd_solve_domain(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_dd_last_smooth();
 */
    
    void p_solver_dd_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_dd_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_dd_end(CkReductionMsg* impl_msg);
 */
    
    void r_solver_dd_end(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_jacobi_continue();
 */
    
    void p_solver_jacobi_continue(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_restrict();
 */
    
    void p_solver_mg0_restrict(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_solve_coarse();
 */
    
    void p_solver_mg0_solve_coarse(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_post_smooth();
 */
    
    void p_solver_mg0_post_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_solver_mg0_last_smooth();
 */
    
    void p_solver_mg0_last_smooth(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_begin_solve(CkReductionMsg* impl_msg) ;

/* DECLS: void r_solver_mg0_barrier(CkReductionMsg* impl_msg);
 */
    
    void r_solver_mg0_barrier(CkReductionMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_prolong_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_prolong_recv(FieldMsg* impl_msg) ;

/* DECLS: void p_solver_mg0_restrict_recv(FieldMsg* impl_msg);
 */
    
    void p_solver_mg0_restrict_recv(FieldMsg* impl_msg) ;

/* DECLS: EnzoBlock(CkMigrateMessage* impl_msg);
 */

};
#define EnzoBlock_SDAG_CODE 
typedef CBaseT1<Block, CProxy_EnzoBlock>CBase_EnzoBlock;

/* DECLS: array IoEnzoReader: IoReader{
IoEnzoReader();
void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
void p_create_level(int level);
void p_init_level(int level);
void p_block_created();
void p_block_ready();
IoEnzoReader(CkMigrateMessage* impl_msg);
};
 */
 class IoEnzoReader;
 class CkIndex_IoEnzoReader;
 class CProxy_IoEnzoReader;
 class CProxyElement_IoEnzoReader;
 class CProxySection_IoEnzoReader;
/* --------------- index object ------------------ */
class CkIndex_IoEnzoReader:public CkIndex_IoReader{
  public:
    typedef IoEnzoReader local_t;
    typedef CkIndex_IoEnzoReader index_t;
    typedef CProxy_IoEnzoReader proxy_t;
    typedef CProxyElement_IoEnzoReader element_t;
    typedef CProxySection_IoEnzoReader section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: IoEnzoReader();
     */
    // Entry point registration at startup
    
    static int reg_IoEnzoReader_void();
    // Entry point index lookup
    
    inline static int idx_IoEnzoReader_void() {
      static int epidx = reg_IoEnzoReader_void();
      return epidx;
    }

    
    static int ckNew() { return idx_IoEnzoReader_void(); }
    
    static void _call_IoEnzoReader_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoEnzoReader_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
     */
    // Entry point registration at startup
    
    static int reg_p_init_root_marshall2();
    // Entry point index lookup
    
    inline static int idx_p_init_root_marshall2() {
      static int epidx = reg_p_init_root_marshall2();
      return epidx;
    }

    
    inline static int idx_p_init_root(void (IoEnzoReader::*)(const std::string &impl_noname_0, const std::string &impl_noname_1, int level) ) {
      return idx_p_init_root_marshall2();
    }


    
    static int p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level) { return idx_p_init_root_marshall2(); }
    
    static void _call_p_init_root_marshall2(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_init_root_marshall2(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_init_root_marshall2(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_init_root_marshall2(PUP::er &p,void *msg);
    /* DECLS: void p_create_level(int level);
     */
    // Entry point registration at startup
    
    static int reg_p_create_level_marshall3();
    // Entry point index lookup
    
    inline static int idx_p_create_level_marshall3() {
      static int epidx = reg_p_create_level_marshall3();
      return epidx;
    }

    
    inline static int idx_p_create_level(void (IoEnzoReader::*)(int level) ) {
      return idx_p_create_level_marshall3();
    }


    
    static int p_create_level(int level) { return idx_p_create_level_marshall3(); }
    
    static void _call_p_create_level_marshall3(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_create_level_marshall3(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_create_level_marshall3(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_create_level_marshall3(PUP::er &p,void *msg);
    /* DECLS: void p_init_level(int level);
     */
    // Entry point registration at startup
    
    static int reg_p_init_level_marshall4();
    // Entry point index lookup
    
    inline static int idx_p_init_level_marshall4() {
      static int epidx = reg_p_init_level_marshall4();
      return epidx;
    }

    
    inline static int idx_p_init_level(void (IoEnzoReader::*)(int level) ) {
      return idx_p_init_level_marshall4();
    }


    
    static int p_init_level(int level) { return idx_p_init_level_marshall4(); }
    
    static void _call_p_init_level_marshall4(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_init_level_marshall4(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_init_level_marshall4(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_init_level_marshall4(PUP::er &p,void *msg);
    /* DECLS: void p_block_created();
     */
    // Entry point registration at startup
    
    static int reg_p_block_created_void();
    // Entry point index lookup
    
    inline static int idx_p_block_created_void() {
      static int epidx = reg_p_block_created_void();
      return epidx;
    }

    
    inline static int idx_p_block_created(void (IoEnzoReader::*)() ) {
      return idx_p_block_created_void();
    }


    
    static int p_block_created() { return idx_p_block_created_void(); }
    
    static void _call_p_block_created_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_block_created_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_block_ready();
     */
    // Entry point registration at startup
    
    static int reg_p_block_ready_void();
    // Entry point index lookup
    
    inline static int idx_p_block_ready_void() {
      static int epidx = reg_p_block_ready_void();
      return epidx;
    }

    
    inline static int idx_p_block_ready(void (IoEnzoReader::*)() ) {
      return idx_p_block_ready_void();
    }


    
    static int p_block_ready() { return idx_p_block_ready_void(); }
    
    static void _call_p_block_ready_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_block_ready_void(void* impl_msg, void* impl_obj);
    /* DECLS: IoEnzoReader(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_IoEnzoReader_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_IoEnzoReader_CkMigrateMessage() {
      static int epidx = reg_IoEnzoReader_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_IoEnzoReader_CkMigrateMessage(); }
    
    static void _call_IoEnzoReader_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoEnzoReader_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_IoEnzoReader : public CProxyElement_IoReader{
  public:
    typedef IoEnzoReader local_t;
    typedef CkIndex_IoEnzoReader index_t;
    typedef CProxy_IoEnzoReader proxy_t;
    typedef CProxyElement_IoEnzoReader element_t;
    typedef CProxySection_IoEnzoReader section_t;

    using array_index_t = CkArrayIndex1D;

    /* TRAM aggregators */

    CProxyElement_IoEnzoReader(void) {
    }
    CProxyElement_IoEnzoReader(const ArrayElement *e) : CProxyElement_IoReader(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_IoReader::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_IoReader::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_IoReader::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_IoReader::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_IoReader::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_IoReader::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_IoReader::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_IoReader::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_IoReader::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_IoReader::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_IoReader::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_IoReader::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_IoReader::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_IoReader::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_IoReader::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_IoReader::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_IoReader::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IoReader::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IoReader::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_IoReader::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_IoReader::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_IoReader::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_IoReader::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_IoReader::ckGetIndex(); }

    IoEnzoReader *ckLocal(void) const
    { return (IoEnzoReader *)CProxyElement_IoReader::ckLocal(); }

    CProxyElement_IoEnzoReader(const CkArrayID &aid,const CkArrayIndex1D &idx,CK_DELCTOR_PARAM)
        :CProxyElement_IoReader(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoEnzoReader(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_IoReader(aid,idx)
    {
}

    CProxyElement_IoEnzoReader(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_IoReader(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoEnzoReader(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_IoReader(aid,idx)
    {
}
/* DECLS: IoEnzoReader();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
    
    void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_create_level(int level);
 */
    
    void p_create_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_init_level(int level);
 */
    
    void p_init_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_created();
 */
    
    void p_block_created(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_ready();
 */
    
    void p_block_ready(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_IoEnzoReader : public CProxy_IoReader{
  public:
    typedef IoEnzoReader local_t;
    typedef CkIndex_IoEnzoReader index_t;
    typedef CProxy_IoEnzoReader proxy_t;
    typedef CProxyElement_IoEnzoReader element_t;
    typedef CProxySection_IoEnzoReader section_t;

    using array_index_t = CkArrayIndex1D;
    CProxy_IoEnzoReader(void) {
    }
    CProxy_IoEnzoReader(const ArrayElement *e) : CProxy_IoReader(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_IoReader::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_IoReader::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_IoReader::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_IoReader::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_IoReader::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_IoReader::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_IoReader::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_IoReader::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_IoReader::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_IoReader::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_IoReader::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_IoReader::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_IoReader::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_IoReader::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_IoReader::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_IoReader::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_IoReader::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IoReader::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IoReader::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_IoReader::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_IoEnzoReader operator [] (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoEnzoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoEnzoReader operator() (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoEnzoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoEnzoReader operator [] (int idx) const 
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoReader operator () (int idx) const 
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxy_IoEnzoReader(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_IoReader(aid,CK_DELCTOR_ARGS) {}
    CProxy_IoEnzoReader(const CkArrayID &aid) 
        :CProxy_IoReader(aid) {}
/* DECLS: IoEnzoReader();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);
    static CkArrayID ckNew(const int s1, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
    
    void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_create_level(int level);
 */
    
    void p_create_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_init_level(int level);
 */
    
    void p_init_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_created();
 */
    
    void p_block_created(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_ready();
 */
    
    void p_block_ready(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_IoEnzoReader : public CProxySection_IoReader{
  public:
    typedef IoEnzoReader local_t;
    typedef CkIndex_IoEnzoReader index_t;
    typedef CProxy_IoEnzoReader proxy_t;
    typedef CProxyElement_IoEnzoReader element_t;
    typedef CProxySection_IoEnzoReader section_t;

    using array_index_t = CkArrayIndex1D;
    CProxySection_IoEnzoReader(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_IoReader::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_IoReader::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_IoReader::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_IoReader::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_IoReader::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_IoReader::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_IoReader::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_IoReader::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_IoReader::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_IoReader::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_IoReader::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_IoReader::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_IoReader::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_IoReader::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_IoReader::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_IoReader::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_IoReader::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IoReader::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IoReader::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_IoReader::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_IoReader::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_IoReader::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_IoReader::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_IoReader::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_IoReader::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_IoReader::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_IoReader::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_IoReader::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_IoReader::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_IoReader::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_IoEnzoReader operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoReader operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoReader operator [] (int idx) const 
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoReader operator () (int idx) const 
        {return CProxyElement_IoEnzoReader(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
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
    CProxySection_IoEnzoReader(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_IoReader(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoReader(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_IoReader(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoReader(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_IoReader(aid,elems,nElems, factor) {}
    CProxySection_IoEnzoReader(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_IoReader(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_IoEnzoReader(const CkSectionID &sid)  
        :CProxySection_IoReader(sid) { ckAutoDelegate(); }
    CProxySection_IoEnzoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_IoReader(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_IoReader(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_IoReader(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_IoEnzoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_IoReader(aid,elems) { ckAutoDelegate(); }
    CProxySection_IoEnzoReader(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_IoReader(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_IoEnzoReader(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_IoReader(aid,elems, factor) { ckAutoDelegate(); }
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
      CProxySection_IoReader::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_IoReader::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_IoReader::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: IoEnzoReader();
 */
    

/* DECLS: void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level);
 */
    
    void p_init_root(const std::string &impl_noname_0, const std::string &impl_noname_1, int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_create_level(int level);
 */
    
    void p_create_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_init_level(int level);
 */
    
    void p_init_level(int level, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_created();
 */
    
    void p_block_created(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_block_ready();
 */
    
    void p_block_ready(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: IoEnzoReader(CkMigrateMessage* impl_msg);
 */

};
#define IoEnzoReader_SDAG_CODE 
typedef CBaseT1<IoReader, CProxy_IoEnzoReader>CBase_IoEnzoReader;

/* DECLS: array IoEnzoWriter: IoWriter{
IoEnzoWriter();
IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
void p_write(EnzoMsgCheck* impl_msg);
IoEnzoWriter(CkMigrateMessage* impl_msg);
};
 */
 class IoEnzoWriter;
 class CkIndex_IoEnzoWriter;
 class CProxy_IoEnzoWriter;
 class CProxyElement_IoEnzoWriter;
 class CProxySection_IoEnzoWriter;
/* --------------- index object ------------------ */
class CkIndex_IoEnzoWriter:public CkIndex_IoWriter{
  public:
    typedef IoEnzoWriter local_t;
    typedef CkIndex_IoEnzoWriter index_t;
    typedef CProxy_IoEnzoWriter proxy_t;
    typedef CProxyElement_IoEnzoWriter element_t;
    typedef CProxySection_IoEnzoWriter section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: IoEnzoWriter();
     */
    // Entry point registration at startup
    
    static int reg_IoEnzoWriter_void();
    // Entry point index lookup
    
    inline static int idx_IoEnzoWriter_void() {
      static int epidx = reg_IoEnzoWriter_void();
      return epidx;
    }

    
    static int ckNew() { return idx_IoEnzoWriter_void(); }
    
    static void _call_IoEnzoWriter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoEnzoWriter_void(void* impl_msg, void* impl_obj);
    /* DECLS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
     */
    // Entry point registration at startup
    
    static int reg_IoEnzoWriter_marshall2();
    // Entry point index lookup
    
    inline static int idx_IoEnzoWriter_marshall2() {
      static int epidx = reg_IoEnzoWriter_marshall2();
      return epidx;
    }

    
    static int ckNew(int num_files, const std::string &ordering, int monitor_iter) { return idx_IoEnzoWriter_marshall2(); }
    
    static void _call_IoEnzoWriter_marshall2(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoEnzoWriter_marshall2(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_IoEnzoWriter_marshall2(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_IoEnzoWriter_marshall2(PUP::er &p,void *msg);
    /* DECLS: void p_write(EnzoMsgCheck* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_p_write_EnzoMsgCheck();
    // Entry point index lookup
    
    inline static int idx_p_write_EnzoMsgCheck() {
      static int epidx = reg_p_write_EnzoMsgCheck();
      return epidx;
    }

    
    inline static int idx_p_write(void (IoEnzoWriter::*)(EnzoMsgCheck* impl_msg) ) {
      return idx_p_write_EnzoMsgCheck();
    }


    
    static int p_write(EnzoMsgCheck* impl_msg) { return idx_p_write_EnzoMsgCheck(); }
    
    static void _call_p_write_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_write_EnzoMsgCheck(void* impl_msg, void* impl_obj);
    /* DECLS: IoEnzoWriter(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_IoEnzoWriter_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_IoEnzoWriter_CkMigrateMessage() {
      static int epidx = reg_IoEnzoWriter_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_IoEnzoWriter_CkMigrateMessage(); }
    
    static void _call_IoEnzoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_IoEnzoWriter_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_IoEnzoWriter : public CProxyElement_IoWriter{
  public:
    typedef IoEnzoWriter local_t;
    typedef CkIndex_IoEnzoWriter index_t;
    typedef CProxy_IoEnzoWriter proxy_t;
    typedef CProxyElement_IoEnzoWriter element_t;
    typedef CProxySection_IoEnzoWriter section_t;

    using array_index_t = CkArrayIndex1D;

    /* TRAM aggregators */

    CProxyElement_IoEnzoWriter(void) {
    }
    CProxyElement_IoEnzoWriter(const ArrayElement *e) : CProxyElement_IoWriter(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxyElement_IoWriter::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxyElement_IoWriter::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxyElement_IoWriter::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxyElement_IoWriter::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxyElement_IoWriter::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxyElement_IoWriter::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxyElement_IoWriter::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxyElement_IoWriter::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxyElement_IoWriter::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxyElement_IoWriter::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxyElement_IoWriter::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxyElement_IoWriter::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxyElement_IoWriter::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxyElement_IoWriter::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxyElement_IoWriter::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxyElement_IoWriter::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxyElement_IoWriter::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IoWriter::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxyElement_IoWriter::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxyElement_IoWriter::ckSetReductionClient(cb); }

    inline void ckInsert(CkArrayMessage *m,int ctor,int onPe)
    { CProxyElement_IoWriter::ckInsert(m,ctor,onPe); }
    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0) const
    { CProxyElement_IoWriter::ckSend(m,ep,opts); }
    inline void *ckSendSync(CkArrayMessage *m, int ep) const
    { return CProxyElement_IoWriter::ckSendSync(m,ep); }
    inline const CkArrayIndex &ckGetIndex() const
    { return CProxyElement_IoWriter::ckGetIndex(); }

    IoEnzoWriter *ckLocal(void) const
    { return (IoEnzoWriter *)CProxyElement_IoWriter::ckLocal(); }

    CProxyElement_IoEnzoWriter(const CkArrayID &aid,const CkArrayIndex1D &idx,CK_DELCTOR_PARAM)
        :CProxyElement_IoWriter(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoEnzoWriter(const CkArrayID &aid,const CkArrayIndex1D &idx)
        :CProxyElement_IoWriter(aid,idx)
    {
}

    CProxyElement_IoEnzoWriter(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_IoWriter(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_IoEnzoWriter(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_IoWriter(aid,idx)
    {
}
/* DECLS: IoEnzoWriter();
 */
    
    void insert(int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
    
    void insert(int num_files, const std::string &ordering, int monitor_iter, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_write(EnzoMsgCheck* impl_msg);
 */
    
    void p_write(EnzoMsgCheck* impl_msg) ;

/* DECLS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_IoEnzoWriter : public CProxy_IoWriter{
  public:
    typedef IoEnzoWriter local_t;
    typedef CkIndex_IoEnzoWriter index_t;
    typedef CProxy_IoEnzoWriter proxy_t;
    typedef CProxyElement_IoEnzoWriter element_t;
    typedef CProxySection_IoEnzoWriter section_t;

    using array_index_t = CkArrayIndex1D;
    CProxy_IoEnzoWriter(void) {
    }
    CProxy_IoEnzoWriter(const ArrayElement *e) : CProxy_IoWriter(e){
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_IoWriter::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_IoWriter::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_IoWriter::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxy_IoWriter::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_IoWriter::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_IoWriter::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_IoWriter::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_IoWriter::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxy_IoWriter::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxy_IoWriter::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxy_IoWriter::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxy_IoWriter::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxy_IoWriter::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxy_IoWriter::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxy_IoWriter::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxy_IoWriter::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxy_IoWriter::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IoWriter::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxy_IoWriter::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxy_IoWriter::ckSetReductionClient(cb); }

    // Generalized array indexing:
    CProxyElement_IoEnzoWriter operator [] (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoEnzoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoEnzoWriter operator() (const CkArrayIndex1D &idx) const
    { return CProxyElement_IoEnzoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_IoEnzoWriter operator [] (int idx) const 
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoWriter operator () (int idx) const 
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), CkArrayIndex1D(idx), CK_DELCTOR_CALL);}
    CProxy_IoEnzoWriter(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_IoWriter(aid,CK_DELCTOR_ARGS) {}
    CProxy_IoEnzoWriter(const CkArrayID &aid) 
        :CProxy_IoWriter(aid) {}
/* DECLS: IoEnzoWriter();
 */
    
    static CkArrayID ckNew(const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);
    static CkArrayID ckNew(const int s1, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
    
    static CkArrayID ckNew(int num_files, const std::string &ordering, int monitor_iter, const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(int num_files, const std::string &ordering, int monitor_iter, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);
    static CkArrayID ckNew(int num_files, const std::string &ordering, int monitor_iter, const int s1, const CkEntryOptions *impl_e_opts=NULL);
    static void ckNew(int num_files, const std::string &ordering, int monitor_iter, const int s1, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_write(EnzoMsgCheck* impl_msg);
 */
    
    void p_write(EnzoMsgCheck* impl_msg) ;

/* DECLS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_IoEnzoWriter : public CProxySection_IoWriter{
  public:
    typedef IoEnzoWriter local_t;
    typedef CkIndex_IoEnzoWriter index_t;
    typedef CProxy_IoEnzoWriter proxy_t;
    typedef CProxyElement_IoEnzoWriter element_t;
    typedef CProxySection_IoEnzoWriter section_t;

    using array_index_t = CkArrayIndex1D;
    CProxySection_IoEnzoWriter(void) {
    }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxySection_IoWriter::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxySection_IoWriter::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxySection_IoWriter::pup(p);
    }

    int ckIsDelegated(void) const
    { return CProxySection_IoWriter::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxySection_IoWriter::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxySection_IoWriter::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxySection_IoWriter::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxySection_IoWriter::ckCheck(); }
    inline operator CkArrayID () const
    { return ckGetArrayID(); }
    inline CkArrayID ckGetArrayID(void) const
    { return CProxySection_IoWriter::ckGetArrayID(); }
    inline CkArray *ckLocalBranch(void) const
    { return CProxySection_IoWriter::ckLocalBranch(); }
    inline CkLocMgr *ckLocMgr(void) const
    { return CProxySection_IoWriter::ckLocMgr(); }

    inline static CkArrayID ckCreateEmptyArray(CkArrayOptions opts = CkArrayOptions())
    { return CProxySection_IoWriter::ckCreateEmptyArray(opts); }
    inline static void ckCreateEmptyArrayAsync(CkCallback cb, CkArrayOptions opts = CkArrayOptions())
    { CProxySection_IoWriter::ckCreateEmptyArrayAsync(cb, opts); }
    inline static CkArrayID ckCreateArray(CkArrayMessage *m,int ctor,const CkArrayOptions &opts)
    { return CProxySection_IoWriter::ckCreateArray(m,ctor,opts); }
    inline void ckInsertIdx(CkArrayMessage *m,int ctor,int onPe,const CkArrayIndex &idx)
    { CProxySection_IoWriter::ckInsertIdx(m,ctor,onPe,idx); }
    inline void doneInserting(void)
    { CProxySection_IoWriter::doneInserting(); }

    inline void ckBroadcast(CkArrayMessage *m, int ep, int opts=0) const
    { CProxySection_IoWriter::ckBroadcast(m,ep,opts); }
    inline void setReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IoWriter::setReductionClient(fn,param); }
    inline void ckSetReductionClient(CkReductionClientFn fn,void *param=NULL) const
    { CProxySection_IoWriter::ckSetReductionClient(fn,param); }
    inline void ckSetReductionClient(CkCallback *cb) const
    { CProxySection_IoWriter::ckSetReductionClient(cb); }

    inline void ckSend(CkArrayMessage *m, int ep, int opts = 0)
    { CProxySection_IoWriter::ckSend(m,ep,opts); }
    inline CkSectionInfo &ckGetSectionInfo()
    { return CProxySection_IoWriter::ckGetSectionInfo(); }
    inline CkSectionID *ckGetSectionIDs()
    { return CProxySection_IoWriter::ckGetSectionIDs(); }
    inline CkSectionID &ckGetSectionID()
    { return CProxySection_IoWriter::ckGetSectionID(); }
    inline CkSectionID &ckGetSectionID(int i)
    { return CProxySection_IoWriter::ckGetSectionID(i); }
    inline CkArrayID ckGetArrayIDn(int i) const
    { return CProxySection_IoWriter::ckGetArrayIDn(i); } 
    inline CkArrayIndex *ckGetArrayElements() const
    { return CProxySection_IoWriter::ckGetArrayElements(); }
    inline CkArrayIndex *ckGetArrayElements(int i) const
    { return CProxySection_IoWriter::ckGetArrayElements(i); }
    inline int ckGetNumElements() const
    { return CProxySection_IoWriter::ckGetNumElements(); } 
    inline int ckGetNumElements(int i) const
    { return CProxySection_IoWriter::ckGetNumElements(i); }    // Generalized array indexing:
    CProxyElement_IoEnzoWriter operator [] (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoWriter operator() (const CkArrayIndex1D &idx) const
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoWriter operator [] (int idx) const 
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
    CProxyElement_IoEnzoWriter operator () (int idx) const 
        {return CProxyElement_IoEnzoWriter(ckGetArrayID(), *(CkArrayIndex1D*)&ckGetArrayElements()[idx], CK_DELCTOR_CALL);}
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
    CProxySection_IoEnzoWriter(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_IoWriter(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoWriter(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_IoWriter(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoWriter(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_IoWriter(aid,elems,nElems, factor) {}
    CProxySection_IoEnzoWriter(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_IoWriter(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_IoEnzoWriter(const CkSectionID &sid)  
        :CProxySection_IoWriter(sid) { ckAutoDelegate(); }
    CProxySection_IoEnzoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_IoWriter(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_IoWriter(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_IoEnzoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_IoWriter(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_IoEnzoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_IoWriter(aid,elems) { ckAutoDelegate(); }
    CProxySection_IoEnzoWriter(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_IoWriter(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_IoEnzoWriter(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
        :CProxySection_IoWriter(aid,elems, factor) { ckAutoDelegate(); }
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
      CProxySection_IoWriter::ckAutoDelegate(opts);
    } 
    void setReductionClient(CkCallback *cb) {
      CProxySection_IoWriter::setReductionClient(cb);
    } 
    void resetSection() {
      CProxySection_IoWriter::resetSection();
    } 
    static void contribute(CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, int userData=-1, int fragSize=-1);
    static void contribute(CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    static void contribute(int dataSize,void *data,CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
    template <typename T>
    static void contribute(std::vector<T> &data, CkReduction::reducerType type, CkSectionInfo &sid, const CkCallback &cb, int userData=-1, int fragSize=-1);
/* DECLS: IoEnzoWriter();
 */
    

/* DECLS: IoEnzoWriter(int num_files, const std::string &ordering, int monitor_iter);
 */
    

/* DECLS: void p_write(EnzoMsgCheck* impl_msg);
 */
    
    void p_write(EnzoMsgCheck* impl_msg) ;

/* DECLS: IoEnzoWriter(CkMigrateMessage* impl_msg);
 */

};
#define IoEnzoWriter_SDAG_CODE 
typedef CBaseT1<IoWriter, CProxy_IoEnzoWriter>CBase_IoEnzoWriter;

/* DECLS: array EnzoLevelArray: ArrayElement{
EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
void p_request_data();
void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
void p_done(const Index &impl_noname_4);
EnzoLevelArray(CkMigrateMessage* impl_msg);
};
 */
 class EnzoLevelArray;
 class CkIndex_EnzoLevelArray;
 class CProxy_EnzoLevelArray;
 class CProxyElement_EnzoLevelArray;
 class CProxySection_EnzoLevelArray;
/* --------------- index object ------------------ */
class CkIndex_EnzoLevelArray:public CkIndex_ArrayElement{
  public:
    typedef EnzoLevelArray local_t;
    typedef CkIndex_EnzoLevelArray index_t;
    typedef CProxy_EnzoLevelArray proxy_t;
    typedef CProxyElement_EnzoLevelArray element_t;
    typedef CProxySection_EnzoLevelArray section_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
     */
    // Entry point registration at startup
    
    static int reg_EnzoLevelArray_marshall1();
    // Entry point index lookup
    
    inline static int idx_EnzoLevelArray_marshall1() {
      static int epidx = reg_EnzoLevelArray_marshall1();
      return epidx;
    }

    
    static int ckNew(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz) { return idx_EnzoLevelArray_marshall1(); }
    
    static void _call_EnzoLevelArray_marshall1(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoLevelArray_marshall1(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_EnzoLevelArray_marshall1(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_EnzoLevelArray_marshall1(PUP::er &p,void *msg);
    /* DECLS: void p_request_data();
     */
    // Entry point registration at startup
    
    static int reg_p_request_data_void();
    // Entry point index lookup
    
    inline static int idx_p_request_data_void() {
      static int epidx = reg_p_request_data_void();
      return epidx;
    }

    
    inline static int idx_p_request_data(void (EnzoLevelArray::*)() ) {
      return idx_p_request_data_void();
    }


    
    static int p_request_data() { return idx_p_request_data_void(); }
    
    static void _call_p_request_data_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_request_data_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
     */
    // Entry point registration at startup
    
    static int reg_p_transfer_data_marshall3();
    // Entry point index lookup
    
    inline static int idx_p_transfer_data_marshall3() {
      static int epidx = reg_p_transfer_data_marshall3();
      return epidx;
    }

    
    inline static int idx_p_transfer_data(void (EnzoLevelArray::*)(const Index &impl_noname_3, int nf, const enzo_float *field_data) ) {
      return idx_p_transfer_data_marshall3();
    }


    
    static int p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data) { return idx_p_transfer_data_marshall3(); }
    
    static void _call_p_transfer_data_marshall3(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_transfer_data_marshall3(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_transfer_data_marshall3(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_transfer_data_marshall3(PUP::er &p,void *msg);
    /* DECLS: void p_done(const Index &impl_noname_4);
     */
    // Entry point registration at startup
    
    static int reg_p_done_marshall4();
    // Entry point index lookup
    
    inline static int idx_p_done_marshall4() {
      static int epidx = reg_p_done_marshall4();
      return epidx;
    }

    
    inline static int idx_p_done(void (EnzoLevelArray::*)(const Index &impl_noname_4) ) {
      return idx_p_done_marshall4();
    }


    
    static int p_done(const Index &impl_noname_4) { return idx_p_done_marshall4(); }
    
    static void _call_p_done_marshall4(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_done_marshall4(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_done_marshall4(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_done_marshall4(PUP::er &p,void *msg);
    /* DECLS: EnzoLevelArray(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_EnzoLevelArray_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_EnzoLevelArray_CkMigrateMessage() {
      static int epidx = reg_EnzoLevelArray_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_EnzoLevelArray_CkMigrateMessage(); }
    
    static void _call_EnzoLevelArray_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_EnzoLevelArray_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
 class CProxyElement_EnzoLevelArray : public CProxyElement_ArrayElement{
  public:
    typedef EnzoLevelArray local_t;
    typedef CkIndex_EnzoLevelArray index_t;
    typedef CProxy_EnzoLevelArray proxy_t;
    typedef CProxyElement_EnzoLevelArray element_t;
    typedef CProxySection_EnzoLevelArray section_t;

    using array_index_t = CkArrayIndexIndex3;

    /* TRAM aggregators */

    CProxyElement_EnzoLevelArray(void) {
    }
    CProxyElement_EnzoLevelArray(const ArrayElement *e) : CProxyElement_ArrayElement(e){
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

    EnzoLevelArray *ckLocal(void) const
    { return (EnzoLevelArray *)CProxyElement_ArrayElement::ckLocal(); }

    CProxyElement_EnzoLevelArray(const CkArrayID &aid,const CkArrayIndexIndex3 &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_EnzoLevelArray(const CkArrayID &aid,const CkArrayIndexIndex3 &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}

    CProxyElement_EnzoLevelArray(const CkArrayID &aid,const CkArrayIndex &idx,CK_DELCTOR_PARAM)
        :CProxyElement_ArrayElement(aid,idx,CK_DELCTOR_ARGS)
    {
}
    CProxyElement_EnzoLevelArray(const CkArrayID &aid,const CkArrayIndex &idx)
        :CProxyElement_ArrayElement(aid,idx)
    {
}
/* DECLS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
    
    void insert(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, int onPE=-1, const CkEntryOptions *impl_e_opts=NULL);
/* DECLS: void p_request_data();
 */
    
    void p_request_data(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
    
    void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_done(const Index &impl_noname_4);
 */
    
    void p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- collective proxy -------------- */
 class CProxy_EnzoLevelArray : public CProxy_ArrayElement{
  public:
    typedef EnzoLevelArray local_t;
    typedef CkIndex_EnzoLevelArray index_t;
    typedef CProxy_EnzoLevelArray proxy_t;
    typedef CProxyElement_EnzoLevelArray element_t;
    typedef CProxySection_EnzoLevelArray section_t;

    using array_index_t = CkArrayIndexIndex3;
    CProxy_EnzoLevelArray(void) {
    }
    CProxy_EnzoLevelArray(const ArrayElement *e) : CProxy_ArrayElement(e){
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

    // Empty array construction
    static CkArrayID ckNew(CkArrayOptions opts = CkArrayOptions()) { return ckCreateEmptyArray(opts); }
    static void      ckNew(CkCallback cb, CkArrayOptions opts = CkArrayOptions()) { ckCreateEmptyArrayAsync(cb, opts); }

    // Generalized array indexing:
    CProxyElement_EnzoLevelArray operator [] (const CkArrayIndexIndex3 &idx) const
    { return CProxyElement_EnzoLevelArray(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxyElement_EnzoLevelArray operator() (const CkArrayIndexIndex3 &idx) const
    { return CProxyElement_EnzoLevelArray(ckGetArrayID(), idx, CK_DELCTOR_CALL); }
    CProxy_EnzoLevelArray(const CkArrayID &aid,CK_DELCTOR_PARAM) 
        :CProxy_ArrayElement(aid,CK_DELCTOR_ARGS) {}
    CProxy_EnzoLevelArray(const CkArrayID &aid) 
        :CProxy_ArrayElement(aid) {}
/* DECLS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
    
    static CkArrayID ckNew(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, const CkArrayOptions &opts = CkArrayOptions(), const CkEntryOptions *impl_e_opts=NULL);
    static void      ckNew(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz, const CkArrayOptions &opts, CkCallback _ck_array_creation_cb, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_request_data();
 */
    
    void p_request_data(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
    
    void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_done(const Index &impl_noname_4);
 */
    
    void p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */

};
/* ---------------- section proxy -------------- */
 class CProxySection_EnzoLevelArray : public CProxySection_ArrayElement{
  public:
    typedef EnzoLevelArray local_t;
    typedef CkIndex_EnzoLevelArray index_t;
    typedef CProxy_EnzoLevelArray proxy_t;
    typedef CProxyElement_EnzoLevelArray element_t;
    typedef CProxySection_EnzoLevelArray section_t;

    using array_index_t = CkArrayIndexIndex3;
    CProxySection_EnzoLevelArray(void) {
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
    CProxyElement_EnzoLevelArray operator [] (const CkArrayIndexIndex3 &idx) const
        {return CProxyElement_EnzoLevelArray(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxyElement_EnzoLevelArray operator() (const CkArrayIndexIndex3 &idx) const
        {return CProxyElement_EnzoLevelArray(ckGetArrayID(), idx, CK_DELCTOR_CALL);}
    CProxySection_EnzoLevelArray(const CkArrayID &aid, CkArrayIndex *elems, int nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoLevelArray(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoLevelArray(const CkArrayID &aid, CkArrayIndex *elems, int nElems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems,nElems, factor) {}
    CProxySection_EnzoLevelArray(const CkArrayID &aid, const std::vector<CkArrayIndex> &elems, int factor=USE_DEFAULT_BRANCH_FACTOR) 
        :CProxySection_ArrayElement(aid,elems, factor) { ckAutoDelegate(); }
    CProxySection_EnzoLevelArray(const CkSectionID &sid)  
        :CProxySection_ArrayElement(sid) { ckAutoDelegate(); }
    CProxySection_EnzoLevelArray(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(n,aid,elems,nElems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoLevelArray(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, CK_DELCTOR_PARAM) 
        :CProxySection_ArrayElement(aid,elems,CK_DELCTOR_ARGS) {}
    CProxySection_EnzoLevelArray(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems) 
        :CProxySection_ArrayElement(n,aid,elems,nElems) { ckAutoDelegate(); }
    CProxySection_EnzoLevelArray(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems) 
        :CProxySection_ArrayElement(aid,elems) { ckAutoDelegate(); }
    CProxySection_EnzoLevelArray(int n, const CkArrayID *aid, CkArrayIndex const * const *elems, const int *nElems, int factor) 
        :CProxySection_ArrayElement(n,aid,elems,nElems, factor) { ckAutoDelegate(); }
    CProxySection_EnzoLevelArray(const std::vector<CkArrayID> &aid, const std::vector<std::vector<CkArrayIndex> > &elems, int factor) 
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
/* DECLS: EnzoLevelArray(const std::string &impl_noname_2, int level_base, int level_array, int level_infer, int nax, int nay, int naz);
 */
    

/* DECLS: void p_request_data();
 */
    
    void p_request_data(const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data);
 */
    
    void p_transfer_data(const Index &impl_noname_3, int nf, const enzo_float *field_data, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: void p_done(const Index &impl_noname_4);
 */
    
    void p_done(const Index &impl_noname_4, const CkEntryOptions *impl_e_opts=NULL) ;

/* DECLS: EnzoLevelArray(CkMigrateMessage* impl_msg);
 */

};
#define EnzoLevelArray_SDAG_CODE 
typedef CBaseT1<ArrayElementT<Index3>, CProxy_EnzoLevelArray>CBase_EnzoLevelArray;





















































































































/* ---------------- method closures -------------- */
class Closure_EnzoSimulation {
  public:



    struct p_get_msg_refine_3_closure;


    struct p_get_msg_check_4_closure;



    struct p_method_balance_check_6_closure;



    struct p_check_done_8_closure;


    struct p_set_io_writer_9_closure;


    struct p_infer_set_array_count_10_closure;


    struct p_infer_array_created_11_closure;


    struct p_infer_done_12_closure;


    struct p_set_io_reader_13_closure;


    struct p_io_reader_created_14_closure;


    struct p_restart_next_level_15_closure;


    struct p_restart_level_created_16_closure;


    struct p_set_level_array_17_closure;


    struct p_fbnet_concatenate_sphere_lists_18_closure;


    struct p_fbnet_done_19_closure;


};

/* ---------------- method closures -------------- */
class Closure_EnzoBlock {
  public:





    struct p_method_feedback_starss_end_5_closure;


    struct p_method_m1_closure_solve_transport_eqn_6_closure;





    struct p_method_balance_migrate_10_closure;


    struct p_method_balance_done_11_closure;


    struct p_method_gravity_continue_12_closure;


    struct p_method_gravity_end_13_closure;


    struct p_method_infer_merge_masks_14_closure;


    struct p_method_infer_count_arrays_15_closure;


    struct p_method_infer_request_data_16_closure;


    struct p_method_infer_update_17_closure;


    struct p_method_infer_exit_18_closure;



    struct p_method_fbnet_exit_20_closure;


    struct p_check_write_first_21_closure;


    struct p_check_write_next_22_closure;


    struct p_check_done_23_closure;


    struct p_restart_refine_24_closure;



    struct p_restart_done_26_closure;


    struct p_method_accretion_end_27_closure;


    struct p_solver_cg_matvec_28_closure;





    struct p_solver_cg_loop_2_32_closure;










    struct p_solver_bicgstab_loop_2_41_closure;


    struct p_solver_bicgstab_loop_3_42_closure;


    struct p_solver_bicgstab_loop_8_43_closure;


    struct p_solver_bicgstab_loop_9_44_closure;


    struct p_dot_recv_parent_45_closure;


    struct p_dot_recv_children_46_closure;




    struct p_solver_dd_solve_coarse_49_closure;


    struct p_solver_dd_solve_domain_50_closure;


    struct p_solver_dd_last_smooth_51_closure;




    struct p_solver_jacobi_continue_54_closure;


    struct p_solver_mg0_restrict_55_closure;


    struct p_solver_mg0_solve_coarse_56_closure;


    struct p_solver_mg0_post_smooth_57_closure;


    struct p_solver_mg0_last_smooth_58_closure;






};

/* ---------------- method closures -------------- */
class Closure_IoEnzoReader {
  public:


    struct p_init_root_2_closure;


    struct p_create_level_3_closure;


    struct p_init_level_4_closure;


    struct p_block_created_5_closure;


    struct p_block_ready_6_closure;


};

/* ---------------- method closures -------------- */
class Closure_IoEnzoWriter {
  public:




};

/* ---------------- method closures -------------- */
class Closure_EnzoLevelArray {
  public:


    struct p_request_data_2_closure;


    struct p_transfer_data_3_closure;


    struct p_done_4_closure;


};

extern void _registerenzo(void);
#endif
