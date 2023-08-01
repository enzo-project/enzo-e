#ifndef _DECL_main_mesh_H_
#define _DECL_main_mesh_H_
#include "charm++.h"
#include "envelope.h"
#include <memory>
#include "sdag.h"
#include "mesh.decl.h"

/* DECLS: readonly CProxy_Main proxy_main;
 */

/* DECLS: mainchare Main: Chare{
Main(CkArgMsg* impl_msg);
void p_exit(int count_blocks);
void p_checkpoint_output(int count, const std::string &dir);
void p_initial_exit();
void p_adapt_enter();
void p_adapt_called();
void p_adapt_end();
void p_adapt_update();
void p_adapt_exit();
void p_compute_enter();
void p_compute_continue();
void p_compute_exit();
void p_output_enter();
void p_output_exit();
void p_stopping_enter();
void p_stopping_balance();
void p_stopping_exit();
void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
void p_exit();
Main(CkMigrateMessage* impl_msg);
};
 */
 class Main;
 class CkIndex_Main;
 class CProxy_Main;
/* --------------- index object ------------------ */
class CkIndex_Main:public CkIndex_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    static int __idx;
    static void __register(const char *s, size_t size);
    /* DECLS: Main(CkArgMsg* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_Main_CkArgMsg();
    // Entry point index lookup
    
    inline static int idx_Main_CkArgMsg() {
      static int epidx = reg_Main_CkArgMsg();
      return epidx;
    }

    
    static int ckNew(CkArgMsg* impl_msg) { return idx_Main_CkArgMsg(); }
    
    static void _call_Main_CkArgMsg(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Main_CkArgMsg(void* impl_msg, void* impl_obj);
    /* DECLS: void p_exit(int count_blocks);
     */
    // Entry point registration at startup
    
    static int reg_p_exit_marshall2();
    // Entry point index lookup
    
    inline static int idx_p_exit_marshall2() {
      static int epidx = reg_p_exit_marshall2();
      return epidx;
    }

    
    inline static int idx_p_exit(void (Main::*)(int count_blocks) ) {
      return idx_p_exit_marshall2();
    }


    
    static int p_exit(int count_blocks) { return idx_p_exit_marshall2(); }
    
    static void _call_p_exit_marshall2(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_exit_marshall2(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_exit_marshall2(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_exit_marshall2(PUP::er &p,void *msg);
    /* DECLS: void p_checkpoint_output(int count, const std::string &dir);
     */
    // Entry point registration at startup
    
    static int reg_p_checkpoint_output_marshall3();
    // Entry point index lookup
    
    inline static int idx_p_checkpoint_output_marshall3() {
      static int epidx = reg_p_checkpoint_output_marshall3();
      return epidx;
    }

    
    inline static int idx_p_checkpoint_output(void (Main::*)(int count, const std::string &dir) ) {
      return idx_p_checkpoint_output_marshall3();
    }


    
    static int p_checkpoint_output(int count, const std::string &dir) { return idx_p_checkpoint_output_marshall3(); }
    
    static void _call_p_checkpoint_output_marshall3(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_checkpoint_output_marshall3(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_checkpoint_output_marshall3(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_checkpoint_output_marshall3(PUP::er &p,void *msg);
    /* DECLS: void p_initial_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_initial_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_initial_exit_void() {
      static int epidx = reg_p_initial_exit_void();
      return epidx;
    }

    
    inline static int idx_p_initial_exit(void (Main::*)() ) {
      return idx_p_initial_exit_void();
    }


    
    static int p_initial_exit() { return idx_p_initial_exit_void(); }
    
    static void _call_p_initial_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_initial_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_enter_void() {
      static int epidx = reg_p_adapt_enter_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_enter(void (Main::*)() ) {
      return idx_p_adapt_enter_void();
    }


    
    static int p_adapt_enter() { return idx_p_adapt_enter_void(); }
    
    static void _call_p_adapt_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_called();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_called_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_called_void() {
      static int epidx = reg_p_adapt_called_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_called(void (Main::*)() ) {
      return idx_p_adapt_called_void();
    }


    
    static int p_adapt_called() { return idx_p_adapt_called_void(); }
    
    static void _call_p_adapt_called_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_called_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_end();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_end_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_end_void() {
      static int epidx = reg_p_adapt_end_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_end(void (Main::*)() ) {
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

    
    inline static int idx_p_adapt_update(void (Main::*)() ) {
      return idx_p_adapt_update_void();
    }


    
    static int p_adapt_update() { return idx_p_adapt_update_void(); }
    
    static void _call_p_adapt_update_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_update_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_adapt_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_adapt_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_adapt_exit_void() {
      static int epidx = reg_p_adapt_exit_void();
      return epidx;
    }

    
    inline static int idx_p_adapt_exit(void (Main::*)() ) {
      return idx_p_adapt_exit_void();
    }


    
    static int p_adapt_exit() { return idx_p_adapt_exit_void(); }
    
    static void _call_p_adapt_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_adapt_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_compute_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_compute_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_compute_enter_void() {
      static int epidx = reg_p_compute_enter_void();
      return epidx;
    }

    
    inline static int idx_p_compute_enter(void (Main::*)() ) {
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

    
    inline static int idx_p_compute_continue(void (Main::*)() ) {
      return idx_p_compute_continue_void();
    }


    
    static int p_compute_continue() { return idx_p_compute_continue_void(); }
    
    static void _call_p_compute_continue_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_compute_continue_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_compute_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_compute_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_compute_exit_void() {
      static int epidx = reg_p_compute_exit_void();
      return epidx;
    }

    
    inline static int idx_p_compute_exit(void (Main::*)() ) {
      return idx_p_compute_exit_void();
    }


    
    static int p_compute_exit() { return idx_p_compute_exit_void(); }
    
    static void _call_p_compute_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_compute_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_output_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_output_enter_void() {
      static int epidx = reg_p_output_enter_void();
      return epidx;
    }

    
    inline static int idx_p_output_enter(void (Main::*)() ) {
      return idx_p_output_enter_void();
    }


    
    static int p_output_enter() { return idx_p_output_enter_void(); }
    
    static void _call_p_output_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_output_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_output_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_output_exit_void() {
      static int epidx = reg_p_output_exit_void();
      return epidx;
    }

    
    inline static int idx_p_output_exit(void (Main::*)() ) {
      return idx_p_output_exit_void();
    }


    
    static int p_output_exit() { return idx_p_output_exit_void(); }
    
    static void _call_p_output_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_output_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_enter();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_enter_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_enter_void() {
      static int epidx = reg_p_stopping_enter_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_enter(void (Main::*)() ) {
      return idx_p_stopping_enter_void();
    }


    
    static int p_stopping_enter() { return idx_p_stopping_enter_void(); }
    
    static void _call_p_stopping_enter_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_enter_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_balance();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_balance_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_balance_void() {
      static int epidx = reg_p_stopping_balance_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_balance(void (Main::*)() ) {
      return idx_p_stopping_balance_void();
    }


    
    static int p_stopping_balance() { return idx_p_stopping_balance_void(); }
    
    static void _call_p_stopping_balance_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_balance_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_stopping_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_stopping_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_stopping_exit_void() {
      static int epidx = reg_p_stopping_exit_void();
      return epidx;
    }

    
    inline static int idx_p_stopping_exit(void (Main::*)() ) {
      return idx_p_stopping_exit_void();
    }


    
    static int p_stopping_exit() { return idx_p_stopping_exit_void(); }
    
    static void _call_p_stopping_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_stopping_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
     */
    // Entry point registration at startup
    
    static int reg_p_text_file_write_marshall18();
    // Entry point index lookup
    
    inline static int idx_p_text_file_write_marshall18() {
      static int epidx = reg_p_text_file_write_marshall18();
      return epidx;
    }

    
    inline static int idx_p_text_file_write(void (Main::*)(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count) ) {
      return idx_p_text_file_write_marshall18();
    }


    
    static int p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count) { return idx_p_text_file_write_marshall18(); }
    
    static void _call_p_text_file_write_marshall18(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_text_file_write_marshall18(void* impl_msg, void* impl_obj);
    
    static int _callmarshall_p_text_file_write_marshall18(char* impl_buf, void* impl_obj_void);
    
    static void _marshallmessagepup_p_text_file_write_marshall18(PUP::er &p,void *msg);
    /* DECLS: void p_exit();
     */
    // Entry point registration at startup
    
    static int reg_p_exit_void();
    // Entry point index lookup
    
    inline static int idx_p_exit_void() {
      static int epidx = reg_p_exit_void();
      return epidx;
    }

    
    inline static int idx_p_exit(void (Main::*)() ) {
      return idx_p_exit_void();
    }


    
    static int p_exit() { return idx_p_exit_void(); }
    
    static void _call_p_exit_void(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_p_exit_void(void* impl_msg, void* impl_obj);
    /* DECLS: Main(CkMigrateMessage* impl_msg);
     */
    // Entry point registration at startup
    
    static int reg_Main_CkMigrateMessage();
    // Entry point index lookup
    
    inline static int idx_Main_CkMigrateMessage() {
      static int epidx = reg_Main_CkMigrateMessage();
      return epidx;
    }

    
    static int ckNew(CkMigrateMessage* impl_msg) { return idx_Main_CkMigrateMessage(); }
    
    static void _call_Main_CkMigrateMessage(void* impl_msg, void* impl_obj);
    
    static void _call_sdag_Main_CkMigrateMessage(void* impl_msg, void* impl_obj);
};
/* --------------- element proxy ------------------ */
class CProxy_Main:public CProxy_Chare{
  public:
    typedef Main local_t;
    typedef CkIndex_Main index_t;
    typedef CProxy_Main proxy_t;
    typedef CProxy_Main element_t;

    CProxy_Main(void) {};
    CProxy_Main(CkChareID __cid) : CProxy_Chare(__cid){  }
    CProxy_Main(const Chare *c) : CProxy_Chare(c){  }

    int ckIsDelegated(void) const
    { return CProxy_Chare::ckIsDelegated(); }
    inline CkDelegateMgr *ckDelegatedTo(void) const
    { return CProxy_Chare::ckDelegatedTo(); }
    inline CkDelegateData *ckDelegatedPtr(void) const
    { return CProxy_Chare::ckDelegatedPtr(); }
    CkGroupID ckDelegatedIdx(void) const
    { return CProxy_Chare::ckDelegatedIdx(); }

    inline void ckCheck(void) const
    { CProxy_Chare::ckCheck(); }
    const CkChareID &ckGetChareID(void) const
    { return CProxy_Chare::ckGetChareID(); }
    operator const CkChareID &(void) const
    { return ckGetChareID(); }

    void ckDelegate(CkDelegateMgr *dTo,CkDelegateData *dPtr=NULL)
    {       CProxy_Chare::ckDelegate(dTo,dPtr); }
    void ckUndelegate(void)
    {       CProxy_Chare::ckUndelegate(); }
    void pup(PUP::er &p)
    {       CProxy_Chare::pup(p);
    }

    void ckSetChareID(const CkChareID &c)
    {      CProxy_Chare::ckSetChareID(c); }
    Main *ckLocal(void) const
    { return (Main *)CkLocalChare(&ckGetChareID()); }
/* DECLS: Main(CkArgMsg* impl_msg);
 */
    static CkChareID ckNew(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);
    static void ckNew(CkArgMsg* impl_msg, CkChareID* pcid, int onPE=CK_PE_ANY);
    CProxy_Main(CkArgMsg* impl_msg, int onPE=CK_PE_ANY);

/* DECLS: void p_exit(int count_blocks);
 */
    
    void p_exit(int count_blocks, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_checkpoint_output(int count, const std::string &dir);
 */
    
    void p_checkpoint_output(int count, const std::string &dir, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_initial_exit();
 */
    
    void p_initial_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_adapt_enter();
 */
    
    void p_adapt_enter(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_adapt_called();
 */
    
    void p_adapt_called(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_adapt_end();
 */
    
    void p_adapt_end(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_adapt_update();
 */
    
    void p_adapt_update(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_adapt_exit();
 */
    
    void p_adapt_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_compute_enter();
 */
    
    void p_compute_enter(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_compute_continue();
 */
    
    void p_compute_continue(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_compute_exit();
 */
    
    void p_compute_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_output_enter();
 */
    
    void p_output_enter(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_output_exit();
 */
    
    void p_output_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_stopping_enter();
 */
    
    void p_stopping_enter(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_stopping_balance();
 */
    
    void p_stopping_balance(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_stopping_exit();
 */
    
    void p_stopping_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count);
 */
    
    void p_text_file_write(int nd, const char *dir, int nf, const char *file, int nl, const char *line, int count, const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: void p_exit();
 */
    
    void p_exit(const CkEntryOptions *impl_e_opts=NULL);

/* DECLS: Main(CkMigrateMessage* impl_msg);
 */

};
#define Main_SDAG_CODE 
typedef CBaseT1<Chare, CProxy_Main>CBase_Main;



/* ---------------- method closures -------------- */
class Closure_Main {
  public:


    struct p_exit_2_closure;


    struct p_checkpoint_output_3_closure;


    struct p_initial_exit_4_closure;


    struct p_adapt_enter_5_closure;


    struct p_adapt_called_6_closure;


    struct p_adapt_end_7_closure;


    struct p_adapt_update_8_closure;


    struct p_adapt_exit_9_closure;


    struct p_compute_enter_10_closure;


    struct p_compute_continue_11_closure;


    struct p_compute_exit_12_closure;


    struct p_output_enter_13_closure;


    struct p_output_exit_14_closure;


    struct p_stopping_enter_15_closure;


    struct p_stopping_balance_16_closure;


    struct p_stopping_exit_17_closure;


    struct p_text_file_write_18_closure;


    struct p_exit_19_closure;


};

extern void _registermain_mesh(void);
extern "C" void CkRegisterMainModule(void);
#endif
