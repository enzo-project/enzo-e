


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_get_msg_refine_2_closure : public SDAG::Closure {
            Index impl_noname_0;


      p_get_msg_refine_2_closure() {
        init();
      }
      p_get_msg_refine_2_closure(CkMigrateMessage*) {
        init();
      }
            Index & getP0() { return impl_noname_0;}
      void pup(PUP::er& __p) {
        __p | impl_noname_0;
        packClosure(__p);
      }
      virtual ~p_get_msg_refine_2_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_get_msg_refine_2_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::s_write_4_closure : public SDAG::Closure {
      

      s_write_4_closure() {
        init();
      }
      s_write_4_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~s_write_4_closure() {
      }
      PUPable_decl(SINGLE_ARG(s_write_4_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::r_write_checkpoint_output_6_closure : public SDAG::Closure {
      

      r_write_checkpoint_output_6_closure() {
        init();
      }
      r_write_checkpoint_output_6_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~r_write_checkpoint_output_6_closure() {
      }
      PUPable_decl(SINGLE_ARG(r_write_checkpoint_output_6_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_restart_enter_8_closure : public SDAG::Closure {
            std::string dir;


      p_restart_enter_8_closure() {
        init();
      }
      p_restart_enter_8_closure(CkMigrateMessage*) {
        init();
      }
            std::string & getP0() { return dir;}
      void pup(PUP::er& __p) {
        __p | dir;
        packClosure(__p);
      }
      virtual ~p_restart_enter_8_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_restart_enter_8_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_output_write_9_closure : public SDAG::Closure {
            int n;
            char *buffer;

      CkMarshallMsg* _impl_marshall;
      char* _impl_buf_in;
      int _impl_buf_size;

      p_output_write_9_closure() {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
      p_output_write_9_closure(CkMigrateMessage*) {
        init();
        _impl_marshall = 0;
        _impl_buf_in = 0;
        _impl_buf_size = 0;
      }
            int & getP0() { return n;}
            char *& getP1() { return buffer;}
      void pup(PUP::er& __p) {
        __p | n;
        packClosure(__p);
        __p | _impl_buf_size;
        bool hasMsg = (_impl_marshall != 0); __p | hasMsg;
        if (hasMsg) CkPupMessage(__p, (void**)&_impl_marshall);
        else PUParray(__p, _impl_buf_in, _impl_buf_size);
        if (__p.isUnpacking()) {
          char *impl_buf = _impl_marshall ? _impl_marshall->msgBuf : _impl_buf_in;
          PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
          impl_buf+=CK_ALIGN(implP.size(),16);
          buffer = (char *)(impl_buf+impl_off_buffer);
        }
      }
      virtual ~p_output_write_9_closure() {
        if (_impl_marshall) CmiFree(UsrToEnv(_impl_marshall));
      }
      PUPable_decl(SINGLE_ARG(p_output_write_9_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_output_start_11_closure : public SDAG::Closure {
            int index_output;


      p_output_start_11_closure() {
        init();
      }
      p_output_start_11_closure(CkMigrateMessage*) {
        init();
      }
            int & getP0() { return index_output;}
      void pup(PUP::er& __p) {
        __p | index_output;
        packClosure(__p);
      }
      virtual ~p_output_start_11_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_output_start_11_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_monitor_performance_13_closure : public SDAG::Closure {
      

      p_monitor_performance_13_closure() {
        init();
      }
      p_monitor_performance_13_closure(CkMigrateMessage*) {
        init();
      }
            void pup(PUP::er& __p) {
        packClosure(__p);
      }
      virtual ~p_monitor_performance_13_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_monitor_performance_13_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY

    struct Closure_Simulation::p_set_block_array_14_closure : public SDAG::Closure {
            CProxy_Block block_array;


      p_set_block_array_14_closure() {
        init();
      }
      p_set_block_array_14_closure(CkMigrateMessage*) {
        init();
      }
            CProxy_Block & getP0() { return block_array;}
      void pup(PUP::er& __p) {
        __p | block_array;
        packClosure(__p);
      }
      virtual ~p_set_block_array_14_closure() {
      }
      PUPable_decl(SINGLE_ARG(p_set_block_array_14_closure));
    };
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */


/* ---------------- method closures -------------- */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */



/* DEFS: readonly CProxy_Simulation proxy_simulation;
 */
extern CProxy_Simulation proxy_simulation;
#ifndef CK_TEMPLATES_ONLY
extern "C" void __xlater_roPup_proxy_simulation(void *_impl_pup_er) {
  PUP::er &_impl_p=*(PUP::er *)_impl_pup_er;
  _impl_p|proxy_simulation;
}
#endif /* CK_TEMPLATES_ONLY */


/* DEFS: group Simulation: IrrGroup{
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
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_Simulation::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(const char *filename, int n);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &impl_noname_0);
 */
void CProxyElement_Simulation::p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
void CProxyElement_Simulation::r_initialize_block_array(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void s_write();
 */
void CProxyElement_Simulation::s_write(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write(CkReductionMsg* impl_msg);
 */
void CProxyElement_Simulation::r_write(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write_checkpoint_output();
 */
void CProxyElement_Simulation::r_write_checkpoint_output(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_start(CkReductionMsg* impl_msg);
 */
void CProxyElement_Simulation::r_restart_start(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_enter(const std::string &dir);
 */
void CProxyElement_Simulation::p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int n, const char *buffer);
 */
void CProxyElement_Simulation::p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int n, const char *buffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
void CProxyElement_Simulation::r_output_barrier(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_start(int index_output);
 */
void CProxyElement_Simulation::p_output_start(int index_output, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int index_output
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
void CProxyElement_Simulation::r_monitor_performance_reduce(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_monitor_performance();
 */
void CProxyElement_Simulation::p_monitor_performance(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_block_array(const CProxy_Block &block_array);
 */
void CProxyElement_Simulation::p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_Block &block_array
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupPe(), ckGetGroupID());
  } else {
    CkSendMsgBranch(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupPe(), ckGetGroupID(),0);
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(const char *filename, int n);
 */
CkGroupID CProxy_Simulation::ckNew(const char *filename, int n, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const char *filename, int n
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_filename, impl_cnt_filename;
  impl_off_filename=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_filename=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_filename,filename,impl_cnt_filename);
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  CkGroupID gId = CkCreateGroup(CkIndex_Simulation::__idx, CkIndex_Simulation::idx_Simulation_marshall1(), impl_msg);
  return gId;
}
  CProxy_Simulation::CProxy_Simulation(const char *filename, int n, const CkEntryOptions *impl_e_opts)
{
  //Marshall: const char *filename, int n
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_filename, impl_cnt_filename;
  impl_off_filename=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_filename=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_off_filename;
    implP|impl_cnt_filename;
    implP|n;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_filename,filename,impl_cnt_filename);
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  ckSetGroupID(CkCreateGroup(CkIndex_Simulation::__idx, CkIndex_Simulation::idx_Simulation_marshall1(), impl_msg));
}

// Entry point registration function
int CkIndex_Simulation::reg_Simulation_marshall1() {
  int epidx = CkRegisterEp("Simulation(const char *filename, int n)",
      reinterpret_cast<CkCallFnPtr>(_call_Simulation_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_Simulation_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_Simulation_marshall1);

  return epidx;
}

void CkIndex_Simulation::_call_Simulation_marshall1(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  new (impl_obj_void) Simulation(filename, std::move(n.t));
}
int CkIndex_Simulation::_callmarshall_Simulation_marshall1(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  new (impl_obj_void) Simulation(filename, std::move(n.t));
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_Simulation_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const char *filename, int n*/
  PUP::fromMem implP(impl_buf);
  int impl_off_filename, impl_cnt_filename;
  implP|impl_off_filename;
  implP|impl_cnt_filename;
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *filename=(char *)(impl_buf+impl_off_filename);
  if (implDestP.hasComments()) implDestP.comment("filename");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*filename))<impl_cnt_filename;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|filename[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &impl_noname_0);
 */
void CProxy_Simulation::p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_get_msg_refine(const Index &impl_noname_0, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
  }
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_get_msg_refine(const Index &impl_noname_0, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const Index &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
  }
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_get_msg_refine_marshall2() {
  int epidx = CkRegisterEp("p_get_msg_refine(const Index &impl_noname_0)",
      reinterpret_cast<CkCallFnPtr>(_call_p_get_msg_refine_marshall2), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_get_msg_refine_marshall2);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_get_msg_refine_marshall2);

  return epidx;
}

void CkIndex_Simulation::_call_p_get_msg_refine_marshall2(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_0;
  implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_refine(std::move(impl_noname_0.t));
}
int CkIndex_Simulation::_callmarshall_p_get_msg_refine_marshall2(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const Index &impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_0;
  implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_get_msg_refine(std::move(impl_noname_0.t));
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_p_get_msg_refine_marshall2(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const Index &impl_noname_0*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<Index> impl_noname_0;
  implP|impl_noname_0;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_0");
  implDestP|impl_noname_0;
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_get_msg_refine_2_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
void CProxy_Simulation::r_initialize_block_array(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_initialize_block_array(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_initialize_block_array(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_initialize_block_array_CkReductionMsg() {
  int epidx = CkRegisterEp("r_initialize_block_array(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_initialize_block_array_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Simulation::_call_r_initialize_block_array_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_initialize_block_array((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void s_write();
 */
void CProxy_Simulation::s_write(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::s_write(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::s_write(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_s_write_void() {
  int epidx = CkRegisterEp("s_write()",
      reinterpret_cast<CkCallFnPtr>(_call_s_write_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Simulation::_call_s_write_void(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->s_write();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Simulation::s_write_4_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write(CkReductionMsg* impl_msg);
 */
void CProxy_Simulation::r_write(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_write(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_write(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_write_CkReductionMsg() {
  int epidx = CkRegisterEp("r_write(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_write_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Simulation::_call_r_write_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_write((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write_checkpoint_output();
 */
void CProxy_Simulation::r_write_checkpoint_output(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_write_checkpoint_output(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_write_checkpoint_output(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_write_checkpoint_output_void() {
  int epidx = CkRegisterEp("r_write_checkpoint_output()",
      reinterpret_cast<CkCallFnPtr>(_call_r_write_checkpoint_output_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Simulation::_call_r_write_checkpoint_output_void(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_write_checkpoint_output();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Simulation::r_write_checkpoint_output_6_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_start(CkReductionMsg* impl_msg);
 */
void CProxy_Simulation::r_restart_start(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_restart_start(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_restart_start(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_restart_start_CkReductionMsg() {
  int epidx = CkRegisterEp("r_restart_start(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_restart_start_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Simulation::_call_r_restart_start_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_restart_start((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_enter(const std::string &dir);
 */
void CProxy_Simulation::p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_restart_enter(const std::string &dir, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_restart_enter(const std::string &dir, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_restart_enter_marshall8() {
  int epidx = CkRegisterEp("p_restart_enter(const std::string &dir)",
      reinterpret_cast<CkCallFnPtr>(_call_p_restart_enter_marshall8), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_restart_enter_marshall8);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_restart_enter_marshall8);

  return epidx;
}

void CkIndex_Simulation::_call_p_restart_enter_marshall8(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_restart_enter(std::move(dir.t));
}
int CkIndex_Simulation::_callmarshall_p_restart_enter_marshall8(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_restart_enter(std::move(dir.t));
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_p_restart_enter_marshall8(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const std::string &dir*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<std::string> dir;
  implP|dir;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("dir");
  implDestP|dir;
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_restart_enter_8_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int n, const char *buffer);
 */
void CProxy_Simulation::p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int n, const char *buffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_output_write(int n, const char *buffer, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: int n, const char *buffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_output_write(int n, const char *buffer, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: int n, const char *buffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_output_write_marshall9() {
  int epidx = CkRegisterEp("p_output_write(int n, const char *buffer)",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_write_marshall9), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_output_write_marshall9);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_output_write_marshall9);

  return epidx;
}

void CkIndex_Simulation::_call_p_output_write_marshall9(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *buffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  impl_obj->p_output_write(std::move(n.t), buffer);
}
int CkIndex_Simulation::_callmarshall_p_output_write_marshall9(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int n, const char *buffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  impl_obj->p_output_write(std::move(n.t), buffer);
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_p_output_write_marshall9(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int n, const char *buffer*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> n;
  implP|n;
  int impl_off_buffer, impl_cnt_buffer;
  implP|impl_off_buffer;
  implP|impl_cnt_buffer;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  char *buffer=(char *)(impl_buf+impl_off_buffer);
  if (implDestP.hasComments()) implDestP.comment("n");
  implDestP|n;
  if (implDestP.hasComments()) implDestP.comment("buffer");
  implDestP.synchronize(PUP::sync_begin_array);
  for (int impl_i=0;impl_i*(sizeof(*buffer))<impl_cnt_buffer;impl_i++) {
    implDestP.synchronize(PUP::sync_item);
    implDestP|buffer[impl_i];
  }
  implDestP.synchronize(PUP::sync_end_array);
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_output_write_9_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
void CProxy_Simulation::r_output_barrier(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_output_barrier(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_output_barrier(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_output_barrier_CkReductionMsg() {
  int epidx = CkRegisterEp("r_output_barrier(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_output_barrier_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Simulation::_call_r_output_barrier_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_output_barrier((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_start(int index_output);
 */
void CProxy_Simulation::p_output_start(int index_output, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int index_output
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_output_start(int index_output, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: int index_output
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
  }
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_output_start(int index_output, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: int index_output
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
  }
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_output_start_marshall11() {
  int epidx = CkRegisterEp("p_output_start(int index_output)",
      reinterpret_cast<CkCallFnPtr>(_call_p_output_start_marshall11), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_output_start_marshall11);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_output_start_marshall11);

  return epidx;
}

void CkIndex_Simulation::_call_p_output_start_marshall11(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index_output*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_output_start(std::move(index_output.t));
}
int CkIndex_Simulation::_callmarshall_p_output_start_marshall11(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int index_output*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_output_start(std::move(index_output.t));
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_p_output_start_marshall11(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int index_output*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> index_output;
  implP|index_output;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("index_output");
  implDestP|index_output;
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_output_start_11_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
void CProxy_Simulation::r_monitor_performance_reduce(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::r_monitor_performance_reduce(CkReductionMsg* impl_msg, int npes, int *pes) {
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::r_monitor_performance_reduce(CkReductionMsg* impl_msg, CmiGroup &grp) {
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_r_monitor_performance_reduce_CkReductionMsg() {
  int epidx = CkRegisterEp("r_monitor_performance_reduce(CkReductionMsg* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_r_monitor_performance_reduce_CkReductionMsg), CMessage_CkReductionMsg::__idx, __idx, 0);
  CkRegisterMessagePupFn(epidx, (CkMessagePupFn)CkReductionMsg::ckDebugPup);
  return epidx;
}

void CkIndex_Simulation::_call_r_monitor_performance_reduce_CkReductionMsg(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->r_monitor_performance_reduce((CkReductionMsg*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_monitor_performance();
 */
void CProxy_Simulation::p_monitor_performance(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_monitor_performance(int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_monitor_performance(CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_monitor_performance_void() {
  int epidx = CkRegisterEp("p_monitor_performance()",
      reinterpret_cast<CkCallFnPtr>(_call_p_monitor_performance_void), 0, __idx, 0);
  return epidx;
}

void CkIndex_Simulation::_call_p_monitor_performance_void(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  impl_obj->p_monitor_performance();
  if(UsrToEnv(impl_msg)->isVarSysMsg() == 0)
    CkFreeSysMsg(impl_msg);
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_monitor_performance_13_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_block_array(const CProxy_Block &block_array);
 */
void CProxy_Simulation::p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_Block &block_array
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
  }
  if (ckIsDelegated()) {
     CkGroupMsgPrep(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID());
     ckDelegatedTo()->GroupBroadcast(ckDelegatedPtr(),CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID());
  } else CkBroadcastMsgBranch(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID(),0);
}
void CProxy_Simulation::p_set_block_array(const CProxy_Block &block_array, int npes, int *pes, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_Block &block_array
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
  }
  CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID(), npes, pes,0);
}
void CProxy_Simulation::p_set_block_array(const CProxy_Block &block_array, CmiGroup &grp, const CkEntryOptions *impl_e_opts) {
  //Marshall: const CProxy_Block &block_array
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
  }
  CkSendMsgBranchGroup(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetGroupID(), grp,0);
}

// Entry point registration function
int CkIndex_Simulation::reg_p_set_block_array_marshall14() {
  int epidx = CkRegisterEp("p_set_block_array(const CProxy_Block &block_array)",
      reinterpret_cast<CkCallFnPtr>(_call_p_set_block_array_marshall14), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_p_set_block_array_marshall14);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_p_set_block_array_marshall14);

  return epidx;
}

void CkIndex_Simulation::_call_p_set_block_array_marshall14(void* impl_msg, void* impl_obj_void)
{
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_Block &block_array*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_Block> block_array;
  implP|block_array;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_block_array(std::move(block_array.t));
}
int CkIndex_Simulation::_callmarshall_p_set_block_array_marshall14(char* impl_buf, void* impl_obj_void) {
  Simulation* impl_obj = static_cast<Simulation*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: const CProxy_Block &block_array*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_Block> block_array;
  implP|block_array;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  impl_obj->p_set_block_array(std::move(block_array.t));
  return implP.size();
}
void CkIndex_Simulation::_marshallmessagepup_p_set_block_array_marshall14(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: const CProxy_Block &block_array*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<CProxy_Block> block_array;
  implP|block_array;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("block_array");
  implDestP|block_array;
}
PUPable_def(SINGLE_ARG(Closure_Simulation::p_set_block_array_14_closure))
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_Simulation::reg_Simulation_CkMigrateMessage() {
  int epidx = CkRegisterEp("Simulation(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_Simulation_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_Simulation::_call_Simulation_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) Simulation((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(const char *filename, int n);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_get_msg_refine(const Index &impl_noname_0);
 */
void CProxySection_Simulation::p_get_msg_refine(const Index &impl_noname_0, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const Index &impl_noname_0
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<Index>::type>::type &)impl_noname_0;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_get_msg_refine_marshall2(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_initialize_block_array(CkReductionMsg* impl_msg);
 */
void CProxySection_Simulation::r_initialize_block_array(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_initialize_block_array_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void s_write();
 */
void CProxySection_Simulation::s_write(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_s_write_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_s_write_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write(CkReductionMsg* impl_msg);
 */
void CProxySection_Simulation::r_write(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_write_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_write_checkpoint_output();
 */
void CProxySection_Simulation::r_write_checkpoint_output(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_write_checkpoint_output_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_restart_start(CkReductionMsg* impl_msg);
 */
void CProxySection_Simulation::r_restart_start(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_restart_start_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_restart_enter(const std::string &dir);
 */
void CProxySection_Simulation::p_restart_enter(const std::string &dir, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const std::string &dir
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<std::string>::type>::type &)dir;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_restart_enter_marshall8(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_write(int n, const char *buffer);
 */
void CProxySection_Simulation::p_output_write(int n, const char *buffer, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int n, const char *buffer
  int impl_off=0;
  int impl_arrstart=0;
  int impl_off_buffer, impl_cnt_buffer;
  impl_off_buffer=impl_off=CK_ALIGN(impl_off,sizeof(char));
  impl_off+=(impl_cnt_buffer=sizeof(char)*(n));
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
    impl_arrstart=CK_ALIGN(implP.size(),16);
    impl_off+=impl_arrstart;
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|n;
    implP|impl_off_buffer;
    implP|impl_cnt_buffer;
  }
  char *impl_buf=impl_msg->msgBuf+impl_arrstart;
  memcpy(impl_buf+impl_off_buffer,buffer,impl_cnt_buffer);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_output_write_marshall9(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_output_barrier(CkReductionMsg* impl_msg);
 */
void CProxySection_Simulation::r_output_barrier(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_output_barrier_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_output_start(int index_output);
 */
void CProxySection_Simulation::p_output_start(int index_output, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: int index_output
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|index_output;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|index_output;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_output_start_marshall11(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
 */
void CProxySection_Simulation::r_monitor_performance_reduce(CkReductionMsg* impl_msg)
{
  ckCheck();
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_r_monitor_performance_reduce_CkReductionMsg(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_monitor_performance();
 */
void CProxySection_Simulation::p_monitor_performance(const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  void *impl_msg = CkAllocSysMsg(impl_e_opts);
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_monitor_performance_void(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: void p_set_block_array(const CProxy_Block &block_array);
 */
void CProxySection_Simulation::p_set_block_array(const CProxy_Block &block_array, const CkEntryOptions *impl_e_opts)
{
  ckCheck();
  //Marshall: const CProxy_Block &block_array
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    //Have to cast away const-ness to get pup routine
    implP|(typename std::remove_cv<typename std::remove_reference<CProxy_Block>::type>::type &)block_array;
  }
  if (ckIsDelegated()) {
     ckDelegatedTo()->GroupSectionSend(ckDelegatedPtr(),CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg, ckGetNumSections(), ckGetSectionIDs());
  } else {
    void *impl_msg_tmp;
    for (int i=0; i<ckGetNumSections(); ++i) {
       impl_msg_tmp= (i<ckGetNumSections()-1) ? CkCopyMsg((void **) &impl_msg):impl_msg;
       CkSendMsgBranchMulti(CkIndex_Simulation::idx_p_set_block_array_marshall14(), impl_msg_tmp, ckGetGroupIDn(i), ckGetNumElements(i), ckGetElements(i),0);
    }
  }
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: Simulation(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_Simulation::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_IrrGroup::__idx);
   CkRegisterGroupIrr(__idx,Simulation::isIrreducible());
  // REG: Simulation(const char *filename, int n);
  idx_Simulation_marshall1();

  // REG: void p_get_msg_refine(const Index &impl_noname_0);
  idx_p_get_msg_refine_marshall2();

  // REG: void r_initialize_block_array(CkReductionMsg* impl_msg);
  idx_r_initialize_block_array_CkReductionMsg();

  // REG: void s_write();
  idx_s_write_void();

  // REG: void r_write(CkReductionMsg* impl_msg);
  idx_r_write_CkReductionMsg();

  // REG: void r_write_checkpoint_output();
  idx_r_write_checkpoint_output_void();

  // REG: void r_restart_start(CkReductionMsg* impl_msg);
  idx_r_restart_start_CkReductionMsg();

  // REG: void p_restart_enter(const std::string &dir);
  idx_p_restart_enter_marshall8();

  // REG: void p_output_write(int n, const char *buffer);
  idx_p_output_write_marshall9();

  // REG: void r_output_barrier(CkReductionMsg* impl_msg);
  idx_r_output_barrier_CkReductionMsg();

  // REG: void p_output_start(int index_output);
  idx_p_output_start_marshall11();

  // REG: void r_monitor_performance_reduce(CkReductionMsg* impl_msg);
  idx_r_monitor_performance_reduce_CkReductionMsg();

  // REG: void p_monitor_performance();
  idx_p_monitor_performance_void();

  // REG: void p_set_block_array(const CProxy_Block &block_array);
  idx_p_set_block_array_marshall14();

  // REG: Simulation(CkMigrateMessage* impl_msg);
  idx_Simulation_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_Simulation_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: group MappingArray: CkArrayMap{
MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
MappingArray(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_MappingArray::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
CkGroupID CProxy_MappingArray::ckNew(int impl_noname_1, int impl_noname_2, int impl_noname_3, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_1, int impl_noname_2, int impl_noname_3
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_1;
    implP|impl_noname_2;
    implP|impl_noname_3;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_1;
    implP|impl_noname_2;
    implP|impl_noname_3;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  CkGroupID gId = CkCreateGroup(CkIndex_MappingArray::__idx, CkIndex_MappingArray::idx_MappingArray_marshall1(), impl_msg);
  return gId;
}
  CProxy_MappingArray::CProxy_MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_1, int impl_noname_2, int impl_noname_3
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_1;
    implP|impl_noname_2;
    implP|impl_noname_3;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_1;
    implP|impl_noname_2;
    implP|impl_noname_3;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  ckSetGroupID(CkCreateGroup(CkIndex_MappingArray::__idx, CkIndex_MappingArray::idx_MappingArray_marshall1(), impl_msg));
}

// Entry point registration function
int CkIndex_MappingArray::reg_MappingArray_marshall1() {
  int epidx = CkRegisterEp("MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingArray_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_MappingArray_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_MappingArray_marshall1);

  return epidx;
}

void CkIndex_MappingArray::_call_MappingArray_marshall1(void* impl_msg, void* impl_obj_void)
{
  MappingArray* impl_obj = static_cast<MappingArray*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_1, int impl_noname_2, int impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_3;
  implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingArray(std::move(impl_noname_1.t), std::move(impl_noname_2.t), std::move(impl_noname_3.t));
}
int CkIndex_MappingArray::_callmarshall_MappingArray_marshall1(char* impl_buf, void* impl_obj_void) {
  MappingArray* impl_obj = static_cast<MappingArray*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int impl_noname_1, int impl_noname_2, int impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_3;
  implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingArray(std::move(impl_noname_1.t), std::move(impl_noname_2.t), std::move(impl_noname_3.t));
  return implP.size();
}
void CkIndex_MappingArray::_marshallmessagepup_MappingArray_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_1, int impl_noname_2, int impl_noname_3*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_1;
  implP|impl_noname_1;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_2;
  implP|impl_noname_2;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_3;
  implP|impl_noname_3;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_1");
  implDestP|impl_noname_1;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_2");
  implDestP|impl_noname_2;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_3");
  implDestP|impl_noname_3;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_MappingArray::reg_MappingArray_CkMigrateMessage() {
  int epidx = CkRegisterEp("MappingArray(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingArray_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_MappingArray::_call_MappingArray_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) MappingArray((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingArray(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_MappingArray::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_CkArrayMap::__idx);
   CkRegisterGroupIrr(__idx,MappingArray::isIrreducible());
  // REG: MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
  idx_MappingArray_marshall1();

  // REG: MappingArray(CkMigrateMessage* impl_msg);
  idx_MappingArray_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_MappingArray_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: group MappingTree: CkArrayMap{
MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
MappingTree(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_MappingTree::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
CkGroupID CProxy_MappingTree::ckNew(int impl_noname_4, int impl_noname_5, int impl_noname_6, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_4, int impl_noname_5, int impl_noname_6
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_4;
    implP|impl_noname_5;
    implP|impl_noname_6;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_4;
    implP|impl_noname_5;
    implP|impl_noname_6;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  CkGroupID gId = CkCreateGroup(CkIndex_MappingTree::__idx, CkIndex_MappingTree::idx_MappingTree_marshall1(), impl_msg);
  return gId;
}
  CProxy_MappingTree::CProxy_MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_4, int impl_noname_5, int impl_noname_6
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_4;
    implP|impl_noname_5;
    implP|impl_noname_6;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_4;
    implP|impl_noname_5;
    implP|impl_noname_6;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  ckSetGroupID(CkCreateGroup(CkIndex_MappingTree::__idx, CkIndex_MappingTree::idx_MappingTree_marshall1(), impl_msg));
}

// Entry point registration function
int CkIndex_MappingTree::reg_MappingTree_marshall1() {
  int epidx = CkRegisterEp("MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingTree_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_MappingTree_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_MappingTree_marshall1);

  return epidx;
}

void CkIndex_MappingTree::_call_MappingTree_marshall1(void* impl_msg, void* impl_obj_void)
{
  MappingTree* impl_obj = static_cast<MappingTree*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_4, int impl_noname_5, int impl_noname_6*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_4;
  implP|impl_noname_4;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_5;
  implP|impl_noname_5;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_6;
  implP|impl_noname_6;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingTree(std::move(impl_noname_4.t), std::move(impl_noname_5.t), std::move(impl_noname_6.t));
}
int CkIndex_MappingTree::_callmarshall_MappingTree_marshall1(char* impl_buf, void* impl_obj_void) {
  MappingTree* impl_obj = static_cast<MappingTree*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int impl_noname_4, int impl_noname_5, int impl_noname_6*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_4;
  implP|impl_noname_4;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_5;
  implP|impl_noname_5;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_6;
  implP|impl_noname_6;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingTree(std::move(impl_noname_4.t), std::move(impl_noname_5.t), std::move(impl_noname_6.t));
  return implP.size();
}
void CkIndex_MappingTree::_marshallmessagepup_MappingTree_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_4, int impl_noname_5, int impl_noname_6*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_4;
  implP|impl_noname_4;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_5;
  implP|impl_noname_5;
  PUP::detail::TemporaryObjectHolder<int> impl_noname_6;
  implP|impl_noname_6;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_4");
  implDestP|impl_noname_4;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_5");
  implDestP|impl_noname_5;
  if (implDestP.hasComments()) implDestP.comment("impl_noname_6");
  implDestP|impl_noname_6;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_MappingTree::reg_MappingTree_CkMigrateMessage() {
  int epidx = CkRegisterEp("MappingTree(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingTree_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_MappingTree::_call_MappingTree_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) MappingTree((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingTree(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_MappingTree::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_CkArrayMap::__idx);
   CkRegisterGroupIrr(__idx,MappingTree::isIrreducible());
  // REG: MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
  idx_MappingTree_marshall1();

  // REG: MappingTree(CkMigrateMessage* impl_msg);
  idx_MappingTree_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_MappingTree_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

/* DEFS: group MappingIo: CkArrayMap{
MappingIo(int impl_noname_7);
MappingIo(CkMigrateMessage* impl_msg);
};
 */
#ifndef CK_TEMPLATES_ONLY
 int CkIndex_MappingIo::__idx=0;
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(int impl_noname_7);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(int impl_noname_7);
 */
CkGroupID CProxy_MappingIo::ckNew(int impl_noname_7, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_7
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_7;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_7;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  CkGroupID gId = CkCreateGroup(CkIndex_MappingIo::__idx, CkIndex_MappingIo::idx_MappingIo_marshall1(), impl_msg);
  return gId;
}
  CProxy_MappingIo::CProxy_MappingIo(int impl_noname_7, const CkEntryOptions *impl_e_opts)
{
  //Marshall: int impl_noname_7
  int impl_off=0;
  { //Find the size of the PUP'd data
    PUP::sizer implP;
    implP|impl_noname_7;
    impl_off+=implP.size();
  }
  CkMarshallMsg *impl_msg=CkAllocateMarshallMsg(impl_off,impl_e_opts);
  { //Copy over the PUP'd data
    PUP::toMem implP((void *)impl_msg->msgBuf);
    implP|impl_noname_7;
  }
  UsrToEnv(impl_msg)->setMsgtype(BocInitMsg);
  ckSetGroupID(CkCreateGroup(CkIndex_MappingIo::__idx, CkIndex_MappingIo::idx_MappingIo_marshall1(), impl_msg));
}

// Entry point registration function
int CkIndex_MappingIo::reg_MappingIo_marshall1() {
  int epidx = CkRegisterEp("MappingIo(int impl_noname_7)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingIo_marshall1), CkMarshallMsg::__idx, __idx, 0+CK_EP_NOKEEP);
  CkRegisterMarshallUnpackFn(epidx, _callmarshall_MappingIo_marshall1);
  CkRegisterMessagePupFn(epidx, _marshallmessagepup_MappingIo_marshall1);

  return epidx;
}

void CkIndex_MappingIo::_call_MappingIo_marshall1(void* impl_msg, void* impl_obj_void)
{
  MappingIo* impl_obj = static_cast<MappingIo*>(impl_obj_void);
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_7*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_7;
  implP|impl_noname_7;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingIo(std::move(impl_noname_7.t));
}
int CkIndex_MappingIo::_callmarshall_MappingIo_marshall1(char* impl_buf, void* impl_obj_void) {
  MappingIo* impl_obj = static_cast<MappingIo*>(impl_obj_void);
  envelope *env = UsrToEnv(impl_buf);
  /*Unmarshall pup'd fields: int impl_noname_7*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_7;
  implP|impl_noname_7;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  new (impl_obj_void) MappingIo(std::move(impl_noname_7.t));
  return implP.size();
}
void CkIndex_MappingIo::_marshallmessagepup_MappingIo_marshall1(PUP::er &implDestP,void *impl_msg) {
  CkMarshallMsg *impl_msg_typed=(CkMarshallMsg *)impl_msg;
  char *impl_buf=impl_msg_typed->msgBuf;
  envelope *env = UsrToEnv(impl_msg_typed);
  /*Unmarshall pup'd fields: int impl_noname_7*/
  PUP::fromMem implP(impl_buf);
  PUP::detail::TemporaryObjectHolder<int> impl_noname_7;
  implP|impl_noname_7;
  impl_buf+=CK_ALIGN(implP.size(),16);
  /*Unmarshall arrays:*/
  if (implDestP.hasComments()) implDestP.comment("impl_noname_7");
  implDestP|impl_noname_7;
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(CkMigrateMessage* impl_msg);
 */

// Entry point registration function
int CkIndex_MappingIo::reg_MappingIo_CkMigrateMessage() {
  int epidx = CkRegisterEp("MappingIo(CkMigrateMessage* impl_msg)",
      reinterpret_cast<CkCallFnPtr>(_call_MappingIo_CkMigrateMessage), 0, __idx, 0);
  return epidx;
}

void CkIndex_MappingIo::_call_MappingIo_CkMigrateMessage(void* impl_msg, void* impl_obj_void)
{
  new (impl_obj_void) MappingIo((CkMigrateMessage*)impl_msg);
}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(int impl_noname_7);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
/* DEFS: MappingIo(CkMigrateMessage* impl_msg);
 */
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
void CkIndex_MappingIo::__register(const char *s, size_t size) {
  __idx = CkRegisterChare(s, size, TypeGroup);
  CkRegisterBase(__idx, CkIndex_CkArrayMap::__idx);
   CkRegisterGroupIrr(__idx,MappingIo::isIrreducible());
  // REG: MappingIo(int impl_noname_7);
  idx_MappingIo_marshall1();

  // REG: MappingIo(CkMigrateMessage* impl_msg);
  idx_MappingIo_CkMigrateMessage();
  CkRegisterMigCtor(__idx, idx_MappingIo_CkMigrateMessage());

}
#endif /* CK_TEMPLATES_ONLY */

#ifndef CK_TEMPLATES_ONLY
void _registersimulation(void)
{
  static int _done = 0; if(_done) return; _done = 1;
  _registermesh();

  CkRegisterReadonly("proxy_simulation","CProxy_Simulation",sizeof(proxy_simulation),(void *) &proxy_simulation,__xlater_roPup_proxy_simulation);

  _registerInitCall(method_close_files_mutex_init,1);

/* REG: group Simulation: IrrGroup{
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
  CkIndex_Simulation::__register("Simulation", sizeof(Simulation));

/* REG: group MappingArray: CkArrayMap{
MappingArray(int impl_noname_1, int impl_noname_2, int impl_noname_3);
MappingArray(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_MappingArray::__register("MappingArray", sizeof(MappingArray));

/* REG: group MappingTree: CkArrayMap{
MappingTree(int impl_noname_4, int impl_noname_5, int impl_noname_6);
MappingTree(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_MappingTree::__register("MappingTree", sizeof(MappingTree));

/* REG: group MappingIo: CkArrayMap{
MappingIo(int impl_noname_7);
MappingIo(CkMigrateMessage* impl_msg);
};
*/
  CkIndex_MappingIo::__register("MappingIo", sizeof(MappingIo));

}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_Simulation::virtual_pup(PUP::er &p) {
    recursive_pup<Simulation>(dynamic_cast<Simulation*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_MappingArray::virtual_pup(PUP::er &p) {
    recursive_pup<MappingArray>(dynamic_cast<MappingArray*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_MappingTree::virtual_pup(PUP::er &p) {
    recursive_pup<MappingTree>(dynamic_cast<MappingTree*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
#ifndef CK_TEMPLATES_ONLY
template <>
void CBase_MappingIo::virtual_pup(PUP::er &p) {
    recursive_pup<MappingIo>(dynamic_cast<MappingIo*>(this), p);
}
#endif /* CK_TEMPLATES_ONLY */
